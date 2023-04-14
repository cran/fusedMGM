require(fastDummies) ;
require(Rcpp) ;
require(parallel) ;
require(bigmemory) ;
require(bigalgebra) ;
require(biganalytics) ;

#' Defining S3 object "MGM"
#
#' @param X data frame or matrix of continuous variables (row: observation, column: variable)
#' @param Y data frame or matrix of discrete variables (row: observation, column: variable)
#' @param g group index, needed for temporary files
#' @return An S3 `MGM` object, containing data, network parameters, and the 1st derivatives
#' @import bigmemory
#' @export
MGM	<- function(X, Y, g) {
	input_missing	<- missing(X) | missing(Y) ;
	
	if (input_missing) {
		stop("Please provide input data")
	}
	
	if (!is.null(X)) {
		tryCatch(X	<- as.matrix(X),
			error = function(e) stop("Continuous data should be able to be converted to numeric matrix")) ;
	}
	if (!is.null(Y)) {
		tryCatch(Y	<- as.matrix(Y),
			error = function(e) stop("Discrete data should be able to be converted to character matrix")) ;
	}
	
	if (!(typeof(X) %in% c("double"))) stop("1st input data should be numeric") ;
	if (!all(apply(Y, 2, function(t) all(is.factor(t) | is.character(t))))) {
		tryCatch(Y 	<- apply(Y, 2, as.character), error = function(e) stop("2nd input data should be characters or factors")) ;
	}
	
	n_x	<- nrow(X) ; p_x	<- ncol(X) ;
	n_y	<- nrow(Y) ; p_y	<- ncol(Y) ;
	
	if (n_x != n_y) stop("Numbers of samples in continuous/discrete data are different") ;
	
	X		<- as.big.matrix(X) ;
	
	# Convert discrete data into a dummy matrix
	# Define list & number of unique levels of discrete variables
    L_y		<- numeric(p_y) ;
	Y_dummy	<- fastDummies::dummy_cols(Y, remove_selected_columns=TRUE) ;
	
	levels_y	<- list() ;
	L_cumsum	<- numeric(p_y + 1) ;
	L_cumsum[1]	<- 0 ;
    for (r in seq(p_y)) {
		varname			<- colnames(Y)[r] ;
		L_y[r]			<- length(unique(Y[,r])) ;
		L_cumsum[r+1]	<- L_cumsum[r] + L_y[r] ;
        levels_y[[r]]	<- gsub(paste0(varname, "_"), "", colnames(Y_dummy)[(L_cumsum[r]+1):(L_cumsum[r+1])], fixed=TRUE) ;
	}
	L_sum	<- sum(L_y) ;
	Y_dummy	<- as.matrix(Y_dummy) ;
	mode(Y_dummy)	<- "double" ;
	Y_dummy	<- as.big.matrix(Y_dummy) ;
	rownames(Y_dummy)	<- rownames(Y) ;
    
	# Convert discrete data into level indices
    Y_lev	<- matrix(0, nrow=n_y, ncol=p_y) ;
    for (n in seq(n_y)) { for (r in seq(p_y)) {
		Y_lev[n,r]	<- match(as.character(Y[n,r]), levels_y[[r]]) ;
	}}
	Y_lev	<- as.big.matrix(Y_lev) ;
	colnames(Y_lev)	<- colnames(Y) ;
	rownames(Y_lev)	<- rownames(Y) ;
	
	gc() ;
	
	# Initialize parameters
	alpha	<- big.matrix(nrow=p_x, ncol=1, init=0) ;
	rownames(alpha)	<- colnames(X) ;
	
	# Cont-cont 
	beta	<- big.matrix(nrow=p_x, ncol=p_x, init=0) ;
	for (r in seq(p_x)) {
		beta[r,r]	<- 1 ;
	}
	rownames(beta)	<- colnames(beta)	<- colnames(X) ;
	
	# Disc-cont
	rho		<- big.matrix(nrow=L_sum, ncol=p_x, init=0) ;
	rownames(rho)	<- colnames(Y_dummy) ;
	colnames(rho)	<- colnames(X) ;
	
	# Cont-cont
	phi		<- big.matrix(nrow=L_sum, ncol=L_sum, init=0) ;
	for (l in seq(L_sum)) {
		phi[l,l]	<- 1 ;
	}
	rownames(phi)	<- colnames(phi)	<- colnames(Y_dummy) ;
	
	# Define linear loss function / predictor
	L_x		<- big.matrix(nrow=n_x, ncol=p_x, init=0) ;
	b_y		<- big.matrix(nrow=n_y, ncol=L_sum, init=0) ;
	Prob_y	<- big.matrix(nrow=n_y, ncol=L_sum, init=0) ;
	
	# Define the slopes
	d_alpha	<- deepcopy(alpha) ;
	d_beta	<- deepcopy(beta) ;
	d_rho	<- deepcopy(rho) ;
	d_phi	<- deepcopy(phi) ;
	
	MGMobj	<- list(Continuous=X, Discrete=Y_lev, Dummy=Y_dummy, Categories=levels_y, Levels=L_y,
					alpha=alpha, beta=beta, rho=rho, phi=phi,
					Cont_loss=L_x, Disc_pred=b_y, Disc_prob=Prob_y,
					d_alpha=d_alpha, d_beta=d_beta, d_rho=d_rho, d_phi=d_phi) ;
	attr(MGMobj, "class")	<- "MGM" ;
	
	return(MGMobj)
}

#' Make MGM lists from input data
#
#' @param X data frame or matrix of continuous variables (row: observation, column: variable)
#' @param Y data frame or matrix of discrete variables (row: observation, column: variable)
#' @param group group variable vector, with the sample names
#' @return A list of MGM objects. The length is equal to the unique number of groups.
#' @export
make_MGM_list	<- function(X, Y, group) {
							
	if (is.null(rownames(X)) | is.null(rownames(Y)) | is.null(names(group))) {
		stop("Please provide sample names for data or group variable!") ;
	}
	
	if (!is.na(suppressWarnings(as.numeric(rownames(X)[1])))) {
		warning("Seems the sample names are literally row names e.g 1, 2, ... \n please check again to prevent the possibly erroneous results") ;
	}
							
	sample_common	<- intersect(rownames(X), rownames(Y)) ;
	sample_common	<- intersect(sample_common, names(group)) ;
	
	ind_y	<- seq(ncol(Y)) ;
	
	if (length(sample_common) == 0) stop("No common samples, please check the sample names") ;
	
	X		<- X[sample_common,,drop=FALSE] ;
	Y		<- Y[sample_common,,drop=FALSE] ;
	group	<- group[sample_common] ;
	
	group_uniq	<- unique(group) ;
	G	<- length(group_uniq) ;
	
	# Split data by group
	X_list	<- rep(list(NULL), G) ;
	Y_list	<- rep(list(NULL), G) ;
	for (g in seq(G)) {
		X_list[[g]]	<- X[which(group[rownames(X)]==group_uniq[g]),,drop=FALSE] ;
		Y_list[[g]]	<- Y[which(group[rownames(Y)]==group_uniq[g]),,drop=FALSE] ;
	}
	names(X_list)	<- group_uniq ;
	names(Y_list)	<- group_uniq ;
	
	# For discrete variables, preserve data only with common levels
	# If there are only 1 or less common level, drop the variable
	
	drop_y	<- c() ;
	drop_sample	<- rep(list(c()), G) ;
	
	while (TRUE) {
	
		for (r in ind_y) {
			# Note: levels with singleton will be ignored
			lev_uniq	<- Reduce(intersect, lapply(Y_list, function(y) {
																tab_uniq	<- table(y[,r]) ;
																levlist		<- names(tab_uniq)[tab_uniq > 1] ;
																return(levlist) ;
																})) ;
			
			if (length(lev_uniq) <= 1) {	# 1 or less common level: drop the variable
				drop_y	<- c(drop_y, r) ;
			} else {						# Else, drop samples not in the common levels
			
				for (g in seq(G)) {
					y	<- Y_list[[g]] ;
					drop_sample[[g]]	<- c(drop_sample[[g]], rownames(y)[!(y[,r] %in% lev_uniq)]) ;
				}
				
			}
		}
		
		# Drop variables / samples
		for (g in seq(G)) {
			drop_sample[[g]]	<- unique(drop_sample[[g]]) ;
		
			if (length(drop_sample[[g]]) > 0) {
				X_list[[g]]	<- X_list[[g]][which(!(rownames(X_list[[g]]) %in% drop_sample[[g]])),,drop=FALSE] ;
				Y_list[[g]]	<- Y_list[[g]][which(!(rownames(Y_list[[g]]) %in% drop_sample[[g]])),,drop=FALSE] ;
			}
			if (length(drop_y) > 0) {
				Y_list[[g]]	<- Y_list[[g]][, -drop_y, drop=FALSE] ;
			}
			
			if (nrow(Y_list[[g]]) == 0) {
				stop(paste0("No sample left in group ", g, "!")) ;
			}
		}
		
		# Count number of dropped samples
		n_dropped	<- sum(sapply(drop_sample, length)) ;
		
		if ((n_dropped == 0) | (ncol(Y_list[[1]]) == 0)) {	# No samples dropped or no discrete data left: end dropping
			break ;
		} else {				# Otherwise, update 
			ind_y		<- if (!is.null(drop_y)) seq_along(ind_y[-drop_y]) else ind_y ;
			drop_sample	<- rep(list(c()), G) ;
		}
	}
	
	gc() ;
	
	# For continuous variables, center & scale
	for (g in seq(G)) {
		X_list[[g]]		<- scale(X_list[[g]], center=TRUE, scale=FALSE) ;
		X_list[[g]]		<- scale(X_list[[g]], scale=apply(X_list[[g]], 2, function(x) sqrt(sum(x^2/nrow(X_list[[g]]))))) ;
	}
	
	# Make MGM for each group
	MGM_list	<- list() ;
	for (g in seq(G)) {
		MGM_list[[g]]	<- MGM(X_list[[g]], Y_list[[g]], g) ;
	}
	
	return(MGM_list) ;
}

#' Main function of fused MGM
#'
#' Infers networks from 2-class mixed data
#'
#' If the value of Lipschitz constant, L, is not provided, the backtracking will be performed
#' @param data Data frame with rows as observations and columns as variables
#' @param ind_disc Indices of discrete variables
#' @param group Group indices, must be provided with the observation names
#' @param t Numeric. Initial value of coefficient that reflect 2 previous iterations in fast proximal gradient method. Default: 1
#' @param L Numeric. Initial guess of Lipschitz constant. Default: missing (use backtracking)
#' @param eta Numeric. Multipliers for L in backtracking. Default: 2
#' @param lambda_intra Vector with 3 numeric variables. Penalization parameters for network edge weights
#' @param lambda_intra_prior Vector with 3 numeric variables. Penalization parameters for network edge weights, applied to the edges with prior information
#' @param lambda_inter Vector with 3 numeric variables. Penalization parameters for network edge weight differences
#' @param with_prior Logical. Is prior information provided? Default: FALSE
#' @param prior_list List of prior information. Each element must be a 3-column data frames, with the 1st and the 2nd columns being variable names and the 3rd column being prior confidence (0,1)
#' @param converge_by_edge Logical. The convergence should be judged by null differences of network edges after iteration. If FALSE, the rooted mean square difference (RMSD) of edge weights is used. Default: TRUE
#' @param tol_edge Integer. Number of consecutive iterations of convergence to stop the iteration. Default: 3
#' @param tol_mgm Numeric. Cutoff of network edge RMSD for convergence. Default: 1e-04
#' @param tol_g Numeric. Cutoff of iternations in prox-grad map calculation. Default: 5e-03
#' @param tol_fpa Numeric. Cutoff for fixed-point approach. Default: 1e-12
#' @param maxit Integer. Maximum number of iterations in fixed-point approach. Default: 1000000
#' @param polish Logical. Should the edges with the weights below the cutoff should be discarded? Default: TRUE
#' @param tol_polish Numeric. Cutoff of polishing the resulting network. Default: 1e-12
#' @param cores Integer. Number of cores to use multi-core utilization. Default: maximum number of available cores
#' @param verbose Logical. If TRUE, the procedures are reported in real-time manner. Default: FALSE
#' @return The resulting networks, in the form of a list of MGMs
#' @import bigmemory
#' @importFrom utils combn
#' @importFrom stats rnorm
#' @importFrom stats pnorm
#' @importFrom parallel mclapply
#' @examples
#' \donttest{
#' data(data_all) ;  # Example 500-by-100 simulation data
#' data(ind_disc) ;
#' 
#' group <- rep(c(1,2), each=250) ;
#' names(group) <- seq(500) ;
#' 
#' if (Sys.info()['sysname'] == 'Windows') {
#'   cores=1
#' } else {
#'   cores=parallel::detectCores() ;
#' }
#' 
#' res_FMGM <- FMGM_mc(data_all, ind_disc, group, 
#'                     lambda_intra=c(0.2,0.15,0.1), lambda_inter=c(0.2,0.15,0.1), 
#'                     cores=cores, verbose=TRUE)
#' }
#' @export
FMGM_mc	<- function(data, ind_disc, group, t=1, L=NULL, eta=2, lambda_intra, lambda_intra_prior=NULL, lambda_inter, 
					with_prior=FALSE, prior_list=NULL, converge_by_edge=TRUE, tol_edge=3,
					tol_mgm=1e-04, tol_g=5e-03, tol_fpa=1e-12, maxit=1000000, 
					polish=TRUE, tol_polish=1e-12, 
					cores=parallel::detectCores(), verbose=FALSE){

	options(bigmemory.allow.dimnames=TRUE) ;

	if (!polish) tol_polish <- 0 ;
	
	sample_common	<- intersect(rownames(data), names(group)) ;
	
	X	<- data[,-ind_disc,drop=FALSE] ;
	Y	<- data[, ind_disc,drop=FALSE] ;
	
	MGM_list	<- make_MGM_list(X, Y, group) ;
	
	X	<- rbind(MGM_list[[1]]$Continuous[], MGM_list[[2]]$Continuous[]) ;
	Y	<- rbind(MGM_list[[1]]$Discrete[],   MGM_list[[2]]$Discrete[]) ;
	p_x	<- ncol(X) ;	p_y	<- ncol(Y) ;	pq	<- p_x + p_y ;
	rownames(Y)	<- rownames(X) ;
	data	<- cbind(X,Y) ;
	
	if (with_prior) {
		
		for (k in seq_along(prior_list)) {
			match_1 <- !is.na(match(prior_list[[k]][[1]], colnames(data))) ;
			match_2 <- !is.na(match(prior_list[[k]][[2]], colnames(data))) ;
			
			prior_list[[k]] <- prior_list[[k]][match_1 & match_2,] ;
                }
		
		# Make information array
		info_tmp	<- make_net_tmp(data, length(prior_list)) ;
		
		# Put information to the big.matrix
		for (k in seq_along(prior_list)) {
			prior_info	<- prior_list[[k]] ;
			for (j in seq(nrow(prior_info))) {
				ind_1	<- match(prior_list[[k]][j,1], colnames(data)) ;
				ind_2	<- match(prior_list[[k]][j,2], colnames(data)) ;
				info_tmp[[k]][ind_1,ind_2]	<- info_tmp[[k]][ind_2,ind_1]	<- prior_list[[k]][j,3] ;
			}
		}
		
		# Make information of edges with prior
		wp	<- big.matrix(nrow=ncol(data), ncol=ncol(data), init=0) ;
		rownames(wp)	<- colnames(wp)	<- colnames(data) ;
		data_ind	<- combn(ncol(data), 2) ;
		mclapply(seq(ncol(data_ind)), function(j) {
			ind_1	<- data_ind[1,j] ;
			ind_2	<- data_ind[2,j] ;
			
			pr_get	<- sapply(info_tmp, function(mat) return(!is.na(mat[ind_1,ind_2]))) ;
			if (any(pr_get))	wp[ind_1,ind_2]	<- wp[ind_2,ind_1]	<- 1 ;
		}, mc.cores=cores) ;
	} else {
	  wp = NULL ;
	}
	
	learn_res	<- learn_tb(MGM_list, t, L, eta, lambda_intra, lambda_intra_prior, lambda_inter,
						with_prior, wp, converge_by_edge, tol_edge,
						tol_mgm, tol_g, tol_fpa, tol_polish, maxit, cores, verbose) ;
	
	if (polish) {
		if (verbose) print("Dropping weak interactions under the given threshold") ;
		polish_MGM_list(learn_res, tol_polish, cores) ;
	}
	
	return(learn_res) ;
}
