require(fastDummies) ;
require(Rcpp) ;
require(parallel) ;
require(bigmemory) ;
require(bigalgebra) ;
require(biganalytics) ;

#' StEPS: train subsamples and calculate edge instabilities
#'
#' From large to small values of candidates, calculate the edge inference instabilities from subsamples
#' The smallest values with the instabilities under the cutoff are chosen.
#' See Sedgewich et al. (2016) for more details
#'
#' @param data Data frame with rows as observations and columns as variables
#' @param ind_disc Indices of discrete variables
#' @param group Group indices, must be provided with the observation names
#' @param lambda_list Vector with numeric variables. Penalization parameter candidates
#' @param with_prior Logical. Is prior information provided? Default: FALSE
#' @param prior_list List of prior information. Each element must be a 3-column data frames, with the 1st and the 2nd columns being variable names and the 3rd column being prior confidence (0,1)
#' @param N Integer. Number of subsamples to use. Default: 20
#' @param b Integer. Number of observations in each subsample. Default: ceiling(10*sqrt(number of total observations))
#' @param gamma Numeric. Instability cutoff. Default: 0.05
#' @param perm Integer. Number of permutations to normalize the prior confidence. Default: 10000
#' @param eps Numeric. Pseudocount to calculate the likelihood of edge detection. Default: 0.05
#' @param tol_polish Numeric. Cutoff of polishing the resulting network. Default: 1e-12
#' @param ... Other arguments sent to fast proximal gradient method
#' @param cores Integer. Number of cores to use multi-core utilization. Default: maximum number of available cores
#' @param verbose Logical. If TRUE, the procedures are reported in real-time manner. Default: FALSE
#' @return The resulting networks, in the form of a list of MGMs
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
#' lambda_list <- 2^seq(log2(.08), log2(.32), length.out=7) ;
#' lambda_list <- sort(lambda_list, decreasing=TRUE) ;
#' 
#' res_steps <- FMGM_StEPS(data_all, ind_disc, group, 
#'                     lambda_list=lambda_list, 
#'                     cores=cores, verbose=TRUE)
#' }
#' @export
FMGM_StEPS	<- function(data, ind_disc, group, lambda_list, with_prior=FALSE, prior_list=NULL, 
						N=20, b=NULL, gamma=.05, perm=10000, eps=.05, tol_polish=1e-12, ..., cores=parallel::detectCores(), verbose=FALSE) {

	if (is.null(rownames(data)) | is.null(names(group))) {
		stop("Please provide sample names for data or group variable!") ;
	}
	
	if (!is.na(suppressWarnings(as.numeric(rownames(data)[1])))) {
		warning("Seems the sample names are literally row names e.g 1, 2, ... \n please check again to prevent the possibly erroneous results") ;
	}
	
	if (with_prior & (is.null(prior_list) | length(prior_list) == 0)) {
		stop("No prior provided") ;
	}

	options(bigmemory.allow.dimnames=TRUE) ;
	
	X	<- data[,-ind_disc,drop=FALSE] ;
	Y	<- data[, ind_disc,drop=FALSE] ;
	
	orig_list	<- make_MGM_list(X, Y, group) ;
			
	G	<- length(orig_list) ;
	inst_temp	<- make_inst_tmp(orig_list) ;
	
	X_use	<- if (G == 2) rbind(orig_list[[1]]$Continuous[], orig_list[[2]]$Continuous[]) else orig_list[[1]]$Continuous[] ;
	Y_use	<- if (G == 2) rbind(orig_list[[1]]$Discrete[],   orig_list[[2]]$Discrete[])   else orig_list[[1]]$Discrete[]  ;
	p_x	<- ncol(X_use) ;	p_y	<- ncol(Y_use) ;	pq	<- p_x + p_y ;
	rownames(Y_use)	<- rownames(X_use) ;
	data_use	<- cbind(X_use,Y_use) ;
	
	sample_common	<- intersect(rownames(data_use), names(group)) ;
	n	<- length(sample_common) ;
	if (!is.null(b) && (b >= length(sample_common))) stop("Provided size of the subsamples provided should be smaller than the number of total samples") ;
	
	X	<- X[sample_common,] ;
	Y	<- Y[sample_common,] ;
	
	if (is.null(b)) {
		b	<- min(ceiling(10*sqrt(n)), ceiling(n*.8)) ;
	}
	
	if (any(lambda_list < 0)) stop("The penalty parameters should be non-negative") ;
	lambda_list	<- sort(lambda_list, decreasing=TRUE) ;
	
	# Make subsamples
	sub_index	<- matrix(0, nrow=N, ncol=b) ;
	
	if (N <= choose(length(sample_common), b)) {
		sub_cnt	<- 0 ;
		while (sub_cnt < N) {
			subsample	<- sample(sample_common, b, replace=FALSE) ;
			sub_list	<- suppressWarnings(make_MGM_list(X[subsample,,drop=FALSE], Y[subsample,,drop=FALSE], group[subsample])) ;
			
			if (length(sub_list) != length(orig_list)) { next ;
			} else if (!isTRUE(all.equal(sub_list[[1]]$Categories, orig_list[[1]]$Categories))) { next ; }
			
			sub_cnt	<- sub_cnt + 1 ;
			rm(sub_list) ; gc() ;
			sub_index[sub_cnt,]	<- subsample ;
		}
	} else {
		sub_index	<- t(combn(sample_common, b)) ;
		exc_ind	<- c() ;
		
		for (n in seq(nrow(sub_index))) {
			subsample	<- sub_index[n,] ;
			sub_list	<- suppressWarnings(make_MGM_list(X[subsample,,drop=FALSE], Y[subsample,,drop=FALSE], group[subsample])) ;
			
			if (length(sub_list) != length(orig_list)) { exc_ind	<- c(exc_ind, n) ;
			} else if (!isTRUE(all.equal(sub_list[[1]]$Categories, orig_list[[1]]$Categories))) { exc_ind	<- c(exc_ind, n) ; }
			
			rm(sub_list) ; gc() ;
		}
		
		sub_index	<- sub_index[-exc_ind,] ;
		N	<- nrow(sub_index) ;
	}

	res_table	<- matrix(0, nrow=length(lambda_list), ncol=if (G>1) 6 else 3) ;
	
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
			
			if (nrow(prior_info) > 0) {
				for (j in seq(nrow(prior_info))) {
					ind_1	<- match(prior_list[[k]][j,1], colnames(data)) ;
					ind_2	<- match(prior_list[[k]][j,2], colnames(data)) ;
					info_tmp[[k]][ind_1,ind_2]	<- info_tmp[[k]][ind_2,ind_1]	<- prior_list[[k]][j,3] ;
				}
			}
		}
	
		# Make information of edges with prior
		wp	<- big.matrix(nrow=ncol(data), ncol=ncol(data), init=0) ;
		rownames(wp)	<- colnames(wp)	<- colnames(data) ;
		data_ind	<- combn(ncol(data), 2) ;
		wp_cnt	<- sum(unlist(mclapply(seq(ncol(data_ind)), function(j) {
			ind_1	<- data_ind[1,j] ;
			ind_2	<- data_ind[2,j] ;
			
			pr_get	<- sapply(info_tmp, function(mat) return(!is.na(mat[ind_1,ind_2]))) ;
			if (any(pr_get))	wp[ind_1,ind_2]	<- wp[ind_2,ind_1]	<- 1 ;
			return(wp[ind_1,ind_2]) ;
		}, mc.cores=cores))) ;
		
		wp_bytype	<- numeric(3) ;
		for (j in seq(ncol(data_ind))) {
			ind_1	<- data_ind[1,j] ;
			ind_2	<- data_ind[2,j] ;
			
			if (wp[ind_1,ind_2] == 1) {
				if ((ind_1 %in% ind_disc) & (ind_2 %in% ind_disc)) {
					wp_bytype[3]	<- wp_bytype[3] + 1 ;
				} else if (((ind_1 %in% ind_disc) & !(ind_2 %in% ind_disc)) | (!(ind_1 %in% ind_disc) & (ind_2 %in% ind_disc))) {
					wp_bytype[2]	<- wp_bytype[2] + 1 ;
				} else {
					wp_bytype[1]	<- wp_bytype[1] + 1 ;
				}
			}
		}
		np_bytype	<- numeric(3) ;
		np_bytype[1]	<- p_x*(p_x-1)/2 - wp_bytype[1] ;
		np_bytype[2]	<- (p_x)*(p_y)   - wp_bytype[2] ;
		np_bytype[3]	<- p_y*(p_y-1)/2 - wp_bytype[3] ;
	} else {
		np_bytype	<- numeric(3) ;
		np_bytype[1]	<- p_x*(p_x-1)/2 ;
		np_bytype[2]	<- (p_x)*(p_y) ;
		np_bytype[3]	<- p_y*(p_y-1)/2 ;
	}
	
	# Make estimation array
	net_tmp_1	<- make_net_tmp(data, length(lambda_list)) ;
	net_tmp_2	<- make_net_tmp(data, length(lambda_list)) ;
	
	for (d in seq_along(lambda_list)) {
		for (n in seq(N)) {
			subsample	<- sub_index[n,] ;
			sub_list	<- suppressWarnings(make_MGM_list(X[subsample,,drop=FALSE], Y[subsample,,drop=FALSE], group[subsample])) ;
			
			lambda	<- lambda_list[d] ;
			sub_res		<- learn_tb(sub_list, lambda_intra=rep(lambda, 3), lambda_inter=rep(lambda, 3), 
										with_prior=FALSE, wp=NULL, ..., cores=cores, verbose=verbose) ;
			
			add_skeleton(sub_res, inst_temp, tol_polish, cores) ;
			add_cnt(sub_res[[1]], net_tmp_1, d, tol_polish, cores) ;
			add_cnt(sub_res[[2]], net_tmp_2, d, tol_polish, cores) ;
			
			rm(sub_list, sub_res) ; gc() ;
		}
		
		if (p_x > 1) {
			beta_ind	<- expand.grid(seq(p_x), seq(p_x)) ;
		
			mclapply(seq(nrow(beta_ind)), function(i) {
				s1	<- beta_ind[i,1] ;
				s2	<- beta_ind[i,2] ;
				
				for (g in seq(G)) {
					inst_temp$beta_sum[[g]][s1,s2]	<- inst_temp$beta_sum[[g]][s1,s2] / N ;
					inst_temp$beta_sum[[g]][s1,s2]	<- inst_temp$beta_sum[[g]][s1,s2] * (1 - inst_temp$beta_sum[[g]][s1,s2]) ;
				}
				
				if (G > 1) {
					for (g2 in seq(2,G)) {
						for (g1 in seq(1,g2-1)) {
							base_tmp	<- sum(seq(0,g2-2)) ;

							inst_temp$betadiff_sum[[base_tmp + g1]][s1,s2]	<- inst_temp$betadiff_sum[[base_tmp + g1]][s1,s2] / N ;
							inst_temp$betadiff_sum[[base_tmp + g1]][s1,s2]	<- inst_temp$betadiff_sum[[base_tmp + g1]][s1,s2] * (1 - inst_temp$betadiff_sum[[base_tmp + g1]][s1,s2]) ;	
						}
					}
				}
				
				return(NULL) ;
			}, mc.cores=cores) ;
			
			beta_ind	<- combn(p_x, 2) ;
			
			res_table[d,1]	<- res_table[d,1] + sum(unlist(mclapply(seq(ncol(beta_ind)), function(j) {
				s1	<- beta_ind[1,j] ;
				s2	<- beta_ind[2,j] ;
				
				var1	<- colnames(X)[s1] ;
				var2	<- colnames(X)[s2] ;
			
				return_tmp	<- 0 ;
				
				for (g in seq(G)) {
					if (!with_prior || wp[var1,var2] != 1) return_tmp	<- return_tmp + inst_temp$beta_sum[[g]][s1,s2] ;
				}
				
				return(return_tmp) ;
			}, mc.cores=cores))) ;

			if (G > 1) {
				res_table[d,4]	<- res_table[d,4] + sum(unlist(mclapply(seq(ncol(beta_ind)), function(j) {
					s1	<- beta_ind[1,j] ;
					s2	<- beta_ind[2,j] ;

					var1	<- colnames(X)[s1] ;
					var2	<- colnames(X)[s2] ;

					return_tmp	<- 0 ;

					for (g2 in seq(2,G)) {
						for (g1 in seq(1,g2-1)) {
							base_tmp	<- sum(seq(0,g2-2)) ;
							return_tmp	<- return_tmp + inst_temp$betadiff_sum[[base_tmp + g1]][s1,s2] ;
						}
					}

					return(return_tmp) ;
				}, mc.cores=cores))) ;
			}
			
			rm(beta_ind) ; gc() ;
		}
		
		if (p_x > 0 & p_y > 0) {
			rho_ind		<- expand.grid(seq(p_y), seq(p_x)) ;
			
			res_table[d,2]	<- res_table[d,2] + sum(unlist(mclapply(seq(nrow(rho_ind)), function(i) {
				r	<- rho_ind[i,1] ;
				s	<- rho_ind[i,2] ;
				
				var1	<- colnames(Y)[r] ;
				var2	<- colnames(X)[s] ;
				
				return_tmp	<- 0 ;
				
				for (g in seq(G)) {
					if (!with_prior || wp[var1,var2] != 1) {
						inst_temp$rho_sum[[g]][r,s]	<- inst_temp$rho_sum[[g]][r,s] / N ;
						inst_temp$rho_sum[[g]][r,s]	<- inst_temp$rho_sum[[g]][r,s] * (1 - inst_temp$rho_sum[[g]][r,s]) ;
					
						return_tmp	<- return_tmp + inst_temp$rho_sum[[g]][r,s] ;
					}
				}

				return(return_tmp) ;
			}, mc.cores=cores))) ;


			if (G > 1) {
				res_table[d,5]	<- res_table[d,5] + sum(unlist(mclapply(seq(nrow(rho_ind)), function(i) {
					r	<- rho_ind[i,1] ;
					s	<- rho_ind[i,2] ;

					var1	<- colnames(Y)[r] ;
					var2	<- colnames(X)[s] ;

					return_tmp	<- 0 ;

					for (g2 in seq(2,G)) {
						for (g1 in seq(1,g2-1)) {
							base_tmp	<- sum(seq(0,g2-2)) ;

							inst_temp$rhodiff_sum[[base_tmp + g1]][r,s]	<- inst_temp$rhodiff_sum[[base_tmp + g1]][r,s] / N ;
							inst_temp$rhodiff_sum[[base_tmp + g1]][r,s]	<- inst_temp$rhodiff_sum[[base_tmp + g1]][r,s] * (1 - inst_temp$rhodiff_sum[[base_tmp + g1]][r,s]) ;	

							return_tmp	<- return_tmp + inst_temp$rhodiff_sum[[base_tmp + g1]][r,s] ;
						}
					}

					return(return_tmp) ;
				}, mc.cores=cores))) ;
			}
				
			rm(rho_ind) ;	gc() ;
		}
		
		if (p_y > 1) {
			phi_ind		<- combn(p_y, 2) ;
			
			res_table[d,3]	<- res_table[d,3] + sum(unlist(mclapply(seq(ncol(phi_ind)), function(j) {
				s1	<- phi_ind[1,j] ;
				s2	<- phi_ind[2,j] ;
				
				var1	<- colnames(Y)[s1] ;
				var2	<- colnames(Y)[s2] ;
				
				return_tmp	<- 0 ;
				
				for (g in seq(G)) {
					if (!with_prior || wp[var1,var2] != 1) {
						inst_temp$phi_sum[[g]][s1,s2]	<- inst_temp$phi_sum[[g]][s1,s2] / N ;
						inst_temp$phi_sum[[g]][s1,s2]	<- inst_temp$phi_sum[[g]][s1,s2] * (1 - inst_temp$phi_sum[[g]][s1,s2]) ;
						
						inst_temp$phi_sum[[g]][s2,s1]	<- inst_temp$phi_sum[[g]][s2,s1] ;
						
						return_tmp	<- return_tmp + inst_temp$phi_sum[[g]][s1,s2] ;
					}
				}

				return(return_tmp) ;
			}, mc.cores=cores))) ;


			if (G > 1) {
				res_table[d,6]	<- res_table[d,6] + sum(unlist(mclapply(seq(ncol(phi_ind)), function(j) {
					s1	<- phi_ind[1,j] ;
					s2	<- phi_ind[2,j] ;

					var1	<- colnames(Y)[s1] ;
					var2	<- colnames(Y)[s2] ;

					return_tmp	<- 0 ;

					for (g2 in seq(2,G)) {
						for (g1 in seq(1,g2-1)) {
							base_tmp	<- sum(seq(0,g2-2)) ;

							inst_temp$phidiff_sum[[base_tmp + g1]][s1,s2]	<- inst_temp$phidiff_sum[[base_tmp + g1]][s1,s2] / N ;
							inst_temp$phidiff_sum[[base_tmp + g1]][s1,s2]	<- inst_temp$phidiff_sum[[base_tmp + g1]][s1,s2] * (1 - inst_temp$phidiff_sum[[base_tmp + g1]][s1,s2]) ;	

							inst_temp$phidiff_sum[[base_tmp + g1]][s2,s1]	<- inst_temp$phidiff_sum[[base_tmp + g1]][s1,s2] ;

							return_tmp	<- return_tmp + inst_temp$phidiff_sum[[base_tmp + g1]][s1,s2] ;
						}
					}

					return(return_tmp) ;
				}, mc.cores=cores))) ;
			}
			
			rm(phi_ind) ;	gc() ;
		}
		
		res_table[d,1]	<- res_table[d,1] / (G*np_bytype[1]) ;
		res_table[d,2]	<- res_table[d,2] / (G*np_bytype[2]) ;
		res_table[d,3]	<- res_table[d,3] / (G*np_bytype[3]) ;
		
		if (G > 1) {
			res_table[d,4]	<- res_table[d,4] / (choose(G,2)*choose(p_x,2)) ;
			res_table[d,5]	<- res_table[d,5] / (choose(G,2)*p_x*p_y) ;
			res_table[d,6]	<- res_table[d,6] / (choose(G,2)*choose(p_y,2)) ;
		}
		
		initialize_instmat(inst_temp, cores) ;
	}
	
	colnames(res_table)	<- if (G > 1) c("Cont-Cont, intra", "Cont-Disc, intra", "Disc-Disc, intra",
							"Cont-Cont, inter", "Cont-Disc, inter", "Disc-Disc, inter") else c("Cont-Cont, intra", "Cont-Disc, intra", "Disc-Disc, intra") ;
	rownames(res_table)	<- as.character(lambda_list) ;
							
	lambda_fin	<- numeric(if (G > 1) 6 else 3) ;
	
	for (j in seq_along(lambda_fin)) {
		if (all(is.nan(res_table[,j]))) next ;
		if (res_table[1,j] > gamma) {
			stop(paste0("Instability cutoff not met for the following edges: ", colnames(res_table)[j], "  Please try higher values of the penalization parameters")) ;
		}
		
		if (all(res_table[,j] <= gamma)) {
			lambda_fin[j]	<- lambda_list[length(lambda_list)] ;
			warning(paste0("All of the parameters met the instability cutoff for: ", colnames(res_table)[j], "  Using smaller values of the parameters or checking the result table may be helpful")) ;
		} else {
			for (d in seq(2, length(lambda_list))) {
				if (res_table[d,j] > gamma | res_table[d,j] < res_table[d-1,j]) {
					lambda_fin[j]	<- lambda_list[d-1] ;
					break ;
				}
			}
		}
	}
	
	names(lambda_fin)	<- colnames(res_table) ;
	
	if (with_prior) {
		pq	<- ncol(data) ;	J	<- length(lambda_list) ;
		net_sum_1	<- big.matrix(nrow=pq, ncol=pq, init=0) ;
		net_sum_2	<- big.matrix(nrow=pq, ncol=pq, init=0) ;
		rownames(net_sum_1)	<- colnames(net_sum_1)	<- colnames(data) ;
		rownames(net_sum_2)	<- colnames(net_sum_2)	<- colnames(data) ;
		
		sumup_cnt(net_tmp_1, net_sum_1, cores) ;
		sumup_cnt(net_tmp_2, net_sum_2, cores) ;
		
		# Calculate the confidence of prior sources
		
		src_conf_1	<- numeric(length=length(prior_list)) ;
		src_conf_2	<- numeric(length=length(prior_list)) ;
		
		for (k in seq_along(prior_list)) {
			prior_info	<- prior_list[[k]] ;
			
			if (nrow(prior_info) > 0) {
				for (j in seq(nrow(prior_info))) {
					ind_1	<- match(prior_info[j,1], colnames(net_sum_1)) ;
					ind_2	<- match(prior_info[j,2], colnames(net_sum_2)) ;

					src_conf_1[k]	<- src_conf_1[k] + abs(N*prior_info[j,3] - net_sum_1[ind_1,ind_2]/J) ;
					src_conf_2[k]	<- src_conf_2[k] + abs(N*prior_info[j,3] - net_sum_2[ind_1,ind_2]/J) ;
				}

				src_conf_1[k]	<- src_conf_1[k] / nrow(prior_info) ;
				src_conf_2[k]	<- src_conf_2[k] / nrow(prior_info) ;
			} else {
				src_conf_1[k]	<- src_conf_2[k]	<- Inf;
			}
		}

		# Permute & re-calculate the confidences

		data_ind	<- combn(pq, 2) ;
		
		src_mat_1	<- matrix(0, nrow=perm, ncol=length(prior_list)) ;
		src_mat_2	<- matrix(0, nrow=perm, ncol=length(prior_list)) ;
		
		for (pp in seq(perm)) {
			for (k in seq_along(prior_list)) {
				prior_info		<- prior_list[[k]] ;
				
				if (nrow(prior_info) > 0) {
					ind_permuted	<- sample(seq(nrow(prior_info)), nrow(prior_info), replace=FALSE) ;

					for (j in seq(nrow(prior_info))) {
		#				ind_1	<- match(prior_info[ind_permuted[j],1], colnames(net_sum_1)) ;
		#				ind_2	<- match(prior_info[ind_permuted[j],2], colnames(net_sum_1)) ;

						ind_1	<- data_ind[1,ind_permuted[j]] ;
						ind_2	<- data_ind[2,ind_permuted[j]] ;

						src_mat_1[pp,k]	<- src_mat_1[pp,k] + abs(N*prior_info[j,3] - net_sum_1[ind_1,ind_2]/J) ;
						src_mat_2[pp,k]	<- src_mat_2[pp,k] + abs(N*prior_info[j,3] - net_sum_2[ind_1,ind_2]/J) ;
					}

					src_mat_1[pp,k]	<- src_mat_1[pp,k] / nrow(prior_info) ;
					src_mat_2[pp,k]	<- src_mat_2[pp,k] / nrow(prior_info) ;
				}
			}
		}
		
		src_conf_1	<- src_conf_1 / colMeans(src_mat_1) ;
		src_conf_2	<- src_conf_2 / colMeans(src_mat_2) ;
		
		ignore_ind	<- logical(length(prior_list)) ;
		for (m in seq(length(prior_list))) {
			if (src_conf_1[m] %in% c(0,Inf) | src_conf_2[m] %in% c(0,Inf)) {
				ignore_ind[m]	<- TRUE ;
			}
		}
		
		weight_1	<- numeric(length(prior_list)) ;
		weight_2	<- numeric(length(prior_list)) ;
		weight_1[!ignore_ind]	<- (1/src_conf_1[!ignore_ind]) / sum(1/src_conf_1[!ignore_ind]) ;
		weight_2[!ignore_ind]	<- (1/src_conf_2[!ignore_ind]) / sum(1/src_conf_2[!ignore_ind]) ;
		
		data_ind	<- combn(pq, 2) ;
		
		mean_KL_1		<- big.matrix(nrow=pq, ncol=pq, init=NULL) ;
		mean_KL_2		<- big.matrix(nrow=pq, ncol=pq, init=NULL) ;
		rownames(mean_KL_1)	<- rownames(mean_KL_2)	<- colnames(mean_KL_1)	<- colnames(mean_KL_2)	<- colnames(data) ;
		mclapply(seq(ncol(data_ind)), function(j) {
			ind_1	<- data_ind[1,j] ;
			ind_2	<- data_ind[2,j] ;
			
			var_1	<- colnames(data)[ind_1] ;
			var_2	<- colnames(data)[ind_2] ;
			
			for (k in seq_along(prior_list)) {
				if (!is.na(info_tmp[[k]][var_1,var_2])) {
					if (is.na(mean_KL_1[ind_1,ind_2]))	mean_KL_1[ind_1,ind_2]	<- mean_KL_1[ind_2,ind_1]	<- 0 ;
					if (is.na(mean_KL_2[ind_1,ind_2]))	mean_KL_2[ind_1,ind_2]	<- mean_KL_2[ind_2,ind_1]	<- 0 ;
					
					mean_KL_1[ind_1,ind_2]	<- mean_KL_1[ind_1,ind_2] + weight_1[k]*N*info_tmp[[k]][var_1,var_2] ;
					mean_KL_2[ind_1,ind_2]	<- mean_KL_2[ind_1,ind_2] + weight_2[k]*N*info_tmp[[k]][var_1,var_2] ;
				}
			}
			
			return(NULL) ;
		}, mc.cores=cores) ;
		
		var_KL_1		<- big.matrix(nrow=pq, ncol=pq, init=NULL) ;
		var_KL_2		<- big.matrix(nrow=pq, ncol=pq, init=NULL) ;
		rownames(var_KL_1)	<- rownames(var_KL_2)	<- colnames(var_KL_1)	<- colnames(var_KL_2)	<- colnames(data) ;
		mclapply(seq(ncol(data_ind)), function(j) {
			ind_1	<- data_ind[1,j] ;
			ind_2	<- data_ind[2,j] ;
			
			var_1	<- colnames(data)[ind_1] ;
			var_2	<- colnames(data)[ind_2] ;
			
			for (k in seq_along(prior_list)) {
				if (!is.na(info_tmp[[k]][var_1,var_2])) {
					if (is.na(var_KL_1[ind_1,ind_2]))	var_KL_1[ind_1,ind_2]	<- var_KL_1[ind_2,ind_1]	<- 0 ;
					if (is.na(var_KL_2[ind_1,ind_2]))	var_KL_2[ind_1,ind_2]	<- var_KL_2[ind_2,ind_1]	<- 0 ;
					
					var_KL_1[ind_1,ind_2]	<- var_KL_1[ind_1,ind_2] + weight_1[k]*(src_conf_1[k]^2 + (mean_KL_1[ind_1,ind_2] - N*info_tmp[[k]][var_1,var_2])^2) ;
					var_KL_2[ind_1,ind_2]	<- var_KL_2[ind_1,ind_2] + weight_2[k]*(src_conf_2[k]^2 + (mean_KL_2[ind_1,ind_2] - N*info_tmp[[k]][var_1,var_2])^2) ; ;
				}
			}
			
			return(NULL) ;
		}, mc.cores=cores) ;
		
		mean_post_1		<- big.matrix(nrow=pq, ncol=pq, init=NULL) ;
		mean_post_2		<- big.matrix(nrow=pq, ncol=pq, init=NULL) ;
		rownames(mean_post_1)	<- rownames(mean_post_2)	<- colnames(mean_post_1)	<- colnames(mean_post_2)	<- colnames(data) ;
		mclapply(seq(ncol(data_ind)), function(j) {
			ind_1	<- data_ind[1,j] ;
			ind_2	<- data_ind[2,j] ;
			
			var_1	<- colnames(data)[ind_1] ;
			var_2	<- colnames(data)[ind_2] ;
			
			var_tmp_1	<- (net_sum_1[ind_1,ind_2] / (N*J)) * (1 - net_sum_1[ind_1,ind_2] / (N*J)) * N ;
			var_tmp_2	<- (net_sum_2[ind_1,ind_2] / (N*J)) * (1 - net_sum_2[ind_1,ind_2] / (N*J)) * N ;

			if (!is.na(mean_KL_1[ind_1,ind_2])) {
				mean_post_1[ind_1,ind_2]	<- (mean_KL_1[ind_1,ind_2]*var_tmp_1 + (net_sum_1[ind_1,ind_2]/J)*var_KL_1[ind_1,ind_2]) / (var_KL_1[ind_1,ind_2] + var_tmp_1) ;
			}
			if (!is.na(mean_KL_2[ind_1,ind_2])) {
				mean_post_2[ind_1,ind_2]	<- (mean_KL_2[ind_1,ind_2]*var_tmp_2 + (net_sum_2[ind_1,ind_2]/J)*var_KL_2[ind_1,ind_2]) / (var_KL_2[ind_1,ind_2] + var_tmp_2) ;
			}
			
			return(NULL) ;
		}, mc.cores=cores) ;
		
		var_post_1		<- big.matrix(nrow=pq, ncol=pq, init=NULL) ;
		var_post_2		<- big.matrix(nrow=pq, ncol=pq, init=NULL) ;
		rownames(var_post_1)	<- rownames(var_post_2)	<- colnames(var_post_1)	<- colnames(var_post_2)	<- colnames(data) ;
		mclapply(seq(ncol(data_ind)), function(j) {
			ind_1	<- data_ind[1,j] ;
			ind_2	<- data_ind[2,j] ;
			
			var_1	<- colnames(data)[ind_1] ;
			var_2	<- colnames(data)[ind_2] ;
			
			var_tmp_1	<- (net_sum_1[ind_1,ind_2] / (N*J)) * (1 - net_sum_1[ind_1,ind_2] / (N*J)) * N ;
			var_tmp_2	<- (net_sum_2[ind_1,ind_2] / (N*J)) * (1 - net_sum_2[ind_1,ind_2] / (N*J)) * N ;
			
			if (!is.na(mean_KL_1[ind_1,ind_2])) {
				var_post_1[ind_1,ind_2]	<- (var_KL_1[ind_1,ind_2] * var_tmp_1) / (var_KL_1[ind_1,ind_2] + var_tmp_1) ;
			}
			if (!is.na(mean_KL_2[ind_1,ind_2])) {
				var_post_2[ind_1,ind_2]	<- (var_KL_2[ind_1,ind_2] * var_tmp_2) / (var_KL_2[ind_1,ind_2] + var_tmp_2) ;
			}
			
			return(NULL) ;
		}, mc.cores=cores) ;
		
		score_wp	<- matrix(0, nrow=length(lambda_list), ncol=3) ;
		
		for (d in seq_along(lambda_list)) {

			score_wp[d,1]	<- sum(unlist(mclapply(seq(ncol(data_ind)), function(j) {
				ind_1	<- data_ind[1,j] ;
				ind_2	<- data_ind[2,j] ;
				
				var_1	<- colnames(data)[ind_1] ;
				var_2	<- colnames(data)[ind_2] ;
				
				if (!(var_1 %in% colnames(X) & var_2 %in% colnames(X))) return(0) ;
				
				return_tmp	<- 0 ;
				
				cnt_1   <- if(is.na(net_tmp_1[[d]][ind_1,ind_2])) 0 else net_tmp_1[[d]][ind_1,ind_2] ;
							cnt_2   <- if(is.na(net_tmp_2[[d]][ind_1,ind_2])) 0 else net_tmp_2[[d]][ind_1,ind_2] ;

				if (!is.na(mean_KL_1[ind_1,ind_2])) {
					a_tmp_1	<- pnorm(cnt_1 + eps, mean=mean_post_1[ind_1,ind_2], sd=sqrt(var_post_1[ind_1,ind_2])) ;
					a_tmp_1	<- a_tmp_1 - pnorm(cnt_1 - eps, mean=mean_post_1[ind_1,ind_2], sd=sqrt(var_post_1[ind_1,ind_2])) ;
					
					a_tmp_2	<- pnorm(cnt_2 + eps, mean=mean_post_2[ind_1,ind_2], sd=sqrt(var_post_2[ind_1,ind_2])) ;
					a_tmp_2	<- a_tmp_2 - pnorm(cnt_2 - eps, mean=mean_post_2[ind_1,ind_2], sd=sqrt(var_post_2[ind_1,ind_2])) ;
				} else {
					return(0) ;
				}
				
				inst_tmp_1	<- 2*(cnt_1/N) * (1 - cnt_1/N) ;
				inst_tmp_2	<- 2*(cnt_2/N) * (1 - cnt_2/N) ;
				return(a_tmp_1*(1-2*inst_tmp_1) + a_tmp_2*(1-2*inst_tmp_2)) ;
			}, mc.cores=cores))) ;
			
			score_wp[d,2]	<- sum(unlist(mclapply(seq(ncol(data_ind)), function(j) {
				ind_1	<- data_ind[1,j] ;
				ind_2	<- data_ind[2,j] ;
				
				var_1	<- colnames(data)[ind_1] ;
				var_2	<- colnames(data)[ind_2] ;
				
				if (!(var_1 %in% colnames(X) & var_2 %in% colnames(Y)) & !(var_1 %in% colnames(Y) & var_2 %in% colnames(X))) return(0) ;
				
				return_tmp	<- 0 ;
				
				cnt_1   <- if(is.na(net_tmp_1[[d]][ind_1,ind_2])) 0 else net_tmp_1[[d]][ind_1,ind_2] ;
							cnt_2   <- if(is.na(net_tmp_2[[d]][ind_1,ind_2])) 0 else net_tmp_2[[d]][ind_1,ind_2] ;

				if (!is.na(mean_KL_1[ind_1,ind_2])) {
					a_tmp_1	<- pnorm(cnt_1 + eps, mean=mean_post_1[ind_1,ind_2], sd=sqrt(var_post_1[ind_1,ind_2])) ;
					a_tmp_1	<- a_tmp_1 - pnorm(cnt_1 - eps, mean=mean_post_1[ind_1,ind_2], sd=sqrt(var_post_1[ind_1,ind_2])) ;
					
					a_tmp_2	<- pnorm(cnt_2 + eps, mean=mean_post_2[ind_1,ind_2], sd=sqrt(var_post_2[ind_1,ind_2])) ;
					a_tmp_2	<- a_tmp_2 - pnorm(cnt_2 - eps, mean=mean_post_2[ind_1,ind_2], sd=sqrt(var_post_2[ind_1,ind_2])) ;
				} else {
					return(0) ;
				}
				
				inst_tmp_1	<- 2*(cnt_1/N) * (1 - cnt_1/N) ;
				inst_tmp_2	<- 2*(cnt_2/N) * (1 - cnt_2/N) ;
				return(a_tmp_1*(1-2*inst_tmp_1) + a_tmp_2*(1-2*inst_tmp_2)) ;
			}, mc.cores=cores))) ;
			
			score_wp[d,3]	<- sum(unlist(mclapply(seq(ncol(data_ind)), function(j) {
				ind_1	<- data_ind[1,j] ;
				ind_2	<- data_ind[2,j] ;
				
				var_1	<- colnames(data)[ind_1] ;
				var_2	<- colnames(data)[ind_2] ;
				
				if (!(var_1 %in% colnames(Y) & var_2 %in% colnames(Y))) return(0);
				
				return_tmp	<- 0 ;
				
				cnt_1   <- if(is.na(net_tmp_1[[d]][ind_1,ind_2])) 0 else net_tmp_1[[d]][ind_1,ind_2] ;
							cnt_2   <- if(is.na(net_tmp_2[[d]][ind_1,ind_2])) 0 else net_tmp_2[[d]][ind_1,ind_2] ;

				if (!is.na(mean_KL_1[ind_1,ind_2])) {
					a_tmp_1	<- pnorm(cnt_1 + eps, mean=mean_post_1[ind_1,ind_2], sd=sqrt(var_post_1[ind_1,ind_2])) ;
					a_tmp_1	<- a_tmp_1 - pnorm(cnt_1- eps, mean=mean_post_1[ind_1,ind_2], sd=sqrt(var_post_1[ind_1,ind_2])) ;
					
					a_tmp_2	<- pnorm(cnt_2 + eps, mean=mean_post_2[ind_1,ind_2], sd=sqrt(var_post_2[ind_1,ind_2])) ;
					a_tmp_2	<- a_tmp_2 - pnorm(cnt_2 - eps, mean=mean_post_2[ind_1,ind_2], sd=sqrt(var_post_2[ind_1,ind_2])) ;
				} else {
					return(0) ;
				}
				
				inst_tmp_1	<- 2*(cnt_1/N) * (1 - cnt_1/N) ;
				inst_tmp_2	<- 2*(cnt_2/N) * (1 - cnt_2/N) ;
				return(a_tmp_1*(1-2*inst_tmp_1) + a_tmp_2*(1-2*inst_tmp_2)) ;
			}, mc.cores=cores))) ;
		}
		
		lambda_fin_prior	<- numeric(3) ;
		
		for (j in seq(3)) {
			if (all(is.nan(score_wp[,j]))) {
				lambda_fin_prior[j] 	<- lambda_list[length(lambda_list)] ;
			} else {
				lambda_fin_prior[j]		<- lambda_list[max(which(score_wp[,j] == min(score_wp[,j])))] ;
			}
		}
		
		colnames(score_wp)	<- names(lambda_fin_prior)	<- colnames(res_table)[1:3] ;
	}
	
	if (!with_prior) {
		return(list(np = list(instability=res_table, lambda_recommended=lambda_fin))) ;
	} else {
		return(list(np = list(instability=res_table, lambda_recommended=lambda_fin), 
					wp = list(score=score_wp, lambda_recommended = lambda_fin_prior))) ;
	}
}

# Make a template for instabilities
make_inst_tmp	<- function(MGM_list) {
	G	<- length(MGM_list) ;

	p_x	<- ncol(MGM_list[[1]]$Continuous) ;	p_y	<- ncol(MGM_list[[1]]$Discrete) ;
	inst_template	<- list() ;
	
	inst_template$beta_sum	<- list() ;
	inst_template$rho_sum	<- list() ;
	inst_template$phi_sum	<- list() ;

	for (g in seq(G)) {
		inst_template$beta_sum[[g]]	<- big.matrix(nrow=p_x, ncol=p_x, init=0) ;
		inst_template$rho_sum[[g]]	<- big.matrix(nrow=p_y, ncol=p_x, init=0) ;
		inst_template$phi_sum[[g]]	<- big.matrix(nrow=p_y, ncol=p_y, init=0) ;
	}

	if (G >= 2) {
		inst_template$betadiff_sum	<- list() ;
		inst_template$rhodiff_sum	<- list() ;
		inst_template$phidiff_sum	<- list() ;
	
		for (g in seq(choose(G,2))) {
			inst_template$betadiff_sum[[g]]	<- big.matrix(nrow=p_x, ncol=p_x, init=0) ;
			inst_template$rhodiff_sum[[g]]	<- big.matrix(nrow=p_y, ncol=p_x, init=0) ;
			inst_template$phidiff_sum[[g]]	<- big.matrix(nrow=p_y, ncol=p_y, init=0) ;
		}
	}
	
	return(inst_template) ;
}

initialize_instmat	<- function(inst_temp, cores) {
	G1	<- length(inst_temp$beta_sum) ;
	G2	<- length(inst_temp$betadiff_sum) ;
	
	for (g1 in seq(G1)) {
		p_x	<- nrow(inst_temp$beta_sum[[g1]]) ;
		p_y	<- nrow(inst_temp$rho_sum[[g1]]) ;
		
		if (p_x > 0) {
			beta_ind	<- expand.grid(seq(p_x), seq(p_x)) ;
			
			mclapply(seq(nrow(beta_ind)), function(i) {
				s1	<- beta_ind[i,1] ;
				s2	<- beta_ind[i,2] ;
				
				inst_temp$beta_sum[[g1]][s1,s2]	<- 0 ;
				
				return(NULL) ;
			}, mc.cores=cores) ;
		}
		
		if (p_x > 0 & p_y > 0) {
			rho_ind		<- expand.grid(seq(p_y), seq(p_x)) ;
			
			mclapply(seq(nrow(rho_ind)), function(i) {
				r	<- rho_ind[i,1] ;
				s	<- rho_ind[i,2] ;
				
				inst_temp$rho_sum[[g1]][r,s]	<- 0 ;
				
				return(NULL) ;
			}, mc.cores=cores) ;
		}
		
		if (p_y > 0) {
			phi_ind		<- expand.grid(seq(p_y), seq(p_y)) ;
			
			mclapply(seq(nrow(phi_ind)), function(i) {
				r1	<- rho_ind[i,1] ;
				r2	<- rho_ind[i,2] ;
				
				inst_temp$phi_sum[[g1]][r1,r2]	<- 0 ;
				
				return(NULL) ;
			}, mc.cores=cores) ;
		}
	}
	
	if (G2 > 0) {
		for (g2 in seq(G2)) {
			p_x	<- nrow(inst_temp$betadiff_sum[[g2]]) ;
			p_y	<- nrow(inst_temp$rhodiff_sum[[g2]]) ;
		
			if (p_x > 0) {
				beta_ind	<- expand.grid(seq(p_x), seq(p_x)) ;
			
				mclapply(seq(nrow(beta_ind)), function(i) {
					s1	<- beta_ind[i,1] ;
					s2	<- beta_ind[i,2] ;
				
					inst_temp$betadiff_sum[[g2]][s1,s2]	<- 0 ;
				
					return(NULL) ;
				}, mc.cores=cores) ;
			}
		
			if (p_x > 0 & p_y > 0) {
				rho_ind		<- expand.grid(seq(p_y), seq(p_x)) ;
			
				mclapply(seq(nrow(rho_ind)), function(i) {
					r	<- rho_ind[i,1] ;
					s	<- rho_ind[i,2] ;
				
					inst_temp$rhodiff_sum[[g2]][r,s]	<- 0 ;
				
					return(NULL) ;
				}, mc.cores=cores) ;
			}
		
			if (p_y > 0) {
				phi_ind		<- expand.grid(seq(p_y), seq(p_y)) ;
			
				mclapply(seq(nrow(phi_ind)), function(i) {
					r1	<- rho_ind[i,1] ;
					r2	<- rho_ind[i,2] ;
				
					inst_temp$phidiff_sum[[g2]][r1,r2]	<- 0 ;
					
					return(NULL) ;
				}, mc.cores=cores) ;
			}
		}
	}
	
	return(NULL) ;
}

make_net_tmp	<- function(data, R) {
	pq	<- ncol(data) ;
	net_tmp	<- rep(list(NULL), R) ;
	
	for (r in seq(R)) {
		net_tmp[[r]]	<- big.matrix(nrow=pq, ncol=pq, init=NA) ;
		rownames(net_tmp[[r]])	<- colnames(net_tmp[[r]])	<- colnames(data) ;
	}
	
	return(net_tmp) ;
}

add_cnt	<- function(MGM_input, net_temp, net_ind, tol_polish=1e-12, cores) {
	X		<- MGM_input$Continuous ;
	Y		<- MGM_input$Discrete ;
	Y_dummy	<- MGM_input$Dummy ;

	p_x	<- ncol(X) ;	p_y	<- ncol(Y) ;
	L_sum		<- ncol(Y_dummy) ;
	L_cumsum	<- c(0, cumsum(MGM_input$Levels)) ;
	
	# Beta
	if (p_x > 1) {
		beta_ind	<- expand.grid(seq(p_x), seq(p_x)) ;
		
		mclapply(seq(nrow(beta_ind)), function(i) {
			s1	<- beta_ind[i,1] ;
			s2	<- beta_ind[i,2] ;
			
			if (abs(MGM_input$beta[s1,s2]) > tol_polish) {
				var_1	<- colnames(X)[s1] ;
				var_2	<- colnames(X)[s2] ;
				
				ind_1	<- match(var_1, colnames(net_temp[[net_ind]])) ;
				ind_2	<- match(var_2, colnames(net_temp[[net_ind]])) ;
				
				if (is.na(net_temp[[net_ind]][ind_1,ind_2])) {
					net_temp[[net_ind]][ind_1,ind_2]	<- 1 ;
				} else {
					net_temp[[net_ind]][ind_1,ind_2]	<- net_temp[[net_ind]][ind_1,ind_2] + 1 ;
				}
			}
			return(NULL) ;
		}, mc.cores=cores) ;
		
		rm(beta_ind) ;	gc() ;
	}
	
	# Rho
	if (p_x > 0 & p_y > 0) {
		rho_ind		<- expand.grid(seq(p_y), seq(p_x)) ;
		
		mclapply(seq(nrow(rho_ind)), function(i) {
			r	<- rho_ind[i,1] ;
			s	<- rho_ind[i,2] ;
			
			ind_tmp		<- (L_cumsum[r]+1):(L_cumsum[r+1]) ;
			rho_tmp		<- MGM_input$rho[ind_tmp, s] ;
			
			if (sqrt(mean(rho_tmp^2)) > tol_polish) {
				var_1	<- colnames(Y)[r] ;
				var_2	<- colnames(X)[s] ;
				
				ind_1	<- match(var_1, colnames(net_temp[[net_ind]])) ;
				ind_2	<- match(var_2, colnames(net_temp[[net_ind]])) ;
				
				if (is.na(net_temp[[net_ind]][ind_1,ind_2])) {
					net_temp[[net_ind]][ind_1,ind_2]	<- 1 ;
				} else {
					net_temp[[net_ind]][ind_1,ind_2]	<-  net_temp[[net_ind]][ind_1,ind_2] + 1 ;
				}
			}
			
			return(NULL) ;
		}, mc.cores=cores) ;
			
		rm(rho_ind) ;	gc() ;
	}
	
	# Phi
	if (p_y > 1) {
		phi_ind		<- combn(p_y, 2) ;
		
		mclapply(seq(ncol(phi_ind)), function(j) {
			s1	<- phi_ind[1,j] ;
			s2	<- phi_ind[2,j] ;

			ind_tmp_1	<- (L_cumsum[s1]+1):(L_cumsum[s1+1]) ;
			ind_tmp_2	<- (L_cumsum[s2]+1):(L_cumsum[s2+1]) ;
				
			phi_tmp		<- MGM_input$phi[ind_tmp_1, ind_tmp_2] ;
			if (sqrt(mean(phi_tmp^2)) > tol_polish) {
				var_1	<- colnames(Y)[s1] ;
				var_2	<- colnames(Y)[s2] ;
				
				ind_1	<- match(var_1, colnames(net_temp[[net_ind]])) ;
				ind_2	<- match(var_2, colnames(net_temp[[net_ind]])) ;
				
				if (is.na(net_temp[[net_ind]][ind_1,ind_2])) {
					net_temp[[net_ind]][ind_1,ind_2]	<- net_temp[[net_ind]][ind_2,ind_1]	<- 1 ;
				} else {
					net_temp[[net_ind]][ind_1,ind_2]	<- net_temp[[net_ind]][ind_1,ind_2] + 1 ;
					net_temp[[net_ind]][ind_2,ind_1]	<- net_temp[[net_ind]][ind_2,ind_1] + 1 ;
				}
			}
			
			return(NULL) ;
		}, mc.cores=cores) ;
		
		rm(phi_ind) ;	gc() ;
	}
	
	return(NULL) ;
}

sumup_cnt	<- function(net_temp, sum_mat, cores) {
	pq	<- nrow(sum_mat) ;
	data_ind	<- combn(pq, 2) ;
	
	mclapply(seq(ncol(data_ind)), function(j) {
		var_1	<- rownames(sum_mat)[data_ind[1,j]] ;
		var_2	<- rownames(sum_mat)[data_ind[2,j]] ;

		ind_1	<- match(var_1, colnames(sum_mat)) ;
		ind_2	<- match(var_2, colnames(sum_mat)) ;
		
		cnt_values	<- sapply(net_temp, function(mat) if (!is.na(mat[ind_1,ind_2])) return(mat[ind_1,ind_2]) else return(mat[ind_2,ind_1] )) ;
		sum_mat[var_1,var_2]	<- sum_mat[var_2,var_1]	<- sum(cnt_values, na.rm=TRUE) ;
		
		return(NULL) ;
	}, mc.cores=cores) ;
	
	return(NULL) ;
}
