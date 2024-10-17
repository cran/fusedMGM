require(plotrix) ;
require(gplots) ;

#' A plot function for a list of MGMs.
#' The output is usually from FMGM main function.
#'
#' This function is written based on R base function 'heatmap'.
#'
#' @param MGM_list A list of graphs from 2 groups. Usually a result of FMGM main function.
#' @param sortby Determines the standard of sorting & dendrograms. Either 1, 2, or "diff" (default). 
#' @param highlight A vector of variable names or indices to highlight
#' @param tol_polish A threshold for the network edge presence
#' @param tol_plot Only network edges above this value will be displayed on the heatmap
#' @param sideColor A named vector determining a sidebar colors. Set NULL to make the colors based on the variable types (discrete/continuous). Default: FALSE (no sidebars)
#' @param distfun A function for the distances between rows/columns
#' @param hclustfun A function for hierarchical clustering
#' @param reorderfun A function of dendrogram and weights for reordering
#' @param margins A numeric vector of 2 numbers for row & column name margins
#' @param cexRow A visual parameter cex for row axis labeling
#' @param cexCol A visual parameter cex for column axis labeling, default to be same as cexRow
#' @param main Main title, default to none
#' @param xlab X-axis title, default to none
#' @param ylab Y-axis title, default to none
#' @param verbose Logical. Should plotting information be printed?
#' @import grDevices
#' @import graphics
#' @import stats
#' @importFrom gplots colorpanel
#' @return None
#' @examples
#' chk <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
#' 
#' if (Sys.info()['sysname'] != 'Linux') {
#'   cores=1L
#' } else {
#'   chk = tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
#'   if (nzchar(chk) && (chk != "false")) {
#'     cores=2L
#'   } else {
#'     cores=parallel::detectCores() - 1 ;
#'   }
#' }
#' 
#' \dontrun{
#' data(data_all) ;  # Example 500-by-100 simulation data
#' data(ind_disc) ;
#' 
#' group <- rep(c(1,2), each=250) ;
#' names(group) <- seq(500) ;
#' 
#' res_FMGM <- FMGM_mc(data_all, ind_disc, group, 
#'                     lambda_intra=c(0.2,0.15,0.1), lambda_inter=c(0.2,0.15,0.1), 
#'                     cores=cores, verbose=TRUE)
#'                     
#' FMGM_plot(res_FMGM)
#' }
#' 
#' \donttest{
#' data(data_mini) ; # Minimal example 500-by-10 simulation data
#' data(ind_disc_mini) ;
#' 
#' group <- rep(c(1,2), each=250) ;
#' names(group) <- rownames(data_mini) ;
#' 
#' res_FMGM_mini <- FMGM_mc(data_mini, ind_disc_mini, group, 
#'                     lambda_intra=c(0.2,0.15,0.1), lambda_inter=c(0.2,0.15,0.1), 
#'                     cores=cores, verbose=TRUE)
#'                     
#' FMGM_plot(res_FMGM_mini)
#' }
#' @export
FMGM_plot <- function(MGM_list, sortby="diff", highlight=c(), tol_polish=1e-12, tol_plot=.01,
					sideColor = FALSE, distfun = dist, hclustfun = hclust,
					reorderfun = function(d,w) reorder(d,w),
					margins = c(2.5,2.5), cexRow = 0.1 + .5/log10(n),
					cexCol = cexRow, main = NULL, xlab = NULL, ylab = NULL,
					verbose = getOption("verbose")) {
					
	# Make each MGM into a single norm matrix
	skel_MGM	<- MGMlist_to_plotskel(MGM_list, tol_polish=tol_polish) ;
	MGM_intra_1	<- skel_MGM$intra[[1]] ;
	MGM_intra_2	<- skel_MGM$intra[[2]] ;
	MGM_inter	<- skel_MGM$inter[[1]] ;
	
	rm(skel_MGM) ;
	
	X	<- MGM_list[[1]]$Continuous ;
	Y	<- MGM_list[[1]]$Discrete ;
	
	p_x	<- ncol(X) ;
	p_y	<- ncol(Y) ;
	n	<- p_x + p_y ;
	varlist	<- c(colnames(X), colnames(Y)) ;
	
	rm(X, Y) ;
	
	if (!isFALSE(sideColor)) {
		if (is.null(sideColor)) {
			sideColor	<- c(rep(1, p_x), rep(2, p_y)) ;
		} else if (length(sideColor != n)) {
			stop("sideColor vector length must match the total number of variables");
		} else  {
			if (is.null(names(sideColor))) {
				warning("sideColor vector is unnamed; the user should be careful with the variable order") ;
			}
			if (!is.numeric(sideColor)) {
				sideColor	<- as.factor(sideColor) ;
				sdInd	<- levels(sideColor) ;
			}
			sideColor	<- as.integer(sideColor) ;
		}
	}
	
	if (identical(as.character(sortby),"1")) {
		hcr		<- hclustfun(distfun(MGM_intra_1)) ;
		Rowv	<- rowMeans(MGM_intra_1, na.rm=TRUE) ;
	} else if (identical(as.character(sortby),"2")) {
		hcr		<- hclustfun(distfun(MGM_intra_2)) ;
		Rowv	<- rowMeans(MGM_intra_2, na.rm=TRUE) ;
	} else if (identical(as.character(sortby),"diff")) {
		hcr		<- hclustfun(distfun(MGM_inter)) ;
		Rowv	<- rowMeans(MGM_inter, na.rm=TRUE) ;
	} else {
		stop ("Argument of 'sortby ' must be one of 1, 2, or 'diff'.") ;
	}
	
	ddr	<- as.dendrogram(hcr) ;
	ddr	<- reorderfun(ddr, Rowv) ;
	rowInd	<- order.dendrogram(ddr) ;
	varlist	<- varlist[rowInd] ;
	
	if (!isFALSE(sideColor)) {
		sideColor	<- sideColor[rowInd] ;
	}
	
	MGM_intra_1	<- MGM_intra_1[rowInd, rowInd] ;
	MGM_intra_2	<- MGM_intra_2[rowInd, rowInd] ;
	MGM_inter	<- MGM_inter[rowInd, rowInd] ;
	
	MGM_intra_1[abs(MGM_intra_1) < tol_plot]        <- 0 ;
	MGM_intra_2[abs(MGM_intra_2) < tol_plot]        <- 0 ;
	MGM_inter[abs(MGM_inter) < tol_plot]            <- 0 ;
	
	lmat	<- rbind(c(NA,4,NA), c(5,1,3), c(NA,2,6)) ;
	lwid	<- c(1,4,4) ;
	lhei	<- c(1 + if(!is.null(main)) 0.2 else 0, 4, 4) ;
	
	if (!isFALSE(sideColor)) {
		lmat	<- rbind(lmat[1,] + 1, c(NA,1,NA), lmat[2:3,] + 1) ;
		lhei	<- c(lhei[1L], 0.2, lhei[2L:3L]) ;
		
		lmat	<- cbind(lmat[,1] + 1, c(NA,NA,1,NA), lmat[,2:3] + 1) ;
		lwid	<- c(lwid[1L], 0.2, lwid[2L:3L]) ;
	}
	
	if (!is.null(main)) {
		lmat	<- rbind(rep(max(lmat, na.rm=TRUE)+1, ncol(lmat)), lmat) ;
		lhei	<- c(0.1, lhei) ;
	}
	
	lmat[is.na(lmat)]	<- 0 ;
	
	dev.hold() ;
	on.exit(dev.flush()) ;
	
	op	<- par(no.readonly = TRUE) ;
	on.exit(par(op), add = TRUE) ;
	
	layout(lmat, widths=lwid, heights=lhei, respect=TRUE) ;
	
	if (!isFALSE(sideColor)) {
		par(mar = c(0.25, 0, 0, 0.25)) ;
		image(rbind(n:1L), col = sideColor, axes=FALSE) ;
		
		par(mar = c(0.25, 0, 0, 0.25)) ;
		image(cbind(1L:n), col = sideColor, axes=FALSE) ;
	}
	
	wr_palette	<- colorpanel(200, "white", "red") ;
	wg_palette	<- colorpanel(200, "white", "green") ;
	wb_palette	<- colorpanel(200, "white", "blue") ;
	scale_1		<- min(1, max(MGM_intra_1)) ;
	scale_2		<- min(1, max(MGM_intra_2)) ;
	scale_diff	<- min(1, max(MGM_inter)) ;
	
	rescale		<- max(scale_1, scale_2, scale_diff) ;
	newind		<- round(length(wr_palette)*rescale) ;
	wr_palette	<- wr_palette[seq(newind)] ;
	wg_palette	<- wg_palette[seq(newind)] ;
	wb_palette	<- wb_palette[seq(newind)] ;
	
	par(mar=c(0.25, 0, 0, 0.25)) ;
	image(1L:n, 1L:n, MGM_inter[,n:1L], 
		xlim = 0.5 + c(0, n), ylim = 0.5 + c(0, n),
		axes = FALSE, xlab = "", ylab = "", col=wr_palette) ;
		
	if (length(highlight) > 0) {
		mat_hl	<- matrix(ncol=2, nrow=0) ;
		for (vv in highlight) {
			ind_hl	<- match(vv, varlist) ;
			if (nrow(mat_hl) == 0 || ind_hl != sum(mat_hl[nrow(mat_hl),])) {
				mat_hl	<- rbind(mat_hl, c(ind_hl, 1L)) ;
			} else {
				mat_hl[nrow(mat_hl), 2]	<- mat_hl[nrow(mat_hl), 2] + 1L ;
			}
		}
		
		for (nn in seq(nrow(mat_hl))) {
			ind_hl	<- mat_hl[nn,1] ;
			num_hl	<- mat_hl[nn,2] ;
			polygon(c(ind_hl -.5 + c(0,num_hl,num_hl,0)),
					c(n+.5,n+.5,0.5,0.5), lwd=.5) ;
		}
	}
		
	par(mar=c(margins[1L],0,0.25,0.25)) ;
	image(1L:n, 1L:n, MGM_intra_1[,n:1L], 
		xlim = 0.5 + c(0, n), ylim = 0.5 + c(0, n),
		axes = FALSE, xlab = "", ylab = "", col=wg_palette) ;
	axis(1, 1L:n, labels=varlist, las=2, line=-.75, 
		tick=0, cex.axis=cexCol) ;
		
	par(mar=c(0.25,0.25,0,margins[2L])) ;
	image(1L:n, 1L:n, MGM_intra_2[,n:1L], 
		xlim = 0.5 + c(0, n), ylim = 0.5 + c(0, n),
		axes = FALSE, xlab = "", ylab = "", col=wb_palette) ;
	axis(4, n:1L, labels=varlist, las=2, line=-.75, 
		tick=0, cex.axis=cexRow) ;
	
#	for (s1 in seq(nrow(MGM_inter))) {
#		for (s2 in seq(ncol(MGM_inter))) {
#			
#		}
#	}

	par(mar=c(0,0,if (!is.null(main)) 1 else 0,0.25)) ;
	plot(rev(ddr), axes=FALSE, xaxs="i", leaflab="none", lwd=.5) ;

	par(mar=c(0.25,0,0,0)) ;
	plot(ddr, horiz=TRUE, axes=FALSE, yaxs="i", leaflab="none", lwd=.5) ;
	
	par(mar=c(margins[1L], 0.25, 0.25, margins[2L])) ;
	image(1L:2L, 1L:2L, matrix(c(0,0,0,0), nrow=2), col="white",
		axes=FALSE, xlab="", ylab="", xlim=c(0,2), ylim=c(0,2)) ;
	polygon(c(0.05,0.05,0.95,0.95), c(0.05,0.95,0.95,0.05)) ;
        polygon(c(1.05,1.05,1.95,1.95), c(0.05,0.95,0.95,0.05)) ;
        polygon(c(0.05,0.05,0.95,0.95), c(1.05,1.95,1.95,1.05)) ;
        polygon(c(1.05,1.05,1.95,1.95), c(1.05,1.95,1.95,1.05)) ;
#       abline(v=1) ; abline(h=1) ;
        text(c(0.5,0.5,1.5), c(1.5,0.5,1.5),
             labels=c("Network\ndifference", "Network in\ngroup 1", "Network in\ngroup 2"), cex=2*cexRow) ;
        legend(1.5, 0.5, if (exists("sdInd")) sdInd else seq(length(unique(sideColor))),
                col = seq(length(unique(sideColor))), pch = 15,
                adj = 0.5, bty = "n", xjust=.5, yjust=.5, cex=1.5*cexRow) ;
	
	
	if (!is.null(main)) {
		par(xpd=NA, mar=c(0,0,0,0)) ;
		image(1L:2L, 1L:2L, matrix(c(0,0,0,0), nrow=2), col="white",
			axes=FALSE, xlab="", ylab="", xlim=c(0.5,2.5), ylim=c(0.5,2.5)) ;
		title(main, cex.main=op[["cex.main"]], line=-.75) ;
	}
}


# Convert an MGM to a skeleton network
# Applied to inter-node elements only
MGM_to_plotskel		<- function(MGM_input, tol_polish=1e-12) {
	X		<- MGM_input$Continuous ;
	Y		<- MGM_input$Discrete ;
	Y_dummy	<- MGM_input$Dummy ;

	p_x	<- ncol(X) ;	p_y	<- ncol(Y) ;
	L_sum		<- ncol(Y_dummy) ;
	L_cumsum	<- c(0, cumsum(MGM_input$Levels)) ;
	
	varlist_X	<- colnames(X) ;
	varlist_Y	<- colnames(Y) ;
	
	# Beta
	if (p_x > 1) {
		skel_beta	<- matrix(0, nrow=p_x, ncol=p_x) ;
		for (s2 in seq(2,p_x)) {
			for (s1 in seq(1,s2-1)) {
				abs_tmp	<- abs(MGM_input$beta[s1,s2]) ;
				if (abs_tmp > tol_polish) {
					skel_beta[s1,s2]	<- skel_beta[s2,s1]	<- abs_tmp ;
				}
			}
		}
	#	skel_beta[abs(MGM_input$beta) > tol_polish]	<- 1 ;
	} else {
		skel_beta	<- NULL ;
	}
	
	# Rho
	if (p_x > 0 & p_y > 0) {
		skel_rho	<- matrix(0, nrow=p_y, ncol=p_x) ;
		for (s in seq(p_x)) {
			for (r in seq(p_y)) {
				ind_tmp		<- (L_cumsum[r]+1):(L_cumsum[r+1]) ;
				
				rho_tmp		<- MGM_input$rho[ind_tmp, s] ;
				rmsd_tmp	<- sqrt(mean(rho_tmp^2)) ;
				if (rmsd_tmp > tol_polish) {
					skel_rho[r,s]	<- rmsd_tmp ;
				}
			}
		}
	} else {
		skel_rho	<- NULL ;
	}
	
	# Phi
	if (p_y > 1) {
		skel_phi	<- matrix(0, nrow=p_y, ncol=p_y) ;
		for (s2 in seq(2,p_y)) {
			for (s1 in seq(1,s2-1)) {
				ind_tmp_1	<- (L_cumsum[s1]+1):(L_cumsum[s1+1]) ;
				ind_tmp_2	<- (L_cumsum[s2]+1):(L_cumsum[s2+1]) ;
				
				phi_tmp		<- MGM_input$phi[ind_tmp_1, ind_tmp_2] ;
				rmsd_tmp	<- sqrt(mean(phi_tmp^2)) ;
				if (rmsd_tmp > tol_polish) {
					skel_phi[s1,s2]	<- skel_phi[s2,s1]	<- rmsd_tmp ;
				}
			}
		}
	}
	
	skel_all	<- matrix(0, nrow=p_x+p_y, ncol = p_x+p_y) ;
	rownames(skel_all)	<- colnames(skel_all)	<- c(varlist_X, varlist_Y) ;
	
	skel_all[1:p_x, 1:p_x]	<- skel_beta ;
	skel_all[1:p_x, (p_x+1):(p_x+p_y)]	<- t(skel_rho) ;
	skel_all[(p_x+1):(p_x+p_y), 1:p_x]	<- skel_rho ;
	skel_all[(p_x+1):(p_x+p_y), (p_x+1):(p_x+p_y)]	<- skel_phi ;
	
	return(skel_all) ;
}

# Indicate edges with differetial interactions between 2 MGMs
MGM_to_plotskel_diff	<- function(MGM_1, MGM_2, tol_polish=1e-12) {
	X		<- MGM_1$Continuous ;
	Y		<- MGM_1$Discrete ;
	Y_dummy	<- MGM_1$Dummy ;

	p_x	<- ncol(X) ;	p_y	<- ncol(Y) ;
	L_sum		<- ncol(Y_dummy) ;
	L_cumsum	<- c(0, cumsum(MGM_1$Levels)) ;
	
	varlist_X	<- colnames(X) ;
	varlist_Y	<- colnames(Y) ;
	
	# Beta
	if (p_x > 1) {
		skel_beta	<- matrix(0, nrow=p_x, ncol=p_x) ;
		for (s2 in seq(2,p_x)) {
			for (s1 in seq(1,s2-1)) {
				abs_tmp	<- abs(MGM_1$beta[s1,s2] - MGM_2$beta[s1,s2]) ;
				if (abs_tmp > tol_polish) {
					skel_beta[s1,s2]	<- skel_beta[s2,s1]	<- abs_tmp ;
				}
			}
		}
	#	skel_beta	<- sign(MGM_1$beta - MGM_2$beta) ;
	#	skel_beta[abs(MGM_1$beta - MGM_2$beta) <= tol_polish]	<- 0 ;
	} else {
		skel_beta	<- NULL ;
	}
	
	# Rho
	if (p_x > 0 & p_y > 0) {
		skel_rho	<- matrix(0, nrow=p_y, ncol=p_x) ;
		for (s in seq(p_x)) {
			for (r in seq(p_y)) {
				ind_tmp		<- (L_cumsum[r]+1):(L_cumsum[r+1]) ;
				
				rho_tmp_1	<- MGM_1$rho[ind_tmp, s] ;
				rho_tmp_2	<- MGM_2$rho[ind_tmp, s] ;
				rmsd_tmp	<- sqrt(mean((rho_tmp_1 - rho_tmp_2)^2)) ;
				if (rmsd_tmp > tol_polish) {
					skel_rho[r,s]	<- rmsd_tmp ;
				}
			}
		}
	} else {
		skel_rho	<- NULL ;
	}
	
	# Phi
	if (p_y > 1) {
		skel_phi	<- matrix(0, nrow=p_y, ncol=p_y) ;
		for (s2 in seq(2,p_y)) {
			for (s1 in seq(1,s2-1)) {
				ind_tmp_1	<- (L_cumsum[s1]+1):(L_cumsum[s1+1]) ;
				ind_tmp_2	<- (L_cumsum[s2]+1):(L_cumsum[s2+1]) ;
				
				phi_tmp_1	<- MGM_1$phi[ind_tmp_1, ind_tmp_2] ;
				phi_tmp_2	<- MGM_2$phi[ind_tmp_1, ind_tmp_2] ;
				
				rmsd_tmp	<- sqrt(mean((phi_tmp_1 - phi_tmp_2)^2)) ;
				if (rmsd_tmp > tol_polish) {
					skel_phi[s1,s2]	<- skel_phi[s2,s1]	<- rmsd_tmp ;
				}
			}
		}
	}
	
	skel_all	<- matrix(0, nrow=p_x+p_y, ncol = p_x+p_y) ;
	rownames(skel_all)	<- colnames(skel_all)	<- c(varlist_X, varlist_Y) ;
	
	skel_all[1:p_x, 1:p_x]	<- skel_beta ;
	skel_all[1:p_x, (p_x+1):(p_x+p_y)]	<- t(skel_rho) ;
	skel_all[(p_x+1):(p_x+p_y), 1:p_x]	<- skel_rho ;
	skel_all[(p_x+1):(p_x+p_y), (p_x+1):(p_x+p_y)]	<- skel_phi ;
	
	return(skel_all) ;
}

# Make skeleton networks for an MGM list
# Intra: skeleton for edges in each of the MGM
# Inter: skeleton for the differences between the MGMs
MGMlist_to_plotskel	<- function(learn_res, tol_polish=1e-12) {
	G	<- length(learn_res) ;
	skeleton_out	<- vector("list", G) ;
	
	for (g in seq(G)) {
		skeleton_out[[g]]	<- MGM_to_plotskel(learn_res[[g]], tol_polish) ;
	}
	
	skel_diff_out	<- vector("list", choose(G,2)) ;
	
	for (g2 in seq(2,G)) {
		for (g1 in seq(1,g2-1)) {
			skel_diff_out[[sum(0,g2-2) + g1]]	<-  
				MGM_to_plotskel_diff(learn_res[[g1]], learn_res[[g2]], tol_polish) ;
		}
	}
	
	return(list(intra=skeleton_out, inter=skel_diff_out)) ;	
}
