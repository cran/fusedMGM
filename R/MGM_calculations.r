require(fastDummies) ;
require(Rcpp) ;
require(parallel) ;
require(bigmemory) ;
require(bigalgebra) ;
require(biganalytics) ;

# Return square sums of differences & # of parameters compared
rmsd_MGM		<- function(MGM_list_1, MGM_list_2, cores) {

	X		<- MGM_list_1[[1]]$Continuous ;
	Y		<- MGM_list_1[[1]]$Discrete ;
	Y_dummy	<- MGM_list_1[[1]]$Dummy ;

	p_x	<- ncol(X) ;	p_y	<- ncol(Y) ;
	L_sum		<- ncol(Y_dummy) ;
	L_cumsum	<- c(0, cumsum(MGM_list_1[[1]]$Levels)) ;
	
	rm(X, Y, Y_dummy) ;
	
	sumdiff		<- 0 ;
	num_param	<- 0 ;
	
	G	<- length(MGM_list_1) ;
	
	for (g in seq(G)) {
	
		MGM_1	<- MGM_list_1[[g]] ;
		MGM_2	<- MGM_list_2[[g]] ;
	
		# Alpha
		if (p_x > 0) {
			sumdiff	<- sumdiff + sum(unlist(mclapply(seq(p_x), function(s) {
				return((MGM_1$alpha[s,1] - MGM_2$alpha[s,1])^2) ;
			}, mc.cores=cores))) ;
			num_param	<- num_param + p_x ;
		}
		
		# Beta
		if (p_x > 0) {
			## Diagonal
			sumdiff	<- sumdiff + sum(unlist(mclapply(seq(p_x), function(s) {
				return((MGM_1$beta[s,s] - MGM_2$beta[s,s])^2) ;
			}, mc.cores=cores))) ;
			num_param	<- num_param + p_x ;
			
			## Non-diagonal
			if (p_x > 1) {
				beta_ind	<- combn(p_x,2) ;
				sumdiff	<- sumdiff + sum(unlist(mclapply(seq(ncol(beta_ind)), function(j) {
					s1	<- beta_ind[1,j] ;
					s2	<- beta_ind[2,j] ;
					
					return((MGM_1$beta[s1,s2] - MGM_2$beta[s1,s2])^2) ;
				}, mc.cores=cores))) ;
				num_param	<- num_param + ncol(beta_ind) ;
			
				rm(beta_ind) ;	gc() ;
			}
		}
		
		# Rho
		if (p_x > 0 & p_y > 0) {
			rho_ind		<- expand.grid(seq(L_sum), seq(p_x)) ;
			sumdiff	<- sumdiff + sum(unlist(mclapply(seq(nrow(rho_ind)), function(i) {
				l	<- rho_ind[i,1] ;
				s	<- rho_ind[i,2] ;
				
				return((MGM_1$rho[l,s] - MGM_2$rho[l,s])^2) ;
			}, mc.cores=cores))) ;
			num_param	<- num_param + nrow(rho_ind) ;
			
			rm(rho_ind) ;	gc() ;
		}
		
		# Phi
		if (p_y > 0) {
			## Diagonal
			sumdiff	<- sumdiff + sum(unlist(mclapply(seq(L_sum), function(r) {
				return((MGM_1$phi[r,r] - MGM_2$phi[r,r])^2) ;
			}, mc.cores=cores))) ;
			num_param	<- num_param + L_sum ;
			
			# Non-diagonal
			if (p_y > 1) {
				phi_ind		<- combn(p_y, 2) ;
			
				sumdiff	<- sumdiff + sum(unlist(mclapply(seq(ncol(phi_ind)), function(j) {
					r1	<- phi_ind[1,j] ;
					r2	<- phi_ind[2,j] ;
					
					ind_tmp1	<- (L_cumsum[r1]+1):(L_cumsum[r1+1]) ;
					ind_tmp2	<- (L_cumsum[r2]+1):(L_cumsum[r2+1]) ;
					
					return(sum((MGM_1$phi[ind_tmp1, ind_tmp2] - MGM_2$phi[ind_tmp1, ind_tmp2])^2)) ;
				}, mc.cores=cores))) ;
			
				rm(phi_ind) ;	gc() ;
				num_param	<- num_param + (L_sum^2 - sum(MGM_1$Levels^2)) / 2 ;
			}
		}	
	}
	
	return(list(sumSq=sumdiff, num=num_param)) ;
}

# Return sum of d_deriv*(After-Deriv)
deltaD		<- function(MGM_deriv, MGM_after, cores) {

	X		<- MGM_after[[1]]$Continuous ;
	Y		<- MGM_after[[1]]$Discrete ;
	Y_dummy	<- MGM_after[[1]]$Dummy ;

	p_x	<- ncol(X) ;	p_y	<- ncol(Y) ;
	L_sum		<- ncol(Y_dummy) ;
	L_cumsum	<- c(0, cumsum(MGM_after[[1]]$Levels)) ;
	
	rm(X, Y, Y_dummy) ;

	deriv_diff	<- 0 ;
	
	G	<- length(MGM_after) ;
	
	for (g in seq(G)) {
		
		# Alpha
		if (p_x > 0) {
			deriv_diff	<- deriv_diff + sum(unlist(mclapply(seq(p_x), function(s) {
				deriv_tmp	<- MGM_deriv[[g]]$d_alpha[s,1] ;
				diff_tmp	<- MGM_after[[g]]$alpha[s,1] - MGM_deriv[[g]]$alpha[s,1] ;
				return(deriv_tmp*diff_tmp) ;
			}, mc.cores=cores))) ;
		}

		# Beta
		if (p_x > 0) {
			## Diagonal
			deriv_diff	<- deriv_diff + sum(unlist(mclapply(seq(p_x), function(s) {
				deriv_tmp	<- MGM_deriv[[g]]$d_beta[s,s] ;
				diff_tmp	<- MGM_after[[g]]$beta[s,s] - MGM_deriv[[g]]$beta[s,s] ;
				return(deriv_tmp*diff_tmp) ;
			}, mc.cores=cores))) ;
			
			## Non-diagonal
			if (p_x > 1) {
				beta_ind	<- combn(p_x,2) ;
				deriv_diff	<- deriv_diff + sum(unlist(mclapply(seq(ncol(beta_ind)), function(j) {
					s1	<- beta_ind[1,j] ;
					s2	<- beta_ind[2,j] ;
					
					deriv_tmp	<- MGM_deriv[[g]]$d_beta[s1,s2] ;
					diff_tmp	<- MGM_after[[g]]$beta[s1,s2] - MGM_deriv[[g]]$beta[s1,s2] ;
					
					return(deriv_tmp*diff_tmp) ;
				}, mc.cores=cores))) ;
			
				rm(beta_ind) ;	gc() ;
			}
		}
		
		# Rho
		if (p_x > 0 & p_y > 0) {
			rho_ind		<- expand.grid(seq(L_sum), seq(p_x)) ;
			deriv_diff	<- deriv_diff + sum(unlist(mclapply(seq(nrow(rho_ind)), function(i) {
				l	<- rho_ind[i,1] ;
				s	<- rho_ind[i,2] ;
				
				deriv_tmp	<- MGM_deriv[[g]]$d_rho[l,s] ;
				diff_tmp	<- MGM_after[[g]]$rho[l,s] - MGM_deriv[[g]]$rho[l,s] ;
				
				return(deriv_tmp*diff_tmp) ;
			}, mc.cores=cores))) ;
			
			rm(rho_ind) ;	gc() ;
		}
		
		# Phi
		if (p_y > 0) {
			## Diagonal
			deriv_diff	<- deriv_diff + sum(unlist(mclapply(seq(L_sum), function(r) {
				deriv_tmp	<- MGM_deriv[[g]]$d_phi[r,r] ;
				diff_tmp	<- MGM_after[[g]]$phi[r,r] - MGM_deriv[[g]]$phi[r,r] ;
				return(deriv_tmp*diff_tmp) ;
			}, mc.cores=cores))) ;
			
			# Non-diagonal
			if (p_y > 1) {
				phi_ind		<- combn(p_y, 2) ;
			
				deriv_diff	<- deriv_diff + sum(unlist(mclapply(seq(ncol(phi_ind)), function(j) {
					r1	<- phi_ind[1,j] ;
					r2	<- phi_ind[2,j] ;
					
					ind_tmp1	<- (L_cumsum[r1]+1):(L_cumsum[r1+1]) ;
					ind_tmp2	<- (L_cumsum[r2]+1):(L_cumsum[r2+1]) ;
					
					deriv_tmp	<- MGM_deriv[[g]]$d_phi[ind_tmp1, ind_tmp2] ;
					diff_tmp	<- MGM_after[[g]]$phi[ind_tmp1, ind_tmp2] - MGM_deriv[[g]]$phi[ind_tmp1, ind_tmp2] ;
					
					return(sum(deriv_tmp*diff_tmp)) ;
				}, mc.cores=cores))) ;
			}
		}
	}
	
	return(deriv_diff) ;
}

# Move points in MGM_init & store the results in MGM_get
# Get = Init - coef*d_init
movepoints	<- function(MGM_get, MGM_init, coef, cores) {

	X		<- MGM_get[[1]]$Continuous ;
	Y		<- MGM_get[[1]]$Discrete ;
	Y_dummy	<- MGM_get[[1]]$Dummy ;

	p_x	<- ncol(X) ;	p_y	<- ncol(Y) ;
	L_sum		<- ncol(Y_dummy) ;
	L_cumsum	<- c(0, cumsum(MGM_get[[1]]$Levels)) ;
	
	rm(X, Y, Y_dummy) ;
	gc() ;
	
	G	<- length(MGM_get) ;
	rat_tmp	<- numeric() ;
	
	for (g in seq(G)) {
	
		MGM_g	<- MGM_get[[g]] ;
		MGM_i	<- MGM_init[[g]] ;
	
		# Alpha
		if (p_x > 0) {
			mclapply(seq(p_x), function(s) {
				MGM_g$alpha[s,1]	<- MGM_i$alpha[s,1] - coef*MGM_i$d_alpha[s,1] ;
				return(NULL) ;
			}, mc.cores=cores) ;
		}
		
		# Beta
		if (p_x > 0) {
			beta_ind	<- expand.grid(seq(p_x), seq(p_x)) ;
			mclapply(seq(nrow(beta_ind)), function(i) {
				s1	<- beta_ind[i,1] ;
				s2	<- beta_ind[i,2] ;
				
				MGM_g$beta[s1,s2]	<- MGM_i$beta[s1,s2] - coef*MGM_i$d_beta[s1,s2] ;
				return(NULL) ;
			}, mc.cores=cores) ;
			
			rm(beta_ind) ;	gc() ;
		}
		
		# Rho
		if (p_x > 0 & p_y > 0) {
			rho_ind		<- expand.grid(seq(L_sum), seq(p_x)) ;
			mclapply(seq(nrow(rho_ind)), function(i) {
				l	<- rho_ind[i,1] ;
				s	<- rho_ind[i,2] ;
				
				MGM_g$rho[l,s]		<- MGM_i$rho[l,s] - coef*MGM_i$d_rho[l,s] ;
				return(NULL) ;
			}, mc.cores=cores) ;
			
			rm(rho_ind) ;	gc() ;
		}
		
		# Phi
		if (p_y > 0) {
			mclapply(seq(L_sum), function(r) {
				MGM_g$phi[r,r]		<- MGM_i$phi[r,r] - coef*MGM_i$d_phi[r,r] ;
				return(NULL) ;
			}, mc.cores=cores) ;
			
			# Non-diagonal
			if (p_y > 1) {
				phi_ind		<- combn(p_y, 2) ;
			
				mclapply(seq(ncol(phi_ind)), function(j) {
					r1	<- phi_ind[1,j] ;
					r2	<- phi_ind[2,j] ;
					
					ind_tmp1	<- (L_cumsum[r1]+1):(L_cumsum[r1+1]) ;
					ind_tmp2	<- (L_cumsum[r2]+1):(L_cumsum[r2+1]) ;
					
					MGM_g$phi[ind_tmp1,ind_tmp2]	<- MGM_i$phi[ind_tmp1,ind_tmp2] - coef*MGM_i$d_phi[ind_tmp1,ind_tmp2] ;
					MGM_g$phi[ind_tmp2,ind_tmp1]	<- MGM_i$phi[ind_tmp2,ind_tmp1] - coef*MGM_i$d_phi[ind_tmp2,ind_tmp1] ;
					
					return(NULL) ;
				}, mc.cores=cores) ;
			
				rm(phi_ind) ;	gc() ;
			}
		}
		
		# Vector for skipping some unnecessary iterations
		if (p_x > 0) {
			rat_tmp	<- c(rat_tmp, unlist(mclapply(seq(p_x), function(s) {
				if (MGM_g$beta[s,s] < 0) {
					return(MGM_i$d_beta[s,s]/((1/coef)*diag(MGM_i$beta[s,s]))) ;
				} else {
					return(NULL) ;
				}
			}, mc.cores=cores))) ;
		} else {
			rat_tmp	<- numeric(0) ;
		}
	}
	
	return(rat_tmp) ;
}

# Add 2 MGM lists
# Get = coef1*Add1 + coef2*Add2
add_MGM			<- function(MGM_get, MGM_add1, MGM_add2, coef1, coef2, cores) {

	X		<- MGM_get[[1]]$Continuous ;
	Y		<- MGM_get[[1]]$Discrete ;
	Y_dummy	<- MGM_get[[1]]$Dummy ;

	p_x	<- ncol(X) ;	p_y	<- ncol(Y) ;
	L_sum		<- ncol(Y_dummy) ;
	L_cumsum	<- c(0, cumsum(MGM_get[[1]]$Levels)) ;
	
	rm(X, Y, Y_dummy) ;
	
	G	<- length(MGM_get) ;
	
	for (g in seq(G)) {
	
		MGM_g	<- MGM_get[[g]] ;
		MGM_1	<- MGM_add1[[g]] ;
		MGM_2	<- MGM_add2[[g]] ;
	
		# Alpha
		if (p_x > 0) {
			mclapply(seq(p_x), function(s) {
				MGM_g$alpha[s,1]	<- coef1*MGM_1$alpha[s,1] + coef2*MGM_2$alpha[s,1] ;
				return(NULL) ;
			}, mc.cores=cores) ;
		}
		
		# Beta
		if (p_x > 0) {
			beta_ind	<- expand.grid(seq(p_x), seq(p_x)) ;
			mclapply(seq(nrow(beta_ind)), function(i) {
				s1	<- beta_ind[i,1] ;
				s2	<- beta_ind[i,2] ;
				
				MGM_g$beta[s1,s2]	<- coef1*MGM_1$beta[s1,s2] + coef2*MGM_2$beta[s1,s2] ;
				return(NULL) ;
			}, mc.cores=cores) ;
			
			rm(beta_ind) ;	gc() ;
		}
		
		# Rho
		if (p_x > 0 & p_y > 0) {
			rho_ind		<- expand.grid(seq(L_sum), seq(p_x)) ;
			mclapply(seq(nrow(rho_ind)), function(i) {
				l	<- rho_ind[i,1] ;
				s	<- rho_ind[i,2] ;
				
				MGM_g$rho[l,s]		<- coef1*MGM_1$rho[l,s] + coef2*MGM_2$rho[l,s] ;
				return(NULL) ;
			}, mc.cores=cores) ;
			
			rm(rho_ind) ;	gc() ;
		}
		
		# Phi
		if (p_y > 0) {
			mclapply(seq(L_sum), function(r) {
				MGM_g$phi[r,r]		<- coef1*MGM_1$phi[r,r] + coef2*MGM_2$phi[r,r] ;
				return(NULL) ;
			}, mc.cores=cores) ;
			
			# Non-diagonal
			if (p_y > 1) {
				phi_ind		<- combn(p_y, 2) ;
			
				mclapply(seq(ncol(phi_ind)), function(j) {
					r1	<- phi_ind[1,j] ;
					r2	<- phi_ind[2,j] ;
					
					ind_tmp1	<- (L_cumsum[r1]+1):(L_cumsum[r1+1]) ;
					ind_tmp2	<- (L_cumsum[r2]+1):(L_cumsum[r2+1]) ;
					
					MGM_g$phi[ind_tmp1,ind_tmp2]	<- coef1*MGM_1$phi[ind_tmp1,ind_tmp2] + coef2*MGM_2$phi[ind_tmp1,ind_tmp2] ;
					MGM_g$phi[ind_tmp2,ind_tmp1]	<- coef1*MGM_1$phi[ind_tmp2,ind_tmp1] + coef2*MGM_2$phi[ind_tmp2,ind_tmp1] ;
					
					return(NULL) ;
				}, mc.cores=cores) ;
			
				rm(phi_ind) ;	gc() ;
			}
		}
	}
	
	return(NULL) ;
}