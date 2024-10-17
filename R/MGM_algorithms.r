require(fastDummies) ;
require(Rcpp) ;
require(parallel) ;
require(bigmemory) ;
require(bigalgebra) ;
require(biganalytics) ;

# log of pseudo likelihood + gradient of the smooth function
f_fun	<- function(MGM_input, cores) {

	X		<- MGM_input$Continuous ;
	Y		<- MGM_input$Discrete ;
	Y_dummy	<- MGM_input$Dummy ;

	n_x	<- nrow(X) ;	p_x	<- ncol(X) ;
	n_y	<- nrow(Y) ;	p_y	<- ncol(Y) ;
	L_sum		<- ncol(Y_dummy) ;
	L_cumsum	<- c(0, cumsum(MGM_input$Levels)) ;
	
	L_x		<- MGM_input$Cont_loss ;
	b_y		<- MGM_input$Disc_pred ;
	Prob_y	<- MGM_input$Disc_prob ;

	ll_return	<- 0 ;
		
	alpha	<- MGM_input$alpha ;
	beta	<- MGM_input$beta ;
	rho		<- MGM_input$rho ;
	phi		<- MGM_input$phi ;
	
	d_alpha	<- MGM_input$d_alpha ;
	d_beta	<- MGM_input$d_beta ;
	d_rho	<- MGM_input$d_rho ;
	d_phi	<- MGM_input$d_phi ;
	
	for (s in seq(p_x)) {
		if (beta[s,s] == 0) {
			return(Inf) ;
		}
	}

	# Define linear loss functions of cont. variables
	L_ind	<- expand.grid(seq(n_x), seq(p_x)) ;
	
	mclapply(seq(nrow(L_ind)), function(i) {
		n	<- L_ind[i,1] ;
		s	<- L_ind[i,2] ;
		
		L_x[n,s]	<- sum(Y_dummy[n,]*rho[,s]) - sum(X[n,]*beta[,s]) ;
		L_x[n,s]	<- L_x[n,s] + alpha[s,1] ;
		L_x[n,s]	<- L_x[n,s] / beta[s,s] ;
		
		return(NULL) ;
	}, mc.cores=cores) ;
	
	gc() ;
	
	# Define linear predictors of discrete variables
	b_ind	<- expand.grid(seq(n_y), seq(L_sum)) ;
	
	mclapply(seq(nrow(b_ind)), function(i) {
		n	<- b_ind[i,1] ;
		l	<- b_ind[i,2] ;
		
		b_y[n,l]	<- sum(X[n,]*rho[l,]) + sum(Y_dummy[n,]*phi[,l]) - Y_dummy[n,l]*phi[l,l] ;
		b_y[n,l]	<- b_y[n,l] + phi[l,l] ;
		
		Prob_y[n,l]	<- b_y[n,l] ;
		
		return(NULL) ;
	}, mc.cores=cores) ;
	
	rm(b_ind) ;	gc() ;
	
	# Define expected probabilities for each of the categories
	Prob_ind	<- expand.grid(seq(n_y), seq(p_y)) ;
	
	mclapply(seq(nrow(Prob_ind)), function(i) {
		n	<- Prob_ind[i,1] ;
		r	<- Prob_ind[i,2] ;
		
		ind_tmp	<- (L_cumsum[r]+1):(L_cumsum[r+1]) ;
		denom	<- logsumexp(b_y[n,ind_tmp]) ;
		Prob_y[n,ind_tmp]	<- exp(Prob_y[n,ind_tmp] - denom) ;
		
		return(NULL) ;
	}, mc.cores=cores) ;
	
	gc() ;
	
	# Log-likelihood for continuous variables
	ll_return	<- ll_return + (1/2) * sum(unlist(mclapply(seq(nrow(L_ind)), function(i) {
		n	<- L_ind[i,1] ;
		s	<- L_ind[i,2] ;
		
		ll_tmp1	<- log(beta[s,s]) ;
		ll_tmp2	<- L_x[n,s]^2*beta[s,s] ;
		
		return(-ll_tmp1 + ll_tmp2) ;
	}, mc.cores=cores))) ;
	
	rm(L_ind) ;	gc() ;
	
	# ll_return	<- ll_return - n_x * (1/2) * (sum(log(beta_diag)) - p_x*log(2*pi)) +
	# 							(1/2) * sum((L_x^2) %*% beta_diag) ;
	
	# Log-likelihood for discrete variables
	ll_return	<- ll_return - sum(unlist(mclapply(seq(nrow(Prob_ind)), function(i) {
		n	<- Prob_ind[i,1] ;
		r	<- Prob_ind[i,2] ;
		
		ind_tmp	<- (L_cumsum[r]+1):(L_cumsum[r+1]) ;
		Prob_y_tmp	<- Prob_y[,ind_tmp] ;
		
		return(log(Prob_y_tmp[n,Y[n,r]])) ;
	}, mc.cores=cores))) ;
	
	rm(Prob_ind) ;	gc() ;
	
	# 1st derivative wrt alpha
	if (p_x > 0) {
		mclapply(seq(p_x), function(s) {
			d_alpha[s,1]	<- sum(L_x[,s]) ;
			return(NULL) ;
		}, mc.cores=cores) ;
	}
	
	# 1st derivative wrt beta
	if (p_x > 1) {
		beta_ind	<- combn(p_x,2) ;
		mclapply(seq(ncol(beta_ind)), function(j) {
			s1	<- beta_ind[1,j] ;
			s2	<- beta_ind[2,j] ;
			
			d_beta[s1,s2]	<- d_beta[s2,s1]	<- -(sum(X[,s1]*L_x[,s2]) + sum(X[,s2]*L_x[,s1])) ;
			
			return(NULL) ;
		}, mc.cores=cores) ;
		
		rm(beta_ind) ;	gc() ;
	}
	
	if (p_x > 0) {
		mclapply(seq(p_x), function(s) {
			d_beta[s,s]	<- - n_x / (2*beta[s,s]) -
							sum(L_x[,s]^2) / 2 -
							sum(L_x[,s]*X[,s]) ;
			
			return(NULL) ;
		}, mc.cores=cores) ;
	}
	
	# 1st derivative wrt rho
	if (p_x > 0 & p_y > 0) {
		rho_ind		<- expand.grid(seq(L_sum), seq(p_x)) ;
		mclapply(seq(nrow(rho_ind)), function(i) {
			l	<- rho_ind[i,1] ;
			s	<- rho_ind[i,2] ;
			
			d_rho[l,s]	<- sum(Y_dummy[,l]*L_x[,s] + (Prob_y[,l]-Y_dummy[,l])*X[,s]) ;
			
			return(NULL) ;
		}, mc.cores=cores) ;
		
		rm(rho_ind) ;	gc() ;
		# d_rho	<- t(Y_dummy) %*% L_x + t(Prob_y - Y_dummy) %*% X ;
	}
	
	# 1st derivative wrt phi
	if (p_y > 0) {
		phi_ind		<- combn(L_sum, 2) ;
		mclapply(seq(ncol(phi_ind)), function(j) {
			r1	<- phi_ind[1,j] ;
			r2	<- phi_ind[2,j] ;
		
			d_phi[r1,r2]	<- d_phi[r2,r1]		<- 
				sum(Y_dummy[,r1]*(Prob_y[,r2]-Y_dummy[,r2]) + Y_dummy[,r2]*(Prob_y[,r1]-Y_dummy[,r1])) ;
			
			return(NULL) ;
		}, mc.cores=cores) ;
		
		rm(phi_ind) ;	gc() ;
	
		mclapply(seq(p_y), function(r) {
			ind_tmp	<- (L_cumsum[r]+1):(L_cumsum[r+1]) ;
			d_phi[ind_tmp,ind_tmp]	<- 0 ;
			rm(ind_tmp) ;
			return(NULL) ;
		}, mc.cores=cores) ;
		
		mclapply(seq(L_sum), function(l) {
			d_phi[l,l]	<- sum(Prob_y[,l] - Y_dummy[,l]) ;
			return(NULL) ;
		}, mc.cores=cores) ;
	}
	
	return(ll_return) ;
}

# Penalty function + proximal operator
# Assumes all MGMs have common variables & common levels of discrete variables
# Input
## MGM list
## t: proximal operator parameter
## lambda: penalty parameters
## tol: tolerance of RMSE for each of the parameters
## tol_fpa: tolerance used in fixed point approach
## maxit: maximum number of iterations in FPA
g_fun	<- function(MGM_list, t, lambda_intra, lambda_intra_prior=NULL, lambda_inter,
						with_prior=FALSE, wp=NULL,
						tol=5e-03, tol_fpa=1e-12, maxit=1000000, cores) {
	
	X		<- MGM_list[[1]]$Continuous ;
	Y		<- MGM_list[[1]]$Discrete ;
	Y_dummy	<- MGM_list[[1]]$Dummy ;

	p_x	<- ncol(X) ;	p_y	<- ncol(Y) ;
	L_sum		<- ncol(Y_dummy) ;
	L_cumsum	<- c(0, cumsum(MGM_list[[1]]$Levels)) ;
	
	G	<- length(MGM_list) ;
	if (G >= 2) grpind_pair	<- combn(G, 2) ;
	
	penalty	<- 0 ;
	
	# Beta
	if (p_x > 1) {
		beta_ind	<- combn(p_x,2) ;
		
		penalty	<- penalty + sum(unlist(mclapply(seq(ncol(beta_ind)), function(j) {
		
			s1	<- beta_ind[1,j] ;
			s2	<- beta_ind[2,j] ;
			
			var1	<- colnames(X)[s1] ;
			var2	<- colnames(X)[s2] ;
			
			init_values	<- sapply(MGM_list, function(MGM) return(MGM$beta[s1,s2])) ;
			new_values	<- prev_values	<- init_values ;

			ind_1	<- match(var1, colnames(wp)) ;
			ind_2	<- match(var2, colnames(wp)) ;
			
			if (with_prior) {
				if (wp[ind_1,ind_2] == 1) {
					l_intra_tmp 	<- lambda_intra_prior[1] ;
				} else {
					l_intra_tmp		<- lambda_intra[1] ;
				}
			} else {
				l_intra_tmp	<- lambda_intra[1] ;
			}
				
			# l_intra_tmp		<- if (with_prior & wp[var1,var2] == 1) lambda_intra_prior[1] else lambda_intra[1]
		
			while (TRUE) {
				
				# Update for each of the groups
				for (g in seq(G)) {
					prev_val		<- init_values[g] ;		# Previous parameter in g-th group
					xvals_tmp		<- new_values ;
					xvals_tmp[g]	<- 0 ;					# Non-smooth
					
					lambda_tmp		<- rep(lambda_inter[1], G) ;
					lambda_tmp[g]	<- l_intra_tmp ;	# Coefficients
					
					# Target function to minimize
					target_ftn	<- function(x) {
						tg	<- (1/(2*t)) * (x - prev_val)^2 ;
						tg	<- tg + sum(lambda_tmp * abs(x - xvals_tmp)) ;
						
						return(tg) ;
					}
					
					fvals_tmp	<- sapply(xvals_tmp, target_ftn) ;	# Function values for each of the non-smooth points
					ind_min		<- which(fvals_tmp == min(fvals_tmp)) ;
					xval_min	<- xvals_tmp[ind_min[1]] ;			# Corresponding x-values
					
					# Left/right 1st derivatives for the point
					slope_tmp	<- t*sum(lambda_tmp*sign(xval_min - xvals_tmp)) + (xval_min - prev_val) ;
					slope_l		<- slope_tmp - t*sum(lambda_tmp[ind_min])
					slope_r		<- slope_tmp + t*sum(lambda_tmp[ind_min])
					
					if (slope_l * slope_r <= 0) {	# Different signs: minimum
						new_values[g]	<- xval_min ;
					} else if (slope_l > 0) {		# Increasing at the point 
						new_values[g]	<- xval_min - slope_l ;
					} else {						# Descreasing at the point
						new_values[g]	<- xval_min - slope_r ;
					}
					
					# Remove used variables to prevent confusion
					rm(prev_val, xvals_tmp, lambda_tmp, target_ftn, fvals_tmp, ind_min, xval_min, slope_tmp, slope_l, slope_r) ;
				}
				
				if (G >= 2) {
					# For diagonal case, perform additional optimization
					if (new_values[1] == new_values[2]) {
						if (2*l_intra_tmp >= (1/t)*abs(sum(init_values))) {
							new_values	<- rep(0,G) ;
						} else {
							mean_tmp	<- mean(init_values) ;
							new_values	<- rep(mean_tmp - (t*l_intra_tmp)*sign(mean_tmp), 2) ;
						}
					}
				}
				
				# If the condition was satisfied, update the parameters
				if (sqrt(mean((prev_values - new_values)^2)) < tol) {
#					print(new_values) ;
					for (g in seq(G)) {
						MGM_list[[g]]$beta[s1,s2]	<-	MGM_list[[g]]$beta[s2,s1]	<- new_values[g] ;
					}
					break ;
				} else {	# Otherwise, put updated values into previous variable
					prev_values	<- new_values ;
				}
			}
			
			pen_tmp	<- l_intra_tmp * sum(abs(new_values)) ;
			
			if (G >= 2) {
				for (i in seq(ncol(grpind_pair))) {
					pen_tmp	<- pen_tmp + lambda_inter[1] * abs(new_values[grpind_pair[1,i]] - new_values[grpind_pair[2,i]]) ;
				}
			}

			# Remove used variables to prevent confusion
			rm(prev_values, new_values) ;
			gc() ;
			
			return(pen_tmp) ;
		}, mc.cores=cores))) ;
	
		rm(beta_ind) ;	gc() ;
	}
	
	# Rho
	if (p_x > 0 & p_y > 0) {
		rho_ind		<- expand.grid(seq(p_y), seq(p_x)) ;
		
		penalty	<- penalty + sum(unlist(mclapply(seq(nrow(rho_ind)), function(i) {
			r	<- rho_ind[i,1] ;
			s	<- rho_ind[i,2] ;

			var1	<- colnames(Y)[r] ;
			var2	<- colnames(X)[s] ;

			# Initialize parameters with previous values
			init_values		<- lapply(MGM_list, function(MGM) return(MGM$rho[(L_cumsum[r]+1):L_cumsum[r+1],s])) ;
			new_values		<- prev_values	<- init_values ;
			
			ind_1	<- match(var1, colnames(wp)) ;
			ind_2	<- match(var2, colnames(wp)) ;

			if (with_prior) {
				if (wp[ind_1,ind_2] == 1) {
					l_intra_tmp 	<- lambda_intra_prior[2] ;
				} else {
					l_intra_tmp		<- lambda_intra[2] ;
				}
			} else {
				l_intra_tmp	<- lambda_intra[2] ;
			}
			
			# l_intra_tmp		<- if (with_prior & wp[var1,var2] == 1) lambda_intra_prior[2] else lambda_intra[2]
		
			while (TRUE) {
			
				# Update for each of the groups
				for (g in seq(G)) {
					prev_val		<- init_values[[g]] ;					# Previous parameter in g-th group
					xvals_tmp		<- new_values ;
					
					xvals_tmp[[g]]	<- rep(0, length(xvals_tmp[[g]])) ;		# Non-smooth (destination) points
					names(xvals_tmp[[g]])	<- names(new_values[[g]]) ;
					
					lambda_tmp		<- rep(lambda_inter[2], G) ;
					lambda_tmp[g]	<- l_intra_tmp ;					# Coefficients
					
					# Target function to minimize
					target_ftn	<- function(x) {
						tg	<- (1/(2*t)) * sum((x - prev_val)^2) ;
						for (g2 in seq(G)) {
							tg	<- tg + lambda_tmp[g2] * sqrt(sum((x - xvals_tmp[[g2]])^2)) ;
						}
						
						return(tg) ;
					}
					
					# Destination points: criterion according to Katz & Vogl (2010)
					conv_flag	<- FALSE ;
					for (g2 in seq(G)) {
						slope_tmp	<- prev_val - xvals_tmp[[g2]] ;
						for (g3 in seq(G)) {
							if ((g3 != g2) & !isTRUE(all.equal(xvals_tmp[[g2]], xvals_tmp[[g3]]))) {
								diff_tmp	<- xvals_tmp[[g3]] - xvals_tmp[[g2]] ;
								slope_tmp	<- slope_tmp + ((t*lambda_tmp[g3])*(diff_tmp) / sqrt(sum(diff_tmp^2))) ;
							}
						}
						
						if (sum(slope_tmp^2) <= (t^2 * lambda_tmp[g2]^2)) {
							new_values[[g]]	<- xvals_tmp[[g2]] ;
							conv_flag	<- TRUE ;
							break ;
						}
					}
					
					# Non-destination points: fixed-point approach
					if (isFALSE(conv_flag)) {
						if (all(sapply(xvals_tmp, FUN = identical, xvals_tmp[[1]]))) {
							init_point	<- rnorm(length(prev_val), mean=prev_val, sd=.1) ;
						} else {
							init_point	<- lambda_tmp[1] * xvals_tmp[[1]] ;
							init_denom	<- lambda_tmp[1] ;
							for (g2 in seq(2, G)) {
								init_point	<- init_point + (lambda_tmp[g2] * xvals_tmp[[g2]]) ;
								init_denom	<- init_denom + lambda_tmp[g2] ;
							}
							if (init_denom == 0) {
								init_point	<- rowMeans(do.call(cbind, xvals_tmp)) ;
							} else {
								init_point	<- init_point / init_denom ;
								if (isTRUE(sapply(xvals_tmp, all.equal, init_point)) | "TRUE" %in% sapply(xvals_tmp, all.equal, init_point)) {
									init_point	<- rowMeans(do.call(cbind, xvals_tmp)) ;
								}
							}
						}
						
						input_ftn	<- function(x) {
							numer 	<- prev_val ;
							denom	<- 1 ;
							
							for (g2 in seq(G)) {
								diff_tmp	<- (x - xvals_tmp[[g2]]) ;
								numer		<- numer + ((t * lambda_tmp[g2] / sqrt(sum(diff_tmp^2))) * xvals_tmp[[g2]]) ;
								denom		<- denom + ((t * lambda_tmp[g2] / sqrt(sum(diff_tmp^2)))) ;
							}
							
							return((1/denom)*numer) ;
						}
						
						fpa_res	<- fpa(init_point, input_ftn, target_ftn, maxit, tol_fpa) ;

						if (!isFALSE(fpa_res$converged)) {
							new_values[[g]]	<- fpa_res$value ;
						} else {
							stop("Fixed-point method did not converge") ;
						}
					}
					
					# Remove used variables to prevent confusion
					if (isFALSE(conv_flag)) rm(init_point, input_ftn, fpa_res) ;
					rm(prev_val, xvals_tmp, lambda_tmp, target_ftn, conv_flag, slope_tmp) ;
				}

				if (G >= 2) {
					# For diagonal case, perform additional optimization
					if (isTRUE(all.equal(new_values[[1]], new_values[[2]]))) {
						sum_tmp	<- rowSums(do.call(cbind, init_values)) ;
						if (2*l_intra_tmp >= (1/t)*sqrt(sum(sum_tmp^2))) {
							new_values	<- rep(list(rep(0,length(init_values[[1]]))), G) ;
						} else {
							mean_tmp	<- rowMeans(do.call(cbind, init_values)) ;
							init_point	<- mean_tmp ;
							
							target_ftn	<- function(x) {
								tg	<- 2*l_intra_tmp*sqrt(sum(x^2)) ;
								for (g2 in seq(G)) {
									tg	<- tg + (1/(2*t)) * sum((x - init_values[[g2]])^2) ;
								}
								
								return(tg) ;
							}
							
							input_ftn	<- function(x) {
								numer	<- mean_tmp ;
								denom	<- 1 + (t*l_intra_tmp) / sqrt(sum(x^2)) ;
								
								return((1/denom)*numer) ;
							}
							
							fpa_res	<- fpa(init_point, input_ftn, target_ftn, maxit, tol_fpa) ;
							new_values	<- rep(list(fpa_res$value), 2) ;
						}
					}
				}
				
				# If the condition was satisfied, update the parameters
				diff_vec	<- prev_values[[1]] - new_values[[1]] ;
				if (G >= 2) {
					for (g in seq(2, G)) {
						diff_vec	<- c(diff_vec, prev_values[[g]] - new_values[[g]]) ;
					}
				}
				if (sqrt(mean(diff_vec^2)) < tol) {
					for (g in seq(G)) {
						MGM_list[[g]]$rho[(L_cumsum[r]+1):L_cumsum[r+1],s]	<- new_values[[g]] ;
					}
					break ;
				} else {	# Otherwise, put updated values into previous variable
					prev_values	<- new_values ;
				}
			}
				
			pen_tmp	<- 0 ;
			for (g1 in seq(G)) {
				pen_tmp	<- pen_tmp + l_intra_tmp * sqrt(sum(new_values[[g1]]^2)) ;
			}
			
			if (G >= 2) {
				for (g2 in seq(ncol(grpind_pair))) {
					pen_tmp	<- pen_tmp + lambda_inter[2] * sqrt(sum((new_values[[grpind_pair[1,g2]]] - new_values[[grpind_pair[2,g2]]])^2)) ;
				}
			}
			
			# Remove used variables to prevent confusion
			rm(prev_values, new_values) ;
			gc() ;
			
			return(pen_tmp) ;
		}, mc.cores=cores))) ;
		
		rm(rho_ind) ;	gc() ;
	}
	
	# Phi: only if there are 2 or more discrete variables
	if (p_y > 1) {
		phi_ind	<- combn(p_y,2) ;
		
		penalty	<- penalty + sum(unlist(mclapply(seq(ncol(phi_ind)), function(j) {
		
			s1	<- phi_ind[1,j] ;
			s2	<- phi_ind[2,j] ;
			
			var1	<- colnames(Y)[s1] ;
			var2	<- colnames(Y)[s2] ;
			
			# Initialize parameters with previous values
			init_values		<- lapply(MGM_list, function(MGM) return(as.vector(MGM$phi[(L_cumsum[s1]+1):L_cumsum[s1+1],(L_cumsum[s2]+1):L_cumsum[s2+1]]))) ;
			new_values		<- prev_values	<- init_values ;
			
			ind_1	<- match(var1, colnames(wp)) ;
			ind_2	<- match(var2, colnames(wp)) ;

			if (with_prior) {
				if (wp[ind_1,ind_2] == 1) {
					l_intra_tmp 	<- lambda_intra_prior[3] ;
				} else {
					l_intra_tmp		<- lambda_intra[3] ;
				}
			} else {
				l_intra_tmp	<- lambda_intra[3] ;
			}
			
			# l_intra_tmp		<- if (with_prior & wp[var1,var2] == 1) lambda_intra_prior[3] else lambda_intra[3]
			
			while (TRUE) {
			
				# Update for each of the groups
				for (g in seq(G)) {
					prev_val		<- init_values[[g]] ;					# Previous parameter in g-th group

					xvals_tmp		<- new_values ;
					xvals_tmp[[g]]	<- rep(0, length(xvals_tmp[[g]])) ;		# Non-smooth (destination) points
					
					lambda_tmp		<- rep(lambda_inter[3], G) ;
					lambda_tmp[g]	<- l_intra_tmp ;					# Coefficients
					
					# Target function to minimize
					target_ftn	<- function(x) {
						tg	<- (1/(2*t)) * sum((x - prev_val)^2) ;
						for (g2 in seq(G)) {
							tg	<- tg + lambda_tmp[g2] * sqrt(sum((x - xvals_tmp[[g2]])^2)) ;
						}
						
						return(tg) ;
					}
					
					# Destination points: criterion according to Katz & Vogl (2010)
					conv_flag	<- FALSE ;
					for (g2 in seq(G)) {
						slope_tmp	<- prev_val - xvals_tmp[[g2]] ;
						for (g3 in seq(G)) {
							if ((g3 != g2) & !isTRUE(all.equal(xvals_tmp[[g2]], xvals_tmp[[g3]]))) {
								diff_tmp	<- xvals_tmp[[g3]] - xvals_tmp[[g2]] ;
								slope_tmp	<- slope_tmp + ((t*lambda_tmp[g3])*(diff_tmp) / sqrt(sum(diff_tmp^2))) ;
							}
						}
						
						if (sum(slope_tmp^2) <= (t^2 * lambda_tmp[g2]^2)) {
							new_values[[g]]	<- xvals_tmp[[g2]] ;
							conv_flag	<- TRUE ;
							break ;
						}
					}
					
					# Non-destination points: fixed-point approach
					if (isFALSE(conv_flag)) {
						if (all(sapply(xvals_tmp, FUN = identical, xvals_tmp[[1]]))) {
							init_point	<- rnorm(length(prev_val), mean=prev_val, sd=.1) ;
						} else {
							init_point	<- lambda_tmp[1] * xvals_tmp[[1]] ;
							init_denom	<- lambda_tmp[1] ;
							for (g2 in seq(2, G)) {
								init_point	<- init_point + (lambda_tmp[g2] * xvals_tmp[[g2]]) ;
								init_denom	<- init_denom + lambda_tmp[g2] ;
							}
							if (init_denom == 0) {
								init_point	<- rowMeans(do.call(cbind, xvals_tmp)) ;
							} else {
								init_point	<- init_point / init_denom ;
								if (isTRUE(sapply(xvals_tmp, all.equal, init_point)) | "TRUE" %in% sapply(xvals_tmp, all.equal, init_point)) {
									init_point	<- rowMeans(do.call(cbind, xvals_tmp)) ;
								}
							}
						}
						
						input_ftn	<- function(x) {
							numer 	<- prev_val ;
							denom	<- 1 ;
							
							for (g2 in seq(G)) {
								diff_tmp	<- (x - xvals_tmp[[g2]]) ;
								numer		<- numer + ((t * lambda_tmp[g2] / sqrt(sum(diff_tmp^2))) * xvals_tmp[[g2]]) ;
								denom		<- denom + ((t * lambda_tmp[g2] / sqrt(sum(diff_tmp^2)))) ;
							}
							
							return((1/denom)*numer) ;
						}
						
						fpa_res	<- fpa(init_point, input_ftn, target_ftn, maxit, tol_fpa) ;
						
						if (!isFALSE(fpa_res$converged)) {
							new_values[[g]]	<- fpa_res$value ;
						} else {
							stop("Fixed-point method did not converge") ;
						}
					}
					
					# Remove used variables to prevent confusion
					if (isFALSE(conv_flag)) rm(init_point, input_ftn) ;
					rm(prev_val, xvals_tmp, lambda_tmp, target_ftn, conv_flag, slope_tmp) ;
				}
				
				if (G >= 2) {	
					# For diagonal case, perform additional optimization
					if (isTRUE(all.equal(new_values[[1]], new_values[[2]], tol=0))) {
						sum_tmp	<- rowSums(do.call(cbind, init_values)) ;
						if (2*l_intra_tmp >= (1/t)*sqrt(sum(sum_tmp^2))) {
							new_values	<- rep(list(rep(0,length(init_values[[1]]))),G) ;
						} else {
							mean_tmp	<- rowMeans(do.call(cbind, init_values)) ;
							init_point	<- mean_tmp ;
							
							target_ftn	<- function(x) {
								tg	<- 2*l_intra_tmp*sqrt(sum(x^2)) ;
								for (g2 in seq(G)) {
									tg	<- tg + (1/(2*t)) * sum((x - init_values[[g2]])^2) ;
								}
								
								return(tg) ;
							}
							
							input_ftn	<- function(x) {
								numer	<- mean_tmp ;
								denom	<- 1 + (t*l_intra_tmp) / sqrt(sum(x^2)) ;
								
								return((1/denom)*numer) ;
							}
							
							fpa_res	<- fpa(init_point, input_ftn, target_ftn, maxit, tol_fpa) ;
							new_values	<- rep(list(fpa_res$value), 2) ;
						}
					}
				}
				
				# If the condition was satisfied, update the parameters
				diff_vec	<- prev_values[[1]] - new_values[[1]] ;
				if (G >= 2) {
					for (g in seq(2, G)) {
						diff_vec	<- c(diff_vec, prev_values[[g]] - new_values[[g]]) ;
					}
				}
				if (sqrt(mean(diff_vec^2)) < tol) {
					for (g in seq(G)) {
						values_put	<- matrix(new_values[[g]], nrow=(L_cumsum[s1+1]-L_cumsum[s1])) ;
						MGM_list[[g]]$phi[(L_cumsum[s1]+1):L_cumsum[s1+1],(L_cumsum[s2]+1):L_cumsum[s2+1]]	<- values_put ;
						MGM_list[[g]]$phi[(L_cumsum[s2]+1):L_cumsum[s2+1],(L_cumsum[s1]+1):L_cumsum[s1+1]]	<- t(values_put) ;
						rm(values_put) ;
					}
					break ;
				} else {	# Otherwise, put updated values into previous variable
					prev_values	<- new_values ;
				}
			}
			
			pen_tmp 	<- 0 ;

			for (g1 in seq(G)) {
				pen_tmp	<- pen_tmp + l_intra_tmp * sqrt(sum(new_values[[g1]]^2)) ;
			}
			
			if (G >= 2) {
				for (g2 in seq(ncol(grpind_pair))) {
					pen_tmp	<- pen_tmp + lambda_inter[3] * sqrt(sum((new_values[[grpind_pair[1,g2]]] - new_values[[grpind_pair[2,g2]]])^2)) ;
				}
			}
			
			# Remove used variables to prevent confusion
			rm(prev_values, new_values) ;
			gc() ;
			
			return(pen_tmp) ;
		}, mc.cores=cores))) ;
		
		rm(phi_ind) ;	gc() ;
	}
	
	return(penalty) ;
}

# Penalty function only
# Input
## MGM list
## lambda: penalty parameters
g_pen	<- function(MGM_list, lambda_intra, lambda_intra_prior=NULL, lambda_inter, 
						with_prior=FALSE, wp=NULL, cores) {
	
	X		<- MGM_list[[1]]$Continuous ;
	Y		<- MGM_list[[1]]$Discrete ;
	Y_dummy	<- MGM_list[[1]]$Dummy ;

	p_x	<- ncol(X) ;	p_y	<- ncol(Y) ;
	L_sum		<- ncol(Y_dummy) ;
	L_cumsum	<- c(0, cumsum(MGM_list[[1]]$Levels)) ;
	
	G	<- length(MGM_list) ;
	if (G >= 2) grpind_pair	<- combn(G, 2) ;
	
	penalty	<- 0 ;
	
	# Beta
	if (p_x > 1) {
		beta_ind	<- combn(p_x,2) ;
		
		penalty	<- penalty + sum(unlist(mclapply(seq(ncol(beta_ind)), function(j) {
		
			s1	<- beta_ind[1,j] ;
			s2	<- beta_ind[2,j] ;
			
			var1	<- colnames(X)[s1] ;
			var2	<- colnames(X)[s2] ;
			
			param_values	<- sapply(MGM_list, function(MGM) return(MGM$beta[s1,s2])) ;
			
			ind_1	<- match(var1, colnames(wp)) ;
			ind_2	<- match(var2, colnames(wp)) ;

			if (with_prior) {
				if (wp[ind_1,ind_2] == 1) {
					l_intra_tmp 	<- lambda_intra_prior[1] ;
				} else {
					l_intra_tmp		<- lambda_intra[1] ;
				}
			} else {
				l_intra_tmp	<- lambda_intra[1] ;
			}
			
			# l_intra_tmp		<- if (with_prior & wp[var1,var2] == 1) lambda_intra_prior[1] else lambda_intra[1]
		
			pen_tmp	<- l_intra_tmp * sum(abs(param_values)) ;
			if (G >= 2) {
				for (i in seq(ncol(grpind_pair))) {
					pen_tmp	<- pen_tmp + lambda_inter[1] * abs(param_values[grpind_pair[1,i]] - param_values[grpind_pair[2,i]]) ;
				}
			}
			
			# Remove used variables to prevent confusion
			rm(param_values) ;
			gc() ;
			
			return(pen_tmp) ;
		}, mc.cores=cores))) ;
		
		rm(beta_ind) ;	gc() ;
	}
	
	# Rho
	if (p_x > 0 & p_y > 0) {
		rho_ind		<- expand.grid(seq(p_y), seq(p_x)) ;
		
		penalty	<- penalty + sum(unlist(mclapply(seq(nrow(rho_ind)), function(i) {
			r	<- rho_ind[i,1] ;
			s	<- rho_ind[i,2] ;

			var1	<- colnames(Y)[r] ;
			var2	<- colnames(X)[s] ;
			
			param_values	<- lapply(MGM_list, function(MGM) return(MGM$rho[(L_cumsum[r]+1):L_cumsum[r+1],s])) ;
			
			ind_1	<- match(var1, colnames(wp)) ;
			ind_2	<- match(var2, colnames(wp)) ;

			if (with_prior) {
				if (wp[ind_1,ind_2] == 1) {
					l_intra_tmp 	<- lambda_intra_prior[2] ;
				} else {
					l_intra_tmp		<- lambda_intra[2] ;
				}
			} else {
				l_intra_tmp	<- lambda_intra[2] ;
			}
			
			# l_intra_tmp		<- if (with_prior & wp[var1,var2] == 1) lambda_intra_prior[2] else lambda_intra[2]
			
			pen_tmp	<- 0 ;

			for (g1 in seq(G)) {
				pen_tmp	<- pen_tmp + l_intra_tmp * sqrt(sum(param_values[[g1]]^2)) ;
			}
			if (G >= 2) {
				for (g2 in seq(ncol(grpind_pair))) {
					pen_tmp	<- pen_tmp + lambda_inter[2] * sqrt(sum((param_values[[grpind_pair[1,g2]]] - param_values[[grpind_pair[2,g2]]])^2)) ;
				}
			}

			# Remove used variables to prevent confusion
			rm(param_values) ;
			gc() ;
			
			return(pen_tmp) ;
		}, mc.cores=cores))) ;
	}
	
	rm(rho_ind) ;	gc() ;
	
	# Phi: only if there are 2 or more discrete variables
	if (p_y > 1) {
		phi_ind	<- combn(p_y,2) ;
		
		penalty	<- penalty + sum(unlist(mclapply(seq(ncol(phi_ind)), function(j) {
		
			s1	<- phi_ind[1,j] ;
			s2	<- phi_ind[2,j] ;

			var1	<- colnames(Y)[s1] ;
			var2	<- colnames(Y)[s2] ;
			
			param_values	<- lapply(MGM_list, function(MGM) return(as.vector(MGM$phi[(L_cumsum[s1]+1):L_cumsum[s1+1],(L_cumsum[s2]+1):L_cumsum[s2+1]]))) ;
			
			ind_1	<- match(var1, colnames(wp)) ;
			ind_2	<- match(var2, colnames(wp)) ;

			if (with_prior) {
				if (wp[ind_1,ind_2] == 1) {
					l_intra_tmp 	<- lambda_intra_prior[3] ;
				} else {
					l_intra_tmp		<- lambda_intra[3] ;
				}
			} else {
				l_intra_tmp	<- lambda_intra[3] ;
			}
			
			# l_intra_tmp		<- if (with_prior & wp[var1,var2] == 1) lambda_intra_prior[3] else lambda_intra[3]
			
			pen_tmp 	<- 0 ;
			for (g1 in seq(G)) {
				pen_tmp	<- pen_tmp + l_intra_tmp * sqrt(sum(param_values[[g1]]^2)) ;
			}
			if (G >= 2) {
				for (g2 in seq(ncol(grpind_pair))) {
					pen_tmp	<- pen_tmp + lambda_inter[3] * sqrt(sum((param_values[[grpind_pair[1,g2]]] - param_values[[grpind_pair[2,g2]]])^2)) ;
				}
			}
			
			# Remove used variables to prevent confusion
			rm(param_values) ;
			gc() ;
			
			return(pen_tmp) ;
		}, mc.cores=cores))) ;
		
		rm(phi_ind) ;	gc() ;
	}
	
	return(penalty) ;
}

# Learn with PGM method
# Input
## MGM list & MGM for backtracking
## t: proximal operator parameter
## L: step size
## eta: step size magnifier
## lambda: penalty parameters
## tol_mgm: tolerance for the convergence in the main iteration
## tol_g: tolerance for the convergence of the proximal operator
## tol_fpa: tolerance for the fixed point approach
## maxit: maximum number of iterations in FPA
## cores
## verbose
learn_tb	<- function(MGM_list, t=1, L=NULL, eta=2, lambda_intra, lambda_intra_prior=NULL, lambda_inter,
			       			with_prior=FALSE, wp=NULL, converge_by_edge=TRUE, tol_edge=3,
						tol_mgm=1e-04, tol_g=5e-03, tol_fpa=1e-12, tol_polish=1e-12, maxit=1000000, cores, verbose=FALSE) {
	L_miss	<- missing(L) || is.null(L) ;
	G	<- length(MGM_list) ;
	n_tot	<- sum(sapply(MGM_list, function(MGM) nrow(MGM$Continuous)))
	
	lambda_intra	<- n_tot * lambda_intra ;
	lambda_inter	<- n_tot * lambda_inter ;
	if (!is.null(lambda_intra_prior)) {
		lambda_intra_prior	<- n_tot * lambda_intra_prior ;
	}
	
	MGM_x_old	<- deepcopy.MGMlist(MGM_list) ;
	MGM_x_new	<- deepcopy.MGMlist(MGM_list) ;
	MGM_y_old	<- deepcopy.MGMlist(MGM_list) ;
	MGM_y_new	<- deepcopy.MGMlist(MGM_list) ;
	MGM_z		<- deepcopy.MGMlist(MGM_list) ;
	t_old	<- t_new	<- t ;
	
	ll_x_old	<- pen_x_old	<- Inf ;
	target_old	<- Inf ;
	edge_diff	<- 0 ;

	if (is.null(L)) {
		L	<- n_tot ;
	} else {
		L	<- L*n_tot ;
	}
	
	while (TRUE) {

		L	<- L*.9 ;

		# L not given: use backtracking method
		if (L_miss) {
			if (verbose) print("Backtracking started") ;
		
			# if (is.null(L))	L	<- n_tot ;
			
			if (verbose) print("   Calculating the gradients & moving the parameters") ;
			ll_y_old	<- 0 ;
			for (g in seq(G)) {
				ll_tmp		<- f_fun(MGM_y_old[[g]], cores) ;	
				ll_y_old	<- ll_y_old + ll_tmp ;
			}
			
			while (TRUE) {
				
				n_x	<- nrow(MGM_z[[1]]$Continuous) ;
				p_x	<- ncol(MGM_z[[1]]$Continuous) ;
				
				n_y		<- nrow(MGM_z[[1]]$Discrete) ;
				p_y		<- ncol(MGM_z[[1]]$Discrete) ;
				L_sum	<- ncol(MGM_z[[1]]$Dummy) ;
				
				# Calculate gradients & move parameters
				rat_tmp	<- movepoints(MGM_z, MGM_y_old, 1/L, cores) ;

				if (length(rat_tmp) > 0) {
					k	<- ceiling(max(log(rat_tmp, base=eta))) ;
					L	<- L*(eta^k) ;
					if (verbose) print(paste0("   Skipping ", k, " steps")) ;
					next ;
				}
				
				rm(rat_tmp) ;

				# Approximate wrt penalty functions
				if (verbose) print("   Calculating approximate optimals") ;
				gg_z	<- g_fun(MGM_z, 1/L, lambda_intra, lambda_intra_prior, lambda_inter, 
									with_prior, wp, tol_g, tol_fpa, maxit, cores) ;
				ll_z	<- 0 ;
				
				# Re-calculate log-likelihoods after the approximation
				if (verbose) print("   Re-calculating the log-likelihoods after the approximation") ;
				for (g in seq(G)) {
					ll_tmp		<- f_fun(MGM_z[[g]], cores) ;
					ll_z		<- ll_z + ll_tmp ;
				}
				
				sq_diff		<- rmsd_MGM(MGM_y_old, MGM_z, cores) ;
				deriv_diff	<- deltaD(MGM_y_old, MGM_z, cores) ;
				
				# Target function
				FF	<- ll_z + gg_z ;
				# Approximated function values wrt 1st derivative
				QQ	<- ll_y_old + deriv_diff + (L/2)*sq_diff$sumSq + gg_z  ;

				if (FF - QQ > 0) {
					copy.MGMlist(MGM_y_old, MGM_z, cores) ;
					L	<- L*eta ;
					if (verbose) print(paste0("   Updating the step size, currently: ", L)) ;
				} else {
					if (verbose) print(paste0("Backtracking ended, step size: ", L)) ;
					break ;
				}
			}
		}
		
		if (!L_miss) {
			for (g in seq(G)) {
				f_fun(MGM_y_old[[g]], cores) ;	
			}
		}
		
		# movepoints(MGM_z, MGM_y_old, 1/L, cores) ;
		
		# if (verbose) print("   Calculating approximate optimals") ;
		# gg_tmp		<- g_fun(MGM_z, 1/L, lambda_intra, lambda_inter, sparse, tol_g, tol_fpa, maxit, cores) ;
		# pen_x_old	<- g_pen(MGM_x_old, lambda_intra, lambda_inter, sparse, cores) ;
		
		# ll_z	<- ll_x_old	<- 0 ;
		# for (g in seq(G)) {
		# 	ll_z	<- ll_z + f_fun(MGM_z[[g]], cores) ;
		# 	ll_x_old	<- ll_x_old + f_fun(MGM_x_old[[g]], cores) ;
		# }
		
		F_z		<- gg_z + ll_z ;
		F_x_old	<- pen_x_old + ll_x_old ;
		
		if (F_z <= F_x_old) {
			copy.MGMlist(MGM_z, MGM_x_new, cores) ;
		} else {
			copy.MGMlist(MGM_x_old, MGM_x_new, cores) ;
		}
		
		t_new	<- (1 + sqrt(1 + 4*t_old^2)) / 2 ;
		
		if (F_z <= F_x_old) {
			add_MGM(MGM_y_new, MGM_x_new, MGM_x_old, 1+(t_old-1)/t_new, -((t_old - 1)/t_new), cores) ;
		} else {
			add_MGM(MGM_y_new, MGM_x_new, MGM_z, 1-(t_old/t_new), t_old/t_new, cores) ;
		}
		
		sq_diff_x	<- rmsd_MGM(MGM_x_new, MGM_x_old, cores) ;
		ll_x_new	<- if (F_z <= F_x_old) ll_z else ll_x_old ;
		pen_x_new	<- if (F_z <= F_x_old) gg_z else pen_x_old ;
		target_new	<- ll_x_new + pen_x_new ;

		if (converge_by_edge) {
			if (compare_MGMlist(MGM_x_new, MGM_x_old, tol_polish, cores)) {
#		if (target_old < Inf & abs(target_new - target_old) / abs(target_old) < 1e-8) {
				edge_diff	<- edge_diff + 1 ;
			} else {
				edge_diff	<- 0 ;
			}
		} else if (sqrt(sq_diff_x$sumSq / sq_diff_x$num) < tol_mgm) {
			edge_diff	<- edge_diff + 1 ;
		} else {
			edge_diff	<- 0 ;
		}
		

		# print(edge_diff) ;
			if (edge_diff >= tol_edge) {
				if (verbose) print("Iteration converged! Returning the results...") ;
				if (verbose) print(paste0("  Final ll: ", ll_x_new / n_tot, ", final penalty: ", pen_x_new / n_tot)) ;
				if (verbose) print(paste0("  Total: ", (ll_x_new + pen_x_new) / n_tot)) ;
				
				break ;
			} else {
				if (verbose) print(paste0("Convergence not achieved, current difference: ", abs(target_new - target_old))) ;
				if (verbose) print(paste0("Current ll: ", ll_x_new / n_tot, ", current penalty: ", pen_x_new / n_tot)) ;
				if (verbose) print(paste0("Total: ", (ll_x_new + pen_x_new) / n_tot)) ;

				copy.MGMlist(MGM_x_new, MGM_x_old, cores) ;
				copy.MGMlist(MGM_y_new, MGM_y_old, cores) ;
				t_old		<- t_new ;
				ll_x_old	<- ll_x_new ;
				pen_x_old	<- pen_x_new ;
				target_old	<- target_new ;
			}

	}
	
	return(MGM_x_new) ;
}

