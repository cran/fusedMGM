require(fastDummies) ;
require(Rcpp) ;
require(parallel) ;
require(bigmemory) ;
require(bigalgebra) ;
require(biganalytics) ;

# Discard unnecessary edges under the given threshold
# Input
## tol_polish: edges under this threshold will be discarded
polish_MGM	<- function(MGM_input, tol_polish=1e-12, cores) {
	X		<- MGM_input$Continuous ;
	Y		<- MGM_input$Discrete ;
	Y_dummy	<- MGM_input$Dummy ;

	p_x	<- ncol(X) ;	p_y	<- ncol(Y) ;
	L_sum		<- ncol(Y_dummy) ;
	L_cumsum	<- c(0, cumsum(MGM_input$Levels)) ;

	# Alpha
	if (p_x > 0) {
		mclapply(seq(p_x), function(s) {
			if (abs(MGM_input$alpha[s,1]) <= tol_polish) { 
				MGM_input$alpha[s,1]		<-  0 
			} ;
			return(NULL) ;
		}, mc.cores=cores) ;
	}
	
	# Beta
	if (p_x > 0) {
		beta_ind	<- expand.grid(seq(p_x), seq(p_x)) ;
		mclapply(seq(nrow(beta_ind)), function(i) {
			s1	<- beta_ind[i,1] ;
			s2	<- beta_ind[i,2] ;
			
			if (abs(MGM_input$beta[s1,s2]) <= tol_polish) {
				MGM_input$beta[s1,s2]	<-  0 ;
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
			
			if (sqrt(mean(rho_tmp^2)) <= tol_polish) {
				MGM_input$rho[ind_tmp,s]	<- 0 ;
			}
			
			return(NULL) ;
		}, mc.cores=cores) ;
			
		rm(rho_ind) ;	gc() ;
	}
	
	# Phi
	if (p_y > 0) {
		mclapply(seq(p_y), function(s) {
			ind_tmp	<- (L_cumsum[s]+1):(L_cumsum[s+1]) ;
			phi_tmp	<- diag(MGM_input$phi[ind_tmp, ind_tmp]);
			
			if (sqrt(mean(phi_tmp^2)) <= tol_polish) {
				for (l in seq(length(ind_tmp))) {
					MGM_input$phi[L_cumsum[s]+l, L_cumsum[s]+l]	<- 0 ;
				}
			}
		}, mc.cores=cores) ;
	
		if (p_y > 1) {
			phi_ind		<- combn(p_y, 2) ;
			
			mclapply(seq(ncol(phi_ind)), function(j) {
				s1	<- phi_ind[1,j] ;
				s2	<- phi_ind[2,j] ;

				ind_tmp_1	<- (L_cumsum[s1]+1):(L_cumsum[s1+1]) ;
				ind_tmp_2	<- (L_cumsum[s2]+1):(L_cumsum[s2+1]) ;
					
				phi_tmp		<- MGM_input$phi[ind_tmp_1, ind_tmp_2] ;
				if (sqrt(mean(phi_tmp^2)) <= tol_polish) {
					MGM_input$phi[ind_tmp_1, ind_tmp_2]	<- MGM_input$phi[ind_tmp_2, ind_tmp_1]	<- 0 ;
				}
				
				return(NULL) ;
			}, mc.cores=cores) ;
			
			rm(phi_ind) ;	gc() ;
		}
	}
	
	return(NULL) ;
}

# Apply the polishing to each of the MGMs in the list
polish_MGM_list	<- function(learn_res, tol_polish=1e-12, cores) {
	G	<- length(learn_res) ;
	
	for (g in seq(G)) {
		polish_MGM(learn_res[[g]], tol_polish, cores) ;
	}
	
	return(learn_res) ;
}

# Convert an MGM to a skeleton network
# Applied to inter-node elements only
MGM_to_skeleton		<- function(MGM_input, inst_temp, inst_ind, tol_polish=1e-12, cores) {
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
				inst_temp$beta_sum[[inst_ind]][s1,s2]	<-  inst_temp$beta_sum[[inst_ind]][s1,s2] + 1 ;
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
				inst_temp$rho_sum[[inst_ind]][r,s]		<- inst_temp$rho_sum[[inst_ind]][r,s] + 1 ;
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
				inst_temp$phi_sum[[inst_ind]][s1,s2]	<- inst_temp$phi_sum[[inst_ind]][s1,s2] + 1 ;
				inst_temp$phi_sum[[inst_ind]][s2,s1]	<- inst_temp$phi_sum[[inst_ind]][s2,s1] + 1 ;
			}
			
			return(NULL) ;
		}, mc.cores=cores) ;
		
		rm(phi_ind) ;	gc() ;
	}
	
	return(NULL) ;
}

# Indicate edges with differetial interactions between 2 MGMs
MGM_to_skel_diff	<- function(MGM_1, MGM_2, inst_temp, inst_ind, tol_polish=1e-12, cores) {
	X		<- MGM_1$Continuous ;
	Y		<- MGM_1$Discrete ;
	Y_dummy	<- MGM_1$Dummy ;

	p_x	<- ncol(X) ;	p_y	<- ncol(Y) ;
	L_sum		<- ncol(Y_dummy) ;
	L_cumsum	<- c(0, cumsum(MGM_1$Levels)) ;
	
	# Beta
	if (p_x > 1) {
		beta_ind	<- expand.grid(seq(p_x), seq(p_x)) ;
		
		mclapply(seq(nrow(beta_ind)), function(i) {
			s1	<- beta_ind[i,1] ;
			s2	<- beta_ind[i,2] ;
			
			if (abs(MGM_1$beta[s1,s2] - MGM_2$beta[s1,s2]) > tol_polish) {
				inst_temp$betadiff_sum[[inst_ind]][s1,s2]	<-  inst_temp$betadiff_sum[[inst_ind]][s1,s2] + 1 ;
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
			rho_tmp_1	<- MGM_1$rho[ind_tmp, s] ;
			rho_tmp_2	<- MGM_2$rho[ind_tmp, s] ;
			
			if (sqrt(mean((rho_tmp_1 - rho_tmp_2)^2)) > tol_polish) {
				inst_temp$rhodiff_sum[[inst_ind]][r,s]	<- inst_temp$rhodiff_sum[[inst_ind]][r,s] + 1 ;
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
				
			phi_tmp_1	<- MGM_1$phi[ind_tmp_1, ind_tmp_2] ;
			phi_tmp_2	<- MGM_2$phi[ind_tmp_1, ind_tmp_2] ;
			if (sqrt(mean((phi_tmp_1 - phi_tmp_2)^2)) > tol_polish) {
				inst_temp$phidiff_sum[[inst_ind]][s1,s2]	<- inst_temp$phidiff_sum[[inst_ind]][s1,s2] + 1 ;
				inst_temp$phidiff_sum[[inst_ind]][s2,s1]	<- inst_temp$phidiff_sum[[inst_ind]][s2,s1] + 1 ;
			}
			
			return(NULL) ;
		}, mc.cores=cores) ;
		
		rm(phi_ind) ;	gc() ;
	}
}

# Make skeleton networks for an MGM list
# Intra: skeleton for edges in each of the MGM
# Inter: skeleton for the differences between the MGMs
add_skeleton	<- function(learn_res, inst_temp, tol_polish=1e-12, cores) {
	G	<- length(learn_res) ;
	
	for (g in seq(G)) {
		MGM_to_skeleton(learn_res[[g]], inst_temp, g, tol_polish, cores) ;
	}
	
	if (G >= 2) {
		for (g2 in seq(2,G)) {
			for (g1 in seq(1,g2-1)) {
				MGM_to_skel_diff(learn_res[[g1]], learn_res[[g2]], inst_temp, 
									sum(0,g2-2) + g1, tol_polish, cores) ;
			}
		}
	}
	
	return(NULL) ;	
}

compare_MGMlist	<- function(MGM_list_1, MGM_list_2, tol_polish=1e-12, cores) {
	G	<- length(MGM_list_1) ;

	X		<- MGM_list_1[[1]]$Continuous ;
	Y		<- MGM_list_1[[1]]$Discrete ;
	Y_dummy	<- MGM_list_1[[1]]$Dummy ;

	p_x	<- ncol(X) ;	p_y	<- ncol(Y) ;
	L_sum		<- ncol(Y_dummy) ;
	L_cumsum	<- c(0, cumsum(MGM_list_1[[1]]$Levels)) ;

	diff_ind 	<- TRUE ;

	for (g in seq(G)) {
		# Beta
		if (p_x > 1) {
			beta_ind	<- combn(p_x, 2) ;
			
			beta_compare	<- unlist(mclapply(seq(ncol(beta_ind)), function(j) {
				s1	<- beta_ind[1,j] ;
				s2	<- beta_ind[2,j] ;
				
				if (abs(MGM_list_1[[g]]$beta[s1,s2]) <= tol_polish) {
					sign_1	<- 0 ;
				} else {
					sign_1	<- sign(MGM_list_1[[g]]$beta[s1,s2]) ;
				}
				
				if (abs(MGM_list_2[[g]]$beta[s1,s2]) <= tol_polish) {
					sign_2	<- 0 ;
				} else {
					sign_2	<- sign(MGM_list_2[[g]]$beta[s1,s2]) ;
				}
				
				return(sign_1 != sign_2) ;
			}, mc.cores=cores)) ;
			
			rm(beta_ind) ;	gc() ;
			
#			print(paste0("beta edge difference: ", sum(beta_compare))) ;

			if (sum(beta_compare) > 0) {
				diff_ind	<- FALSE ;
			}
		}
		
		# Rho
		if (p_x > 0 & p_y > 0) {
			rho_ind		<- expand.grid(seq(p_y), seq(p_x)) ;
			
			rho_compare	<- unlist(mclapply(seq(nrow(rho_ind)), function(i) {
				r	<- rho_ind[i,1] ;
				s	<- rho_ind[i,2] ;
				
				ind_tmp		<- (L_cumsum[r]+1):(L_cumsum[r+1]) ;
				rho_tmp_1	<- MGM_list_1[[g]]$rho[ind_tmp, s] ;
				rho_tmp_2	<- MGM_list_2[[g]]$rho[ind_tmp, s] ;
				
				if (sqrt(mean(rho_tmp_1^2)) <= tol_polish) {
					sign_1	<- 0 ;
				} else {
					sign_1	<- 1 ;
				}
				
				if (sqrt(mean(rho_tmp_2^2)) <= tol_polish) {
					sign_2	<- 0 ;
				} else {
					sign_2	<- 1 ;
				}
				
				return(sign_1 != sign_2) ;
			}, mc.cores=cores)) ;
				
			rm(rho_ind) ;	gc() ;
			
#			print(paste0("rho edge difference: ", sum(rho_compare))) ;

			if (sum(rho_compare) > 0) {
				diff_ind	<- FALSE ;
			}
		}
		
		# Phi
		if (p_y > 1) {
			phi_ind		<- combn(p_y, 2) ;
			
			phi_compare	<- unlist(mclapply(seq(ncol(phi_ind)), function(j) {
				s1	<- phi_ind[1,j] ;
				s2	<- phi_ind[2,j] ;

				ind_tmp_1	<- (L_cumsum[s1]+1):(L_cumsum[s1+1]) ;
				ind_tmp_2	<- (L_cumsum[s2]+1):(L_cumsum[s2+1]) ;
					
				phi_tmp_1	<- MGM_list_1[[g]]$phi[ind_tmp_1, ind_tmp_2] ;
				phi_tmp_2	<- MGM_list_2[[g]]$phi[ind_tmp_1, ind_tmp_2] ;
				
				if (sqrt(mean(phi_tmp_1^2)) <= tol_polish) {
					sign_1	<- 0 ;
				} else {
					sign_1	<- 1 ;
				}
				
				if (sqrt(mean(phi_tmp_2^2)) <= tol_polish) {
					sign_2	<- 0 ;
				} else {
					sign_2	<- 1 ;
				}
				
				return(sign_1 != sign_2) ;
			}, mc.cores=cores)) ;
			
			rm(phi_ind) ;	gc() ;

#			print(paste0("phi edge difference: ", sum(phi_compare))) ;
			
			if (sum(phi_compare) > 0) {
				diff_ind	<- FALSE ;
			}
		}
	}
	
	if (G > 1) {
		for (g2 in seq(2,G)) {
			for (g1 in seq(1,g2-1)) {
				# Beta
				if (p_x > 1) {
					beta_ind	<- combn(p_x, 2) ;

					beta_compare	<- unlist(mclapply(seq(ncol(beta_ind)), function(j) {
						s1	<- beta_ind[1,j] ;
						s2	<- beta_ind[2,j] ;

						if (abs(MGM_list_1[[g1]]$beta[s1,s2] - MGM_list_1[[g2]]$beta[s1,s2]) <= tol_polish) {
							sign_1	<- 0 ;
						} else {
							sign_1	<- sign(MGM_list_1[[g1]]$beta[s1,s2] - MGM_list_1[[g2]]$beta[s1,s2]) ;
						}

						if (abs(MGM_list_2[[g1]]$beta[s1,s2] - MGM_list_2[[g2]]$beta[s1,s2]) <= tol_polish) {
							sign_2	<- 0 ;
						} else {
							sign_2	<- sign(MGM_list_2[[g1]]$beta[s1,s2] - MGM_list_2[[g2]]$beta[s1,s2]) ;
						}

						return(sign_1 != sign_2) ;
					}, mc.cores=cores)) ;

					rm(beta_ind) ;	gc() ;

	#				print(paste0("beta diff difference: ", sum(beta_compare))) ;

					if (sum(beta_compare) > 0) {
						diff_ind	<- FALSE ;
					}
				}

				# Rho
				if (p_x > 0 & p_y > 0) {
					rho_ind		<- expand.grid(seq(p_y), seq(p_x)) ;

					rho_compare	<- unlist(mclapply(seq(nrow(rho_ind)), function(i) {
						r	<- rho_ind[i,1] ;
						s	<- rho_ind[i,2] ;

						ind_tmp		<- (L_cumsum[r]+1):(L_cumsum[r+1]) ;
						rho_tmp_11	<- MGM_list_1[[g1]]$rho[ind_tmp, s] ;
						rho_tmp_12	<- MGM_list_1[[g2]]$rho[ind_tmp, s] ;
						rho_tmp_21	<- MGM_list_2[[g1]]$rho[ind_tmp, s] ;
						rho_tmp_22	<- MGM_list_2[[g2]]$rho[ind_tmp, s] ;

						if (sqrt(mean((rho_tmp_11 - rho_tmp_12)^2)) <= tol_polish) {
							sign_1	<- 0 ;
						} else {
							sign_1	<- 1 ;
						}

						if (sqrt(mean((rho_tmp_21 - rho_tmp_22)^2)) <= tol_polish) {
							sign_2	<- 0 ;
						} else {
							sign_2	<- 1 ;
						}

						return(sign_1 != sign_2) ;
					}, mc.cores=cores)) ;

					rm(rho_ind) ;	gc() ;

	# 				print(paste0("rho diff difference: ", sum(rho_compare))) ;

					if (sum(rho_compare) > 0) {
						diff_ind	<- FALSE ;
					}
				}

				# Phi
				if (p_y > 1) {
					phi_ind		<- combn(p_y, 2) ;

					phi_compare	<- unlist(mclapply(seq(ncol(phi_ind)), function(j) {
						s1	<- phi_ind[1,j] ;
						s2	<- phi_ind[2,j] ;

						ind_tmp_1	<- (L_cumsum[s1]+1):(L_cumsum[s1+1]) ;
						ind_tmp_2	<- (L_cumsum[s2]+1):(L_cumsum[s2+1]) ;

						phi_tmp_11	<- MGM_list_1[[g1]]$phi[ind_tmp_1, ind_tmp_2] ;
						phi_tmp_12	<- MGM_list_1[[g2]]$phi[ind_tmp_1, ind_tmp_2] ;
						phi_tmp_21	<- MGM_list_2[[g1]]$phi[ind_tmp_1, ind_tmp_2] ;
						phi_tmp_22	<- MGM_list_2[[g2]]$phi[ind_tmp_1, ind_tmp_2] ;

						if (sqrt(mean((phi_tmp_11 - phi_tmp_12)^2)) <= tol_polish) {
							sign_1	<- 0 ;
						} else {
							sign_1	<- 1 ;
						}

						if (sqrt(mean((phi_tmp_21 - phi_tmp_22)^2)) <= tol_polish) {
							sign_2	<- 0 ;
						} else {
							sign_2	<- 1 ;
						}

						return(sign_1 != sign_2) ;
					}, mc.cores=cores)) ;

					rm(phi_ind) ;	gc() ;

	#				print(paste0("phi diff difference: ", sum(phi_compare))) ;

					if (sum(phi_compare) > 0) {
						diff_ind	<- FALSE ;
					}
				}
			}
		}
	}
	
	return(diff_ind) ;
}

# Deep copying MGM list
deepcopy.MGMlist	<- function(MGM_list) {
	G	<- length(MGM_list) ;
	MGM_copy	<- rep(list(NULL), G) ;
	
	for (g in seq(G)) {
		MGM_copy[[g]]	<- deepcopy.MGM(MGM_list[[g]]) ;
	}
	
	return(MGM_copy) ;
}

# Deep copying MGM
deepcopy.MGM	<- function(MGM_in) {
	MGM_copy	<- MGM_in ;
	
	MGM_copy$alpha		<- deepcopy(MGM_in$alpha) ;
	MGM_copy$beta		<- deepcopy(MGM_in$beta) ;
	MGM_copy$rho		<- deepcopy(MGM_in$rho) ;
	MGM_copy$phi		<- deepcopy(MGM_in$phi) ;
	
	MGM_copy$Cont_loss	<- deepcopy(MGM_in$Cont_loss) ;
	MGM_copy$Disc_pred	<- deepcopy(MGM_in$Disc_pred) ;
	MGM_copy$Disc_prob	<- deepcopy(MGM_in$Disc_prob) ;
	
	MGM_copy$d_alpha	<- deepcopy(MGM_in$d_alpha) ;
	MGM_copy$d_beta		<- deepcopy(MGM_in$d_beta) ;
	MGM_copy$d_rho		<- deepcopy(MGM_in$d_rho) ;
	MGM_copy$d_phi		<- deepcopy(MGM_in$d_phi) ;	

	return(MGM_copy) ;
}

copy.MGMlist	<- function(MGM_list_send, MGM_list_get, cores) {
	G	<- length(MGM_list_send) ;
	
	for (g in seq(G)) {
		copy.MGM(MGM_list_send[[g]], MGM_list_get[[g]], cores) ;
	}
	
	return(NULL) ;
}

# Copy values in MGM to another
copy.MGM	<- function(MGM_send, MGM_get, cores) {

	X		<- MGM_send$Continuous ;
	Y		<- MGM_send$Discrete ;
	Y_dummy	<- MGM_send$Dummy ;

	n_x	<- nrow(X) ;	p_x	<- ncol(X) ;
	n_y	<- nrow(Y) ;	p_y	<- ncol(Y) ;
	L_sum		<- ncol(Y_dummy) ;
	L_cumsum	<- c(0, cumsum(MGM_send$Levels)) ;

	if (p_x > 0) {
		X_ind	<- expand.grid(seq(n_x), seq(p_x)) ;
		mclapply(seq(nrow(X_ind)), function(i) {
			n	<- X_ind[i,1] ;
			s	<- X_ind[i,2] ;
			MGM_get$Continuous[n,s]	<- MGM_send$Continuous[n,s] ;
			MGM_get$Cont_loss[n,s]	<- MGM_send$Cont_loss[n,s] ;
			return(NULL) ;
		}, mc.cores=cores) ;
		
		rm(X_ind) ;	gc() ;
	}
	
	if (p_y > 0) {
		Y_ind	<- expand.grid(seq(n_y), seq(p_y)) ;
		mclapply(seq(nrow(Y_ind)), function(i) {
			n	<- Y_ind[i,1] ;
			r	<- Y_ind[i,2] ;
			MGM_get$Discrete[n,r]	<- MGM_send$Discrete[n,r] ;
			return(NULL) ;
		}, mc.cores=cores) ;
		
		rm(Y_ind) ;	gc() ;
	}
	
	if (p_y > 0) {
		D_ind	<- expand.grid(seq(n_y), seq(L_sum)) ;
		mclapply(seq(nrow(D_ind)), function(i) {
			n	<- D_ind[i,1] ;
			l	<- D_ind[i,2] ;
			MGM_get$Dummy[n,l]	<- MGM_send$Dummy[n,l] ;
			MGM_get$Disc_pred[n,l]	<- MGM_send$Disc_pred[n,l] ;
			MGM_get$Disc_prob[n,l]	<- MGM_send$Disc_prob[n,l] ;
			return(NULL) ;
		}, mc.cores=cores) ;
		
		rm(D_ind) ;	gc() ;
	}
	
	if (p_x > 0) {
		mclapply(seq(p_x), function(s) {
			MGM_get$alpha[s,1]		<- MGM_send$alpha[s,1] ;
			MGM_get$d_alpha[s,1]	<- MGM_send$d_alpha[s,1] ;
			return(NULL) ;
		}, mc.cores=cores) ;
		
		beta_ind	<- expand.grid(seq(p_x), seq(p_x)) ;
		mclapply(seq(nrow(beta_ind)), function(j) {
			s1	<- beta_ind[j,1] ;
			s2	<- beta_ind[j,2] ;
			
			MGM_get$beta[s1,s2]		<- MGM_send$beta[s1,s2] ;
			MGM_get$d_beta[s1,s2]	<- MGM_send$d_beta[s1,s2] ;
			return(NULL) ;
		}, mc.cores=cores) ;

		rm(beta_ind) ;	gc() ;
	}
	
	if (p_x > 0 & p_y > 0) {
		rho_ind		<- expand.grid(seq(L_sum), seq(p_x)) ;
		mclapply(seq(nrow(rho_ind)), function(i) {
			l	<- rho_ind[i,1] ;
			s	<- rho_ind[i,2] ;
			
			MGM_get$rho[l,s]	<- MGM_send$rho[l,s] ;
			MGM_get$d_rho[l,s]	<- MGM_send$d_rho[l,s] ;
			return(NULL) ;
		}, mc.cores=cores) ;
		
		rm(rho_ind) ;	gc() ;
	}
	
	if (p_y > 0) {
		phi_ind		<- expand.grid(seq(L_sum), seq(L_sum)) ;
		mclapply(seq(nrow(phi_ind)), function(j) {
			r1	<- phi_ind[j,1] ;
			r2	<- phi_ind[j,2] ;
		
			MGM_get$phi[r1,r2]		<- MGM_send$phi[r1,r2] ;
			MGM_get$d_phi[r1,r2]	<- MGM_send$d_phi[r1,r2] ;
			return(NULL) ;
		}, mc.cores=cores) ;
		
		rm(phi_ind) ;	gc() ;
	}
	
	return(NULL) ;
}

