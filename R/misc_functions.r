# Exponential sum
logsumexp	<- function(vec) {
	maxVal	<- max(vec) ;
	return(log(sum(exp(vec-maxVal))) + maxVal) ;
}

# Fixed point approach
## init: initial point
## fun: function for iteration
## target: target function to be converged
## maxit: maximum number of iterations
## tol: tolerance of convergence
fpa	<- function(init, fun, target, maxit=1000000, tol=1e-12) {
	prevval		<- init ;
	newval		<- 0 ;
	conv_flag	<- FALSE ;
	
	for (i in seq(maxit)) {
		newval	<- fun(prevval) ;
		if (((sqrt(mean((prevval - newval)^2)) < tol)) | (abs(target(prevval) - target(newval)) < tol)) {
			conv_flag	<- TRUE ;
			break ;
		} else {
			prevval	<- newval ;
		}
	}
	
	return(list(value=newval, converged=conv_flag)) ;
}

