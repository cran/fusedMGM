% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MGM.r
\name{FMGM_mc}
\alias{FMGM_mc}
\title{Main function of fused MGM}
\usage{
FMGM_mc(
  data,
  ind_disc,
  group,
  t = 1,
  L = NULL,
  eta = 2,
  lambda_intra,
  lambda_intra_prior = NULL,
  lambda_inter,
  with_prior = FALSE,
  prior_list = NULL,
  converge_by_edge = TRUE,
  tol_edge = 3,
  tol_mgm = 1e-04,
  tol_g = 0.005,
  tol_fpa = 1e-12,
  maxit = 1e+06,
  polish = TRUE,
  tol_polish = 1e-12,
  cores = parallel::detectCores(),
  verbose = FALSE
)
}
\arguments{
\item{data}{Data frame with rows as observations and columns as variables}

\item{ind_disc}{Indices of discrete variables}

\item{group}{Group indices, must be provided with the observation names}

\item{t}{Numeric. Initial value of coefficient that reflect 2 previous iterations in fast proximal gradient method. Default: 1}

\item{L}{Numeric. Initial guess of Lipschitz constant. Default: missing (use backtracking)}

\item{eta}{Numeric. Multipliers for L in backtracking. Default: 2}

\item{lambda_intra}{Vector with 3 numeric variables. Penalization parameters for network edge weights}

\item{lambda_intra_prior}{Vector with 3 numeric variables. Penalization parameters for network edge weights, applied to the edges with prior information}

\item{lambda_inter}{Vector with 3 numeric variables. Penalization parameters for network edge weight differences}

\item{with_prior}{Logical. Is prior information provided? Default: FALSE}

\item{prior_list}{List of prior information. Each element must be a 3-column data frames, with the 1st and the 2nd columns being variable names and the 3rd column being prior confidence (0,1)}

\item{converge_by_edge}{Logical. The convergence should be judged by null differences of network edges after iteration. If FALSE, the rooted mean square difference (RMSD) of edge weights is used. Default: TRUE}

\item{tol_edge}{Integer. Number of consecutive iterations of convergence to stop the iteration. Default: 3}

\item{tol_mgm}{Numeric. Cutoff of network edge RMSD for convergence. Default: 1e-04}

\item{tol_g}{Numeric. Cutoff of iternations in prox-grad map calculation. Default: 5e-03}

\item{tol_fpa}{Numeric. Cutoff for fixed-point approach. Default: 1e-12}

\item{maxit}{Integer. Maximum number of iterations in fixed-point approach. Default: 1000000}

\item{polish}{Logical. Should the edges with the weights below the cutoff should be discarded? Default: TRUE}

\item{tol_polish}{Numeric. Cutoff of polishing the resulting network. Default: 1e-12}

\item{cores}{Integer. Number of cores to use multi-core utilization. Default: maximum number of available cores}

\item{verbose}{Logical. If TRUE, the procedures are reported in real-time manner. Default: FALSE}
}
\value{
The resulting networks, in the form of a list of MGMs
}
\description{
Infers networks from 2-class mixed data
}
\details{
If the value of Lipschitz constant, L, is not provided, the backtracking will be performed
}
\examples{
\donttest{
data(data_all) ;  # Example 500-by-100 simulation data
data(ind_disc) ;

group <- rep(c(1,2), each=250) ;
names(group) <- seq(500) ;

if (Sys.info()['sysname'] == 'Windows') {
  cores=1
} else {
  cores=parallel::detectCores() ;
}

res_FMGM <- FMGM_mc(data_all, ind_disc, group, 
                    lambda_intra=c(0.2,0.15,0.1), lambda_inter=c(0.2,0.15,0.1), 
                    cores=cores, verbose=TRUE)
}
}
