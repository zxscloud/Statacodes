/*************************************************
** Author: Joseph Terza				**
** Date: December 8, 2016                  	**
** Purpose: Mullahy (1997) Birth Weight model.	**
** Estimation of the model using the 2SRI. 	**
** The outcome variable (birth weight) is 	**
** non-negative and continuous, and the		**
** endogenous variable (cigarette smoking) is 	**
** non-negative and a follow a two-part modeling**
** specification.								**
*************************************************/
/*************************************************
** Preliminary Stuff.            		**
*************************************************/
clear
clear mata
clear matrix
macro drop _all
set more off
capture log close
clear programs

/*************************************************
** Set up default directory.			**
*************************************************/

/*************************************************
** Set up the output file.			**
*************************************************/
sjlog using terza-two-stage1, replace           
/*************************************************
** Read in the data.				**
*************************************************/
use mullahy-birthweight-data-lbs-not-oz

/*************************************************
** Generate the binary smoking variable.	**
*************************************************/
generate ANYCIGS=CIGSPREG>0

/*************************************************
** 2SRI first-stage first-part probit estimates.**
*************************************************/
/*Step c*/
probit ANYCIGS PARITY WHITE MALE EDFATHER EDMOTHER ///
	FAMINCOM CIGTAX88
	
/*************************************************
** Save the 2SRI first-stage first-part probit	**
** predicted values for use in calculating	**
** the first stage residuals.			**
*************************************************/
/*Step d*/
predict CIGPROB

/*************************************************
** Save the first-stage first-part probit 	**
** estimates and estimated covariance matrix.	**
*************************************************/
/*Step e*/
mata: alpha1hat=st_matrix("e(b)")'
mata: Valpha1hat=st_matrix("e(V)")

/*************************************************
** 2SRI first-stage second-part probit NLS 	**
** estimates.					**
*************************************************/
/*Step c*/
glm CIGSPREG PARITY WHITE MALE EDFATHER EDMOTHER ///
	FAMINCOM CIGTAX88 if ANYCIGS==1, ///
	family(gaussian) link(log) vce(robust)
	
/*************************************************
** Save the 2SRI first-stage second-part NLS 	**
** (glm) predicted values for use in calculating**
** the first-stage residuals.			**
*************************************************/
/*Step d*/
predict CIGMEAN

/*************************************************
** Generate the first-stage residuals.		**
*************************************************/
/*Step d*/
generate Xuhat=CIGSPREG-CIGPROB*CIGMEAN

/*************************************************
** Save the first-stage second-part NLS		**
** estimates and estimated covariance matrix    **
*************************************************/
/*Step e*/
mata: alpha2hat=st_matrix("e(b)")'
mata: Valpha2hat=st_matrix("e(V)")

/*************************************************
** Descriptive statistics.			**
*************************************************/
summ

/*************************************************
** 2SRI second-stage NLS estimates.		**
*************************************************/
/*Step f*/
glm BIRTHWTLB CIGSPREG PARITY WHITE MALE Xuhat, ///
	family(gaussian) link(log) vce(robust)
	
/*************************************************
** Save second-stage estimates and covariance 	**
** matrix.  Single out the coefficient estimate	**
** for Xu.					**
*************************************************/
/*Step g*/
mata: betahat=st_matrix("e(b)")'
mata: Vbetahat=st_matrix("e(V)")
mata: Bu=betahat[5]

/*************************************************
** Send the requisite variables to Mata as	**
** vectors.				 	**
*************************************************/
/*Step h*/
putmata BIRTHWTLB CIGSPREG ANYCIGS PARITY WHITE ///
	MALE EDFATHER EDMOTHER FAMINCOM CIGTAX88 Xuhat 
	
/*************************************************
** Use these vectors to concatenate the needed	**
** matrices.					**
*************************************************/
/*Step h*/
mata: X=CIGSPREG, PARITY, WHITE, MALE, ///
	Xuhat, J(rows(PARITY),1,1)
mata: W=PARITY, WHITE, MALE, EDFATHER, EDMOTHER, ///
	FAMINCOM, CIGTAX88, J(rows(PARITY),1,1)

/*************************************************
** Set up the two gradient matrices for the ACSE**
*************************************************/
/*Step i*/
mata: gradbeta=exp(X*betahat):*X
mata: gradalpha1=-Bu:*exp(X*betahat):*normalden(W*alpha1hat):*exp(W*alpha2hat):*W
mata: gradalpha2=-Bu:*exp(X*betahat):*normal(W*alpha1hat):*exp(W*alpha2hat):*W
mata: gradalpha=gradalpha1,gradalpha2

/*************************************************
** Set up the B1 and B2 matrices for the ACSE.	**
*************************************************/
/*Step j*/
mata: B1=gradbeta'*gradbeta
mata: B2=gradbeta'*gradalpha

/*************************************************
** Set up the full estimated asymptotic		**
** covariance matrix for alpha (first-stage 	**
** two-part model covariance matrix estimator as**
** output by Stata).				**
*************************************************/
/*Step k*/
mata: Valphahat=blockdiag(Valpha1hat,Valpha2hat)
	
/*************************************************
** Construct the covariance matrix of the	**
** second-stage Beta estimates.			**
*************************************************/
/*Step k*/
mata: Dhat=invsym(B1)*B2*Valphahat*B2'*invsym(B1)+Vbetahat

/*************************************************
** Extract the vector of asymptotically correct	**
** standard errors for betahat.		**
*************************************************/
/*Step l*/
mata: ACSE=sqrt(diagonal(Dhat))

/*************************************************
** Calculate the corresponding vector of	**
** asymptotically correct t-stats.		**
*************************************************/
/*Step m*/
mata: ACtstats=betahat:/ACSE

/*************************************************
** Compute the corresponding vector of p-values.
*************************************************/
mata: ACpvalues=2:*(1:-normal(abs(ACtstats)))

/*************************************************
** Display results.				**
*************************************************/
mata: header="Variable","Estimate","ACSE","AC t-stat","pvalue"
mata: varnames="CIGSPREG", "PARITY", "WHITE", "MALE", "Xuhat","Constant"
mata: results=betahat,ACSE,ACtstats,ACpvalues
mata: resview=strofreal(results)
mata: "2SRI Results with ACSE"
mata: header \ (varnames',resview)
sjlog close, replace

