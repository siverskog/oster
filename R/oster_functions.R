################################################################################
##### oster                                      ###############################
##### "Bias-adjusted treatments effects under    ###############################
#####  proportional selection on observables     ###############################
#####  and unobservables: an R adaption of       ###############################
#####  Stata's psacalc"                          ###############################
#####                                            ###############################
##### Version: 2020-06-03                        ###############################
##### Jonathan Siverskog                         ###############################
##### Linkoping University                       ###############################
################################################################################
################################################################################
#
# ...
#
################################################################################
##### MAIN FUNCTION ############################################################
################################################################################

#' Oster bias-adjusted beta and bounds for maximum r-squared and delta
#'
#' This function calculates the bias-adjusted beta proposed by Oster (2019).
#' Adapted for R from Stata's psacalc (Oster, 2013)
#'
#' @param o.fit object of class \code{lm} with fit for short regression
#' @param tilde.fit object of class \code{lm} with fit for intermediate regression
#' @param varname character specifying the name of the treatment variable, e.g. "x" or "log(x)"
#' @param rm numeric specifying maximum r-squared. If \code{rm<=1} then \code{rmax=rm},
#'        if \code{rm>1} then \code{rmax = r_tilde*rm}, if \code{NULL} (default) then \code{rmax = 2*r_tilde - r_o}
#' @param delta numeric specifying the proportional selection on observables and unobservables. Default is 1.
#' @param beta numeric used to solve for bounds on maximum r-squared and delta. Default is 0.
#' @param within logical specifying whether overall r-sqaured (default) or within r-squared is used
#' @param fe character vector of names of variables used for calculating within r-squared. Default is to use \code{mcontrols}.
#' @param wvar logical specifying whether variances should be weighted. Ignored if
#'        \code{o.fit} and \code{tilde.fit} do not include weights. Default is \code{FALSE};
#'        this action is not performed by Stata's \code{psacalc}.
#' @return list containing bias-adjusted beta, bounds on delta and maximum r-squared
#'         for beta (and delta) supplied, and other inputs and output
#' @details If covariates are included in both \code{o.fit} and \code{tilde.fit},
#'          the treatment variable is residualised with respect to these;
#'          compare to \code{mcontrols} in Stata's psacalc. There are 1 to 3
#'          solutions for beta. All solution are reported in \code{allroots}.
#'          \code{dist} is the (squared) distance of each solution to \code{beta_tilde}
#'          and \code{mark} is 0 if the solution has the same direction of bias
#'          compared to \code{beta_tilde} as \code{beta_tilde} versus \code{beta_o}.
#'          \code{input} contains user supplied inputs, and inputs extracted from
#'          the fits supplied.
#' @references Oster, E. (2019). Unobservable selection and coefficient
#'             stability: Theory and evidence. Journal of Business & Economic
#'             Statistics, 37(2), 187-204
#' @references Oster, E. (2013). PSACALC: Stata module to calculate treatment
#'             effects and relative degree of selection under proportional
#'             selection of observables and unobservables. Statistical Software
#'             Components S457677, Boston College Department of Economics,
#'             revised 18 Dec 2016.
#' @examples
#' 
#' # EXAMPLE REPRODUCIBLE IN STATA'S PSACALC
#' 
#' data(nunn)
#'  
#' fit_a <- lm(trust_relatives ~ ln_export_area + age + male + as.factor(education)
#'             + as.factor(isocode), data = nunn)
#'  
#' nunn$infit_a <- is.element(rownames(nunn), names(fit_a$residuals))
#' 
#' fit_b <- lm(trust_relatives ~ ln_export_area + as.factor(isocode),
#' data = nunn, subset = infit_a)
#'
#' oster(fit_b, fit_a, varname = "ln_export_area", rm = 0.2)
#' 
#' # REPRODUCE ROW 1 OF TABLE 5 IN OSTER (2019)
#' 
#' f01c <- trust_relatives ~ ln_export_area + age + age2 + male + urban_dum +
#'  as.factor(education) + as.factor(occupation) + as.factor(religion) +
#'  as.factor(living_conditions) + district_ethnic_frac +
#'  frac_ethnicity_in_district + as.factor(isocode) + malaria_ecology +
#'  total_missions_area + explorer_contact + railway_contact + cities_1400_dum +
#'  as.factor(v30) + v33_alt + ln_init_pop_density
#'
#' f01u <- trust_relatives ~ ln_export_area + as.factor(isocode)
#' 
#' fit01c <- lm(f01c, data = nunn)
#' nunn$infit01c <- is.element(rownames(nunn), names(fit01c$residuals))
#' fit01u <- lm(f01u, data = nunn, subset = infit01c)
#' 
#' z <- oster(fit01u, fit01c, "ln_export_area")
#' b13 <- oster(fit01u, fit01c, "ln_export_area", rm = 1.3)$beta
#' c(z$input$beta_o, z$input$beta_tilde, z$beta, b13, z$rmax)
#' 
#' @export
oster <- function(o.fit, tilde.fit, varname, rm = NULL, delta = 1, beta = 0, within = FALSE, fe = NULL, wvar = FALSE) {

  check_cases <- all(is.element(case.names(o.fit), case.names(tilde.fit))) & all(is.element(case.names(tilde.fit), case.names(o.fit)))
  if(!check_cases) {warning("short and intermediate fit do not include the same cases")}

  ### EXTRACT BETAS ###

  beta_o <- o.fit$coef[varname]
  beta_tilde <- tilde.fit$coef[varname]

  ### RESIDUALISE X WRT CONTROLS IN INTERMEDIATE (TILDE) FIT ###

  tau.fit <- update(o.fit, update(formula(tilde.fit), as.formula(paste(varname, "~ . -", varname))))
  
  if(wvar & length(tau.fit$weights)>0) {
    t_x <- weighted.var(tau.fit$residuals, tau.fit$weights)
  } else {
    t_x <- var(tau.fit$residuals)
  }

  ### CHECK FOR COVARIATES IN BOTH FITS AND RESIDUALISE X WRT TO THESE ###

  m <- intersect(colnames(o.fit$model), colnames(tilde.fit$model))
  m <- m[!is.element(m, c(varname, "(weights)", colnames(o.fit$model)[1]))]

  if(length(m)>0) {

    x.fit <- update(o.fit, as.formula(paste(varname, "~", paste(m, collapse = " + "))))
    
    if(wvar & length(x.fit$weights)>0) {
      sigma_xx <- weighted.var(x.fit$residuals, x.fit$weights)
    } else {
      sigma_xx <- var(x.fit$residuals)
    }
    
  } else {
    
    if(wvar & length(o.fit$weights)>0) {
      sigma_xx <- weighted.var(o.fit$model[,varname], o.fit$model[,"(weights)"])
    } else {
      sigma_xx <- var(o.fit$model[,varname])
    }

  }
  
  ### VARIANCE OF OUTCOME ###
  
  if(wvar & length(o.fit$weights)>0) {
    sigma_yy <- weighted.var(o.fit$model[,1], o.fit$model[,"(weights)"])
  } else {
    sigma_yy <- var(o.fit$model[,1])
  }

  ### R2 OR WITHIN R2 ###
  
  if(within & length(fe)>0) {
    
    r_o <- within.rsq(o.fit, fe)
    r_tilde <- within.rsq(tilde.fit, fe)
    
  } else if(within & length(m)>0) {
    
    r_o <- within.rsq(o.fit, all.vars(formula(x.fit))[-1])
    r_tilde <- within.rsq(tilde.fit, all.vars(formula(x.fit))[-1])
    
  } else if(within) {
    
    stop("If within=TRUE, then fe must be specified or o.fit and tilde.fit must have covariates in common")
    
  } else {
    
    r_o <- summary(o.fit)$r.squared
    r_tilde <- summary(tilde.fit)$r.squared
    
  }

  ### SET R2-MAX ###

  if(is.null(rm)) {
    rmax <- 2*r_tilde - r_o
  } else if(rm>1) {
    rmax <- r_tilde*rm
  } else {
    rmax <- rm
  }

  ### DEFINE A FEW TERMS ###

  bo_m_bt <- beta_o - beta_tilde
  rt_m_ro_t_syy <- (r_tilde - r_o)*sigma_yy
  rm_m_rt_t_syy <- (rmax - r_tilde)*sigma_yy

  ### CALCULATION ###

  if(delta==1) {
    bout <- d1quadsol(rm_m_rt_t_syy, rt_m_ro_t_syy, bo_m_bt, sigma_xx, t_x, beta_o, beta_tilde)
  } else if (is.numeric(delta)) {
    bout <- dnot1cubsol(bo_m_bt, sigma_xx, delta, t_x, rm_m_rt_t_syy, rt_m_ro_t_syy, beta_tilde, beta_o)
  } else {
    stop("delta must be numeric to calculate bias-adjusted beta")
  }

  if(is.numeric(beta) & is.numeric(delta)) {
    dout <- bound(bo_m_bt, rt_m_ro_t_syy, rm_m_rt_t_syy, beta_tilde, beta, t_x, sigma_xx)
    rout <- rbound(bo_m_bt, rt_m_ro_t_syy, beta_tilde, beta, t_x, sigma_xx, sigma_yy, delta, r_tilde)
  } else {
    stop("beta and delta must be numeric to calculate bounds")
  }

  boutx <- bout$beta[1]
  names(bout)[1] <- "allroots"
  input <- c(beta_o, beta_tilde, r_o, r_tilde, rmax, delta, beta)
  names(input) <- c("beta_o", "beta_tilde", "r_o", "r_tilde", "rmax", "delta", "beta")
  input <- as.list(input)

  ret <- list(beta = boutx, delta = dout, rmax = rout, output = bout, input = input)

  return(ret)

}

############################################################################
##### SOLUTION DELTA BOUND #################################################
############################################################################

#' Delta bound
#'
#' Solves for delta. Used in \code{oster()}
#'
#' @param ... arguments passed via \code{oster()}
#' @return numeric for delta
#' @examples
#' ...
bound <- function(bo_m_bt, rt_m_ro_t_syy, rm_m_rt_t_syy, beta_tilde, beta, t_x, sigma_xx) {

  bt_m_b = beta_tilde - beta

  num1 <- (bt_m_b)*rt_m_ro_t_syy*t_x
  num2 <- (bt_m_b)*sigma_xx*t_x*(bo_m_bt)^2
  num3 <- 2*(bt_m_b^2)*(t_x*bo_m_bt*sigma_xx)
  num4 <- (bt_m_b^3)*(t_x*sigma_xx-(t_x^2))

  den1 <- rm_m_rt_t_syy*bo_m_bt*sigma_xx
  den2 <- bt_m_b*rm_m_rt_t_syy*(sigma_xx-t_x)
  den3 <- (bt_m_b^2)*(t_x*bo_m_bt*sigma_xx)
  den4 <- (bt_m_b^3)*(t_x*sigma_xx-(t_x^2))

  num <- num1+num2+num3+num4
  den <- den1+den2+den3+den4
  ret <- num/den
  names(ret) <- NULL

  return(ret)

}

############################################################################
##### SOLUTION RSQ BOUND ###################################################
############################################################################

#' R-squared bound
#'
#' Solves for maximum r-squared. Used in \code{oster()}
#'
#' @param ... arguments passed via \code{oster()}
#' @return numeric for maximum r-squared
#' @examples
#' ...
rbound <- function(bo_m_bt, rt_m_ro_t_syy, beta_tilde, beta, t_x, sigma_xx, sigma_yy, delta, r_tilde) {

  bt_m_b = beta_tilde - beta

  num1 <- (bt_m_b)*rt_m_ro_t_syy*t_x
  num2 <- (bt_m_b)*sigma_xx*t_x*(bo_m_bt)^2
  num3 <- 2*(bt_m_b^2)*(t_x*bo_m_bt*sigma_xx)
  num4 <- (bt_m_b^3)*(t_x*sigma_xx-(t_x^2))
  num5 <- (bt_m_b^2)*(t_x*bo_m_bt*sigma_xx)*(-delta)
  num6 <- (bt_m_b^3)*(t_x*sigma_xx-(t_x^2))*(-delta)

  den1 <- delta*sigma_yy*bo_m_bt*sigma_xx
  den2 <- delta*bt_m_b*sigma_yy*(sigma_xx-t_x)

  num <- num1+num2+num3+num4+num5+num6
  den <- den1+den2
  ret <- (num/den) + r_tilde
  names(ret) <- NULL

  return(ret)

}

############################################################################
##### SOLUTION IF DELTA = 1 ################################################
############################################################################

#' Solve for beta with delta eq to 1
#'
#' Solves for beta. Used in \code{oster()}
#'
#' @param ... arguments passed via \code{oster()}
#' @return list containing all solutions, their distances from beta_tilde and
#'         indicator of whether the solution goes in opposite direction of
#'         observed bias
#' @examples
#' ...
d1quadsol <- function(rm_m_rt_t_syy, rt_m_ro_t_syy, bo_m_bt, sigma_xx, t_x, beta_o, beta_tilde) {

  cap_theta <- rm_m_rt_t_syy*(sigma_xx-t_x) - rt_m_ro_t_syy*t_x-sigma_xx*t_x*(bo_m_bt^2)
  d1_1 <- 4*rm_m_rt_t_syy*(bo_m_bt^2)*(sigma_xx^2)*t_x
  d1_2 <- -2*t_x*bo_m_bt*sigma_xx

  sol1 <- (-1*cap_theta - sqrt((cap_theta^2)+d1_1))/(d1_2)
  sol2 <- (-1*cap_theta + sqrt((cap_theta^2)+d1_1))/(d1_2)

  beta1 <- beta_tilde - sol1
  beta2 <- beta_tilde - sol2

  if ( (beta1-beta_tilde)^2 < (beta2-beta_tilde)^2) {
    betax <- beta1
    altsol1 <- beta2
  } else {
    betax <- beta2
    altsol1 <- beta1
  }

  if ( sign(betax-beta_tilde)!=sign(beta_tilde-beta_o) ) {
    solc <- betax
    betax <- altsol1
    altsol1 <- solc
  }

  beta <- c(betax, altsol1)

  if ( sign(betax-beta_tilde)!=sign(beta_tilde-beta_o) ) { markx <- 1 } else { markx <- 0 }
  if ( sign(altsol1-beta_tilde)!=sign(beta_tilde-beta_o) ) { mark1 <- 1 }  else { mark1 <- 0 }

  mark <- c(markx, mark1)
  dist <- c((betax - beta_tilde)^2, (altsol1 - beta_tilde)^2)

  names(beta) = names(dist) = names(mark) = NULL

  return(list(beta = beta, dist = dist, mark = mark))

}

############################################################################
##### SOLUTION IF DELTA =! 1 ###############################################
############################################################################

#' Solve for beta with delta neq to 1
#'
#' Solves for beta. Used in \code{oster()}
#'
#' @param ... arguments passed via \code{oster()}
#' @return list containing all solutions, their distances from beta_tilde and
#'         indicator of whether the solution goes in opposite direction of
#'         observed bias
#' @examples
#' ...
dnot1cubsol <- function(bo_m_bt, sigma_xx, delta, t_x, rm_m_rt_t_syy, rt_m_ro_t_syy,beta_tilde, beta_o) {

  A <- (t_x*bo_m_bt*sigma_xx*(delta-2))/((delta-1)*(t_x*sigma_xx-t_x^2))
  B <- (delta*rm_m_rt_t_syy*(sigma_xx-t_x)-rt_m_ro_t_syy*t_x-sigma_xx*t_x*bo_m_bt^2)/((delta-1)*(t_x*sigma_xx-t_x^2)) # bo_m_bt'^2 ?
  C <- (rm_m_rt_t_syy*delta*bo_m_bt*sigma_xx)/((delta-1)*(t_x*sigma_xx-t_x^2))

  Q <- (A^2-3*B)/9
  R <- (2*A^3-9*A*B+27*C)/54
  D <- R^2-Q^3
  discrim <- R^2-Q^3

  if (discrim<0) {

    theta <-  acos(R/sqrt(Q^3))

    sol1 <- -2*sqrt(Q)*cos(theta/3)-(A/3)
    sol2 <- -2*sqrt(Q)*cos((theta+2*pi)/3)-(A/3)
    sol3 <- -2*sqrt(Q)*cos((theta-2*pi)/3)-(A/3)

    sols <- beta_tilde - c(sol1,sol2,sol3)
    dists <- (sols - beta_tilde)^2
    marks <- as.numeric(sign(sols-beta_tilde)!=sign(beta_tilde-beta_o))

    ### ORDER SOLUTIONS ###

    sols <- sols[order(dists)]
    marks <- marks[order(dists)]
    dists <- dists[order(dists)]

    beta <- sols[order(marks)]
    dist <- dists[order(marks)]
    mark <- marks[order(marks)]

  } else {

    t1=-1*R+sqrt(D)
    t2=-1*R-sqrt(D)

    crt1 <- sign(t1) * abs(t1)^(1/3)
    crt2 <- sign(t2) * abs(t2)^(1/3)

    sol1 <- crt1+crt2-(A/3)
    betax <- beta_tilde-sol1

    dist <- (betax - beta_tilde)^2
    mark <- as.numeric(sign(betax-beta_tilde)!=sign(beta_tilde-beta_o))
    beta <- betax

  }

  names(beta) = names(dist) = names(mark) = NULL
  return(list(beta = beta, dist = dist, mark = mark))

}

############################################################################
##### WEIGHTED VARIANCE ####################################################
############################################################################

#' Compute Weighted Variance
#'
#' what it says
#'
#' @param x a numeric vector of values
#' @param w a numeric vector of weights
#' @return numeric for weighted variance
#' @examples
#' ...
weighted.var <- function(x, w, na.rm = FALSE) {
  
  if(length(w)!=length(n)) {warning("weights and values of unequal length")}
  w <- (w/sum(w))*length(w)
  ret <- sum(w*((x-weighted.mean(x, w))^2))/(length(x)-1)
  return(ret)
  
}

############################################################################
##### WITHIN R-SQUARED #####################################################
############################################################################

#' Compute within R-squared
#'
#' what it says
#'
#' @param lm.fit an object of class \code{lm}.
#' @param fe character vector specifying the variable names for the fixed effects
#' @return numeric for within r-squared.
#' @examples
#' ...
within.rsq <- function(lm.fit, fe) {
  
  frm <- model.frame(lm.fit)
  
  y <- model.extract(frm, component = "response")
  X <- model.matrix(lm.fit)
  w <- model.extract(frm, component = "weights")
  
  f <- formula(lm.fit)
  f <- update(f, as.formula(paste(". ~", paste(paste("as.factor(", fe, ")", sep = ""), collapse = " + "))))
  tmp <- update(lm.fit, f)
  femat <- as.matrix(as.data.frame(model.matrix(tmp)))
  cases <- case.names(lm.fit)
  femat <- femat[is.element(rownames(femat), cases), , drop = FALSE]
  
  if(length(w)==0) {
    
    xfit <- lm.fit(x = femat, y = X)
    yfit <- lm.fit(x = femat, y = y)
    
  } else {
    
    xfit <- lm.wfit(x = femat, y = X, w = w)
    yfit <- lm.wfit(x = femat, y = y, w = w)
    
  }
  
  newX <- as.matrix(as.data.frame(xfit$residuals))
  newX <- cbind(newX[,!is.element(colnames(newX), colnames(femat)), drop = FALSE], femat)
  newy <- yfit$residuals
  
  if(length(w)==0) {
    
    newfit <- lm.fit(x = newX, y = newy)
    sst <- var(newy)
    sse <- var(newfit$fitted.values)
    
  } else {
    
    newfit <- lm.wfit(x = newX, y = newy, w = w)
    sst <- weighted.var(newy, w)
    sse <- weighted.var(newfit$fitted.values, w)
    
  }
  
  return(sse/sst)
  
}
