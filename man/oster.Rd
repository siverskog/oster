% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oster_functions.R
\name{oster}
\alias{oster}
\title{Oster bias-adjusted beta and bounds for maximum r-squared and delta}
\usage{
oster(
  o.fit,
  tilde.fit,
  varname,
  rm = NULL,
  delta = 1,
  beta = 0,
  within = FALSE,
  fe = NULL,
  wvar = FALSE
)
}
\arguments{
\item{o.fit}{object of class \code{lm} with fit for short regression}

\item{tilde.fit}{object of class \code{lm} with fit for intermediate regression}

\item{varname}{character specifying the name of the treatment variable, e.g. "x" or "log(x)"}

\item{rm}{numeric specifying maximum r-squared. If \code{rm<=1} then \code{rmax=rm},
if \code{rm>1} then \code{rmax = r_tilde*rm}, if \code{NULL} (default) then \code{rmax = 2*r_tilde - r_o}}

\item{delta}{numeric specifying the proportional selection on observables and unobservables. Default is 1.}

\item{beta}{numeric used to solve for bounds on maximum r-squared and delta. Default is 0.}

\item{within}{logical specifying whether overall r-sqaured (default) or within r-squared is used}

\item{fe}{character vector of names of variables used for calculating within r-squared. Default is to use \code{mcontrols}.}

\item{wvar}{logical specifying whether variances should be weighted. Ignored if
\code{o.fit} and \code{tilde.fit} do not include weights. Default is \code{FALSE};
this action is not performed by Stata's \code{psacalc}.}
}
\value{
A list with components
\item{beta}{numeric; bias-adjusted beta}
\item{delta}{numeric; bound on delta for beta and rmax supplied}
\item{rmax}{numeric; bound on rmax for beta and delta supplied}
\item{output}{list; allroots contains a vector with all solutions for beta,
dist contains the (squared) distance from beta_tilde for each solution,
mark contains an indicator for each solution, 1 if the solution has
the opposite direction of bias compared to beta_tilde as beta_tilde
versus beta_o}
\item{input}{list; contains user supplied inputs, and inputs extracted from
the fits supplied}
}
\description{
This function calculates the bias-adjusted beta proposed by Oster (2019).
Adapted for R from Stata's psacalc (Oster, 2013)
}
\details{
If covariates are included in both \code{o.fit} and \code{tilde.fit},
the treatment variable is residualised with respect to these;
compare to \code{mcontrols} in Stata's psacalc. There are 1 to 3
solutions for beta. \code{oster()} returns the solution with minimum
\code{dist} where \code{mark} is 0.
}
\examples{

# EXAMPLE REPRODUCIBLE IN STATA'S PSACALC

data(nunn)
 
fit_a <- lm(trust_relatives ~ ln_export_area + age + male + as.factor(education)
            + as.factor(isocode), data = nunn)
 
nunn$infit_a <- is.element(rownames(nunn), names(fit_a$residuals))

fit_b <- lm(trust_relatives ~ ln_export_area + as.factor(isocode),
data = nunn, subset = infit_a)

oster(fit_b, fit_a, varname = "ln_export_area", rm = 0.2)

# REPRODUCE ROW 1 OF TABLE 5 IN OSTER (2019)

f01c <- trust_relatives ~ ln_export_area + age + age2 + male + urban_dum +
 as.factor(education) + as.factor(occupation) + as.factor(religion) +
 as.factor(living_conditions) + district_ethnic_frac +
 frac_ethnicity_in_district + as.factor(isocode) + malaria_ecology +
 total_missions_area + explorer_contact + railway_contact + cities_1400_dum +
 as.factor(v30) + v33_alt + ln_init_pop_density

f01u <- trust_relatives ~ ln_export_area + as.factor(isocode)

fit01c <- lm(f01c, data = nunn)
nunn$infit01c <- is.element(rownames(nunn), names(fit01c$residuals))
fit01u <- lm(f01u, data = nunn, subset = infit01c)

z <- oster(fit01u, fit01c, "ln_export_area")
b13 <- oster(fit01u, fit01c, "ln_export_area", rm = 1.3)$beta
c(z$input$beta_o, z$input$beta_tilde, z$beta, b13, z$rmax)

}
\references{
Oster, E. (2019). Unobservable selection and coefficient
stability: Theory and evidence. Journal of Business & Economic
Statistics, 37(2), 187-204

Oster, E. (2013). PSACALC: Stata module to calculate treatment
effects and relative degree of selection under proportional
selection of observables and unobservables. Statistical Software
Components S457677, Boston College Department of Economics,
revised 18 Dec 2016.
}
