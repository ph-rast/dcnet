<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- knit with rmarkdown::render("README.Rmd", output_format = "md_document") -->
<!-- <\!-- badges: start -\-> -->
<!--   [![R-CMD-check](https://github.com/ph-rast/dcnet/workflows/R-CMD-check/badge.svg)](https://github.com/ph-rast/dcnet/actions) -->
<!-- [![R-CMD-check](https://github.com/ph-rast/dcnet/workflows/R-CMD-check/badge.svg)](https://github.com/ph-rast/dcnet/actions) -->
<!-- <\!-- badges: end -\-> -->

## Installation

This package is in development. The current version can be insalled with

    if (!requireNamespace("remotes")) { 
      install.packages("remotes")   
    }   
    remotes::install_github("ph-rast/dcnet")
    #> Downloading GitHub repo ph-rast/dcnet@HEAD
    #> ── R CMD build ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    #>      checking for file ‘/tmp/RtmpBZXlgw/remotes561524b8bf0f7/ph-rast-dcnet-8bedc93/DESCRIPTION’ ...  ✔  checking for file ‘/tmp/RtmpBZXlgw/remotes561524b8bf0f7/ph-rast-dcnet-8bedc93/DESCRIPTION’
    #>   ─  preparing ‘dcnet’:
    #>    checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
    #> ─  cleaning src
    #>   ─  installing the package to process help pages
    #>      Loading required namespace: dcnet
    #>      Loading required package: cmdstanr
    #>      This is cmdstanr version 0.8.1
    #>    - CmdStanR documentation and vignettes: mc-stan.org/cmdstanr
    #>    - CmdStan path: /home/philippe/.cmdstan/cmdstan-2.35.0
    #>    - CmdStan version: 2.35.0
    #>      Compiling Stan models...
    #>      Stan models compiled successfully.
    #>   ─  saving partial Rd database
    #>   ─  cleaning src
    #>   ─  checking for LF line-endings in source and make files and shell scripts
    #>   ─  checking for empty or unneeded directories
    #>        NB: this package now depends on R (>= 3.5.0)
    #>        WARNING: Added dependency on R >= 3.5.0 because serialized objects in
    #>      serialize/load version 3 cannot be read in older versions of R.
    #>      File(s) containing such objects:
    #>        ‘dcnet/data/ema-fitbit.rda’
    #> ─  building ‘dcnet_0.0.0.9000.tar.gz’
    #>      
    #> 
    #> Installing package into '/home/philippe/R/x86_64-pc-linux-gnu-library/4.4'
    #> (as 'lib' is unspecified)

## Example: Dynamic Conditional Networks with EMA Data

This section illustrates the use with a dataset that contains synthetic
data from an ecological momentary assessment. The data is derived from a
100 Day study conducted at UC Davis (O’Laughlin et al. 2020, Willams et
al. 2020).

After loading the package, the `ema_fitbit` data is available as a list
of 38 matrices. Each matrix contains the four variables of interest:
`active`, `excited`, `totalDistance`, and `interested`. The rows
represent 98 daily observations, and the columns correspond to these
four variables. Each matrix is for one person (N=38).

    library(dcnet)

    ## Head of the data matrix for the first person:
    head(ema_fitbit[[1]])
    #>        active  excited totalDistance interested
    #> [1,] 43.27581 43.99161      6.163058   49.14676
    #> [2,] 64.61238 69.30582      6.382989   89.82558
    #> [3,] 23.61651 39.73926      5.398924   65.60202
    #> [4,] 63.72430 34.37902      6.345286   62.79230
    #> [5,] 85.38356 43.15638      9.275960   70.69611
    #> [6,] 71.11017 14.57427      7.890428   59.73196

Illustrate first person’s data. Given that `ema_fitbit` is a list, we
first need to extract person 1 and add a time index for `ggplot`:

    library(ggplot2)

    df1 <- data.frame(ema_fitbit[[1]], time = 1:98)

    ## Manually reshape the data to long format using base R
    df1_long <- do.call(rbind, lapply(names(df1)[-5], \(var) {
      data.frame(time = df1$time, variable = var, value = df1[[var]])
    }))

    ## Convert the 'variable' column to a factor for consistent ordering in facets
    df1_long <- df1_long |>
      within({ variable <- factor(variable, levels = names(df1)[-5]) })

    ## Plot with ggplot2
    plt <- ggplot(df1_long, aes(x = time, y = value)) +
      geom_smooth(show.legend = FALSE, se = FALSE) +
      geom_point(alpha = 0.4) +
      facet_wrap(~ variable, scales = "free_y") +
      labs(x = "Time", y = "Value", title = "Time Series for Variables for Person 1") +
      theme_minimal()

    plt
    #> `geom_smooth()` using method = 'loess' and formula = 'y ~ x'

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

## Fit Dynamic Conditional Correlation (DCC) Network Model

Next we will fit the multilevel DCC model for the shocks, the residual
conditional covariances, and a multilevel VAR structure for means
structure.

The function `dcnet` reads the data in the structure given in
`ema_fitbit`. The `parameterization = "DCCr"` selects the DCC with
random effects on all components, the garch parameters for the standard
deviations as well as all garch parameters for the conditional
correlations. In order to shorten computation time, we will use
`sampling_algorithm = 'variational'`. `variational` employs Stan’s
variational inference algorithm that approximates the posterior via the
“meanfield” algorithm. For full sampling, one could use `hmc`, but this
option will need to run hours and days to converge. The `init = 0.5`
argument is optional, but works well to create good starting values that
don’t stray too far.

## Fitting the Multilevel DCC Model

Next, we will fit the **multilevel DCC model** to estimate the shocks,
the residual conditional covariances, and a **multilevel VAR structure**
for the mean.

The `dcnet()` function reads the data provided in the structure of
`ema_fitbit`. Setting the argument `parameterization = "DCCr"` selects
the **DCC model with random effects** on all components. This includes
random effects for both the **GARCH parameters** governing the standard
deviations and the **GARCH parameters** for the conditional
correlations.

To shorten computation time, we will use
`sampling_algorithm = "variational"`. The `variational` algorithm
employs **Stan’s variational inference** to approximate the posterior
distribution using the **meanfield** approximation. For full Bayesian
sampling, one could set `sampling_algorithm = "hmc"`. However, this
option is computationally expensive and may take several hours or even
days to converge.

The argument `init = 0.5` is optional but often useful, as it provides
reasonable starting values that help the algorithm converge faster
without straying too far from the posterior space.

    fit <- dcnet(data = ema_fitbit, parameterization = "DCCr",
                 sampling_algorithm = 'variational', init = 0.5)
    #> ------------------------------------------------------------ 
    #> EXPERIMENTAL ALGORITHM: 
    #>   This procedure has not been thoroughly tested and may be unstable 
    #>   or buggy. The interface is subject to change. 
    #> ------------------------------------------------------------ 
    #> Gradient evaluation took 0.037546 seconds 
    #> 1000 transitions using 10 leapfrog steps per transition would take 375.46 seconds. 
    #> Adjust your expectations accordingly! 
    #> Begin eta adaptation. 
    #> Iteration:   1 / 250 [  0%]  (Adaptation) 
    #> Iteration:  50 / 250 [ 20%]  (Adaptation) 
    #> Iteration: 100 / 250 [ 40%]  (Adaptation) 
    #> Iteration: 150 / 250 [ 60%]  (Adaptation) 
    #> Iteration: 200 / 250 [ 80%]  (Adaptation) 
    #> Iteration: 250 / 250 [100%]  (Adaptation) 
    #> Success! Found best value [eta = 0.1]. 
    #> Begin stochastic gradient ascent. 
    #>   iter             ELBO   delta_ELBO_mean   delta_ELBO_med   notes  
    #>    100   -328291365.241             1.000            1.000 
    #>    200    -12304634.727            13.340           25.680 
    #>    300    -10712594.041             8.943            1.000 
    #>    400    -12715202.222             6.747            1.000 
    #>    500     -1061938.068             7.592            1.000 
    #>    600     -1054785.358             6.328            1.000 
    #>    700      -405834.872             5.652            1.000 
    #>    800      -221394.906             5.050            1.000 
    #>    900      -130308.433             4.566            0.833 
    #>   1000      -107944.588             4.131            0.833 
    #>   1100       -90826.374             3.772            0.699   MAY BE DIVERGING... INSPECT ELBO 
    #>   1200       -84699.818             3.464            0.699   MAY BE DIVERGING... INSPECT ELBO 
    #>   1300       -78582.181             3.203            0.207   MAY BE DIVERGING... INSPECT ELBO 
    #>   1400       -75811.563             2.977            0.207   MAY BE DIVERGING... INSPECT ELBO 
    #>   1500       -72869.152             2.781            0.188   MAY BE DIVERGING... INSPECT ELBO 
    #>   1600       -70536.552             2.610            0.188   MAY BE DIVERGING... INSPECT ELBO 
    #>   1700       -68544.729             2.458            0.157   MAY BE DIVERGING... INSPECT ELBO 
    #>   1800       -67247.556             2.322            0.157   MAY BE DIVERGING... INSPECT ELBO 
    #>   1900       -66112.193             2.201            0.149   MAY BE DIVERGING... INSPECT ELBO 
    #>   2000       -64840.390             2.092            0.149   MAY BE DIVERGING... INSPECT ELBO 
    #>   2100       -63803.124             1.993            0.078   MAY BE DIVERGING... INSPECT ELBO 
    #>   2200       -62919.943             1.903            0.078   MAY BE DIVERGING... INSPECT ELBO 
    #>   2300       -62411.470             1.821            0.072   MAY BE DIVERGING... INSPECT ELBO 
    #>   2400       -61852.541             1.745            0.072   MAY BE DIVERGING... INSPECT ELBO 
    #>   2500       -61368.940             1.676            0.040   MAY BE DIVERGING... INSPECT ELBO 
    #>   2600       -60834.327             1.612            0.040   MAY BE DIVERGING... INSPECT ELBO 
    #>   2700       -60173.532             1.552            0.037   MAY BE DIVERGING... INSPECT ELBO 
    #>   2800       -59678.631             1.497            0.037   MAY BE DIVERGING... INSPECT ELBO 
    #>   2900       -59489.045             1.446            0.033   MAY BE DIVERGING... INSPECT ELBO 
    #>   3000       -59129.428             1.398            0.033   MAY BE DIVERGING... INSPECT ELBO 
    #>   3100       -58819.912             1.365            0.029   MAY BE DIVERGING... INSPECT ELBO 
    #>   3200       -58358.420             0.509            0.020   MAY BE DIVERGING... INSPECT ELBO 
    #>   3300       -58097.770             0.504            0.019   MAY BE DIVERGING... INSPECT ELBO 
    #>   3400       -57699.919             0.499            0.017 
    #>   3500       -57548.685             0.133            0.016 
    #>   3600       -57318.440             0.133            0.016 
    #>   3700       -57002.322             0.080            0.014 
    #>   3800       -56648.696             0.053            0.011 
    #>   3900       -56462.389             0.029            0.009   MEDIAN ELBO CONVERGED 
    #> Drawing a sample of size 1000 from the approximate posterior...  
    #> COMPLETED. 
    #> Finished in  213.4 seconds.

Once the model is fit, we can request a summary (check your RAM – this
step might kill R if you are close to maxing out your memory):

    summary(fit)
    #> Model: mlVAR-DCCr
    #> Basic Specification: H_t = D_t R D_t
    #>  diag(D_t) = sqrt(h_ii,t) = c_h + a_h*y^2_[t-1] + b_h*h_[ii,t-1]
    #>  R_t = Q^[-1]_t Q_t Q^[-1]_t = ( 1 - a_q - b_q)S + a_q(u_[t-1]u'_[t-1]) + b_q(Q_[t-1])
    #> 
    #> Sampling Algorithm:  variational
    #> Distribution:  Gaussian
    #> ---
    #> Iterations:  30000
    #> Chains:  4
    #> Date:  Fri Oct 25 21:04:57 2024
    #> Elapsed time (min):  3.56
    #> 
    #> GARCH(1,1)  estimates for conditional variance on D:
    #> 
    #>            mean   sd  mdn 2.5% 97.5%
    #> c_h_tau_ac 0.61 0.05 0.61 0.52  0.72
    #> c_h_tau_ex 1.00 0.06 0.99 0.88  1.12
    #> c_h_tau_tD 0.38 0.05 0.38 0.29  0.48
    #> c_h_tau_in 0.48 0.05 0.49 0.39  0.59
    #> a_h_tau_ac 0.29 0.05 0.29 0.21  0.39
    #> a_h_tau_ex 0.28 0.05 0.28 0.20  0.39
    #> a_h_tau_tD 0.17 0.04 0.16 0.10  0.25
    #> a_h_tau_in 0.07 0.02 0.07 0.04  0.12
    #> b_h_tau_ac 0.37 0.12 0.35 0.19  0.64
    #> b_h_tau_ex 0.75 0.22 0.71 0.41  1.28
    #> b_h_tau_tD 0.30 0.10 0.29 0.15  0.52
    #> b_h_tau_in 0.30 0.08 0.28 0.17  0.49
    #> c_h_ac     7.83 0.85 7.79 6.28  9.57
    #> c_h_ex     4.88 0.56 4.83 3.95  6.13
    #> c_h_tD     1.80 0.21 1.78 1.41  2.22
    #> c_h_in     9.17 0.90 9.17 7.50 11.12
    #> a_h_ac     0.19 0.01 0.19 0.16  0.21
    #> a_h_ex     0.19 0.01 0.19 0.17  0.22
    #> a_h_tD     0.32 0.02 0.31 0.27  0.37
    #> a_h_in     0.15 0.01 0.15 0.13  0.17
    #> b_h_ac     0.80 0.01 0.80 0.77  0.82
    #> b_h_ex     0.78 0.02 0.78 0.75  0.81
    #> b_h_tD     0.62 0.03 0.62 0.56  0.67
    #> b_h_in     0.81 0.01 0.81 0.79  0.84
    #> 
    #> 
    #> GARCH(1,1) estimates for conditional variance on Q:
    #> 
    #>             mean   sd  mdn 2.5% 97.5%
    #> l_a_q_sigma 0.53 0.06 0.53 0.42  0.64
    #> l_b_q_sigma 1.30 0.12 1.29 1.08  1.57
    #> a_q_fixed   0.12 0.01 0.12 0.10  0.13
    #> b_q_fixed   0.77 0.02 0.77 0.74  0.80
    #> 
    #> 
    #> Unconditional correlation 'S' in Q:
    #> 
    #>         mean   sd  mdn 2.5% 97.5%
    #> S_ex-ac 0.48 0.05 0.49 0.38  0.57
    #> S_tD-ac 0.43 0.06 0.43 0.30  0.55
    #> S_in-ac 0.58 0.04 0.58 0.50  0.65
    #> S_tD-ex 0.22 0.07 0.22 0.09  0.35
    #> S_in-ex 0.65 0.04 0.65 0.57  0.72
    #> S_in-tD 0.31 0.07 0.31 0.18  0.43
    #> 
    #> 
    #> VAR(1) estimates on the location:
    #> 
    #>                mean   sd  mdn  2.5% 97.5%
    #> phi0_fixed_ac  1.42 0.74 1.42 -0.04  2.84
    #> phi0_fixed_ex  1.50 0.74 1.50  0.00  2.94
    #> phi0_fixed_tD  0.57 0.59 0.60 -0.57  1.76
    #> phi0_fixed_in  1.58 0.76 1.57  0.10  3.07
    #> phi0_tau_ac    0.50 0.19 0.47  0.22  0.97
    #> phi0_tau_ex    0.21 0.10 0.19  0.08  0.45
    #> phi0_tau_tD    0.50 0.16 0.47  0.25  0.89
    #> phi0_tau_in    0.70 0.27 0.65  0.30  1.32
    #> phi_fixed_acac 0.31 0.04 0.31  0.23  0.39
    #> phi_fixed_exac 0.17 0.04 0.17  0.09  0.25
    #> phi_fixed_tDac 0.04 0.02 0.04 -0.01  0.08
    #> phi_fixed_inac 0.17 0.03 0.17  0.10  0.23
    #> phi_fixed_acex 0.24 0.04 0.24  0.16  0.32
    #> phi_fixed_exex 0.36 0.05 0.36  0.26  0.45
    #> phi_fixed_tDex 0.02 0.03 0.02 -0.03  0.07
    #> phi_fixed_inex 0.30 0.03 0.30  0.24  0.36
    #> phi_fixed_actD 0.74 0.16 0.75  0.43  1.04
    #> phi_fixed_extD 0.49 0.18 0.49  0.15  0.83
    #> phi_fixed_tDtD 0.17 0.11 0.17 -0.05  0.40
    #> phi_fixed_intD 0.42 0.13 0.43  0.18  0.67
    #> phi_fixed_acin 0.25 0.04 0.26  0.17  0.34
    #> phi_fixed_exin 0.32 0.04 0.32  0.25  0.41
    #> phi_fixed_tDin 0.05 0.03 0.05  0.00  0.10
    #> phi_fixed_inin 0.42 0.04 0.42  0.35  0.48
    #> phi_tau_acac   0.10 0.01 0.10  0.07  0.13
    #> phi_tau_exac   0.05 0.01 0.05  0.03  0.07
    #> phi_tau_tDac   0.09 0.01 0.08  0.06  0.12
    #> phi_tau_inac   0.02 0.00 0.02  0.01  0.03
    #> phi_tau_acex   0.05 0.01 0.05  0.04  0.07
    #> phi_tau_exex   0.09 0.02 0.09  0.06  0.13
    #> phi_tau_tDex   0.04 0.01 0.04  0.03  0.05
    #> phi_tau_inex   0.06 0.01 0.06  0.04  0.08
    #> phi_tau_actD   0.15 0.03 0.15  0.09  0.23
    #> phi_tau_extD   0.07 0.02 0.07  0.04  0.13
    #> phi_tau_tDtD   0.12 0.03 0.12  0.07  0.19
    #> phi_tau_intD   0.08 0.02 0.08  0.05  0.13
    #> phi_tau_acin   0.07 0.01 0.07  0.05  0.10
    #> phi_tau_exin   0.05 0.01 0.04  0.03  0.07
    #> phi_tau_tDin   0.07 0.01 0.07  0.05  0.10
    #> phi_tau_inin   0.04 0.01 0.04  0.03  0.06
    #> 
    #> 
    #> Log density posterior estimate:
    #> 
    #>      mean        sd       mdn      2.5%     97.5% 
    #> -60953.68    737.40 -60787.05 -62760.11 -59977.53

## References

-   O’Laughlin, K. D., Liu, S., & Ferrer, E. (2020). Use of Composites
    in Analysis of Individual Time Series: Implications for
    Person-Specific Dynamic Parameters. Multivariate Behavioral
    Research, 1–18. <https://doi.org/10.1080/00273171.2020.1716673>
-   Williams, D. R., Martin, S. R., Liu, S., & Rast, P. (2020). Bayesian
    Multivariate Mixed-Effects Location Scale Modeling of Longitudinal
    Relations Among Affective Traits, States, and Physical Activity.
    European Journal of Psychological Assessment, 36 (6), 981–997.
    <https://doi.org/10.1027/1015-5759/a000624>
