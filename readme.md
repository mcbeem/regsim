`regsim` is an `R` package to make simulation quick and easy. It can be installed from github with the following `R` code.

``` r
# the following two lines need only to be run once 
install.packages("devtools")
devtools::install_github("mcbeem/regsim")

# add regsim's functions to R namespace so they can be used
library(regsim)
```

Example 1: Confounders and mediators
====================================

A researcher envisions the following true data-generating process:

![](readme_files/figure-markdown_github/unnamed-chunk-4-1.png)

And is interested in estimating the effect of *X* on *Y* using a linear regression model of the following form:

*Y*<sub>*i*</sub> = *b*<sub>0</sub> + *b*<sub>1</sub>(*X*<sub>*i*</sub>)+...+*e*<sub>*i*</sub>

where ... potentially contains the effect of a set of *covariates* or *adjustment variables*. In this case, those covariates could include any combination of *A*, *B*, or *M*. What are the consequences of different choices? Simulation is a useful technique for exploring the results of analytic decisions. One reason it is so useful is that the true model, and the correct effect of *X* → *Y*, is known. Thus, the estimates produced by linear regression can be compared to the known, true values.

Defining the Data Generating (True) Model
-----------------------------------------

The process of specifying the data generating model is simple.

1.  **Sketch the path diagram (or SEM) describing the true model**. The process begins by sketching out a path diagram of the relationships between the variables. The arrows indicate not only causation, but in this case also imply a linear functional form. Each arrow must be assigned a *path coefficient*; a regression slope describing how much the downstream variable will change in response to an isolated-one unit change of the upstream variable.

2.  **Classify the variables**. Next, the variables in the diagram are classified as either **exogenous** or **endogenous**. An exogenous variable has no arrow directed towards it; it has no in-system causes. An endogeous variable has at least one arrow directed toward it. In our example, variables *A* and *B* are exogenous, while *X*, *M*, and *Y* are endogenous.

3.  **Write models for the endogenous variables**. Each exogenous variable must be described in terms of its causes. The syntax is similar to `R`'s model formula used in `lm()` and many other model-fitting routines. In our example path diagram, *X*, *M*, and *Y* were endogenous, so a model for each one will need to be provided.

Let's begin by examining *X*, which is caused by *A* and *B* according to our diagram. In `R`'s formula syntax, we would write:

``` r
X ~ A + B
```

Where the `~` character is read "predicted by" or "regressed on." This formula corresponds to a regression model of the form:

*X*<sub>*i*</sub> = *b*<sub>0</sub> + *b*<sub>1</sub>(*A*<sub>*i*</sub>)+*b*<sub>2</sub>(*B*<sub>*i*</sub>)+*e*<sub>*i*</sub>

Since we are specifying the true model, we must provide values for parameters *b*<sub>1</sub> and *b*<sub>2</sub>. The intercept, *b*<sub>0</sub>, and the residual variance *e*<sub>*i*</sub> are nuisance parameters in many simulation application. We do not need to choose values for them when using `regsim()` to simulate data. The `regsim()` function uses the `simulateSEM()` function from the `lavaan` package to do the data generation.

Let's choose values of 0.5 for both *b*<sub>1</sub> and *b*<sub>2</sub>. The model for *X*, expressed in the proper syntax for `regsim()`, is as follows.

``` r
X ~ .5*A + .5*B
```

We must do the same for the endogenous variables *M* and *Y*. Variable *M* has only one arrow pointing to it, so its model has only one predictor variable, *X*. Choosing 0.5 as the path coefficient for *X* → *M*, *M*'s model is written:

``` r
M ~ .5*X
```

Finally, *Y* is caused by *A*, *B*, and *M*, so its model must include these variables. Sticking with 0.5 as our arbitrary value for all the path coefficients, the model for *Y* is written:

``` r
Y ~ .5*A + .5*B + .5*M
```

We will create an object called `true.model` that will be passed to `regsim()` to describe the model. It will consist of a quoted string including all three of the models specified above. These must be separated with line breaks. The full `lavaan` syntax is available for data generation, including latent variables.

``` r
true.model <- "X ~ .5*A + .5*B
               M ~ .5*X
               Y ~ .5*M + .5*A + .5*B"
```

By default, `regsim()` will sample the exogenous as well as the residuals for the endogenous variables from independent normal distributions. You can alter the distribution, to some extent, via the optional `kurtosis=` and `skewness=` arguments.

Defining the Fitted Model
-------------------------

The fitted model is the regression model that will be fitted to the simulated data sampled from the true model. This model is specified using `R` formula syntax. Here, we will specify a simple regression of *Y* on *X* with no covariates, and will save it an object called `fit.model`. **Note**: the fitted model's syntax must be a quoted string.

``` r
fit.model <- "Y ~ X"
```

Setting the other simulation parameters
---------------------------------------

We must make a few final decisions before running `regsim()`. They are:

-   `reps`: The number of times that the simulation-estimation process will be repeated. Here there is a tradeoff between the accuracy of the results and computation time. I suggest at least 1,000 reps, hopefully more, before interpreting the results too seriously.

-   `n`: The sample size generated for each rep.

-   `targetparm`: Which parameter's estimate and standard error should be scrutinized? This is expressed as a quoted string. It must be a variable from the right-hand side of the fitted model.

-   `targetval`: The true value of the relationship between the target parameter and the response variable as specified in the data generating model. In the path diagram considered in this example, the effect of *X* on *Y* is completely mediated by *M*.

By path tracing rules, the total effect of one variable on another when the effect is mediated is equal to the direct effect plus the indirect (mediated) effect. In this case, the direct effect is zero, because there is no direct *X* → *Y* path in the data generating model. The indirect effect is the product of all the path coefficients in the chain *X* → *M* → *Y*. In this case, the indirect effect is .5 × .5 = .25. We will use this value as our target value for *X*'s effect on *Y*.

Doing the simulation
--------------------

The simulation is performed using the `regsim()` function. The `true.model` and `fit.model` objects, previously created, are passed as arguments. The simulation results are assigned to the `result_1a` object. A random number seed has been specified for reproducibility.

Note that, by default, `regsim()` does not impose a standardized metric on the generated variables, so users do not have to worry about implying variable correlations greater than one. In other words, when `standardized=FALSE`, path coefficients can be freely chosen.

``` r
set.seed(123)

result_1a <- regsim(reps=1000, n=100, true.model=true.model, 
                 fit.model=fit.model, targetparm="X",
                 targetval=.25)
```

Before examining the results, let's verify that we specified the data generating model correctly. Running `plot()` on the `regsim()` output with the argument `plot="model"` will produce a path diagram corresponding to this model.

``` r
plot(result_1a, type="model")
```

![](readme_files/figure-markdown_github/unnamed-chunk-12-1.png)

This plot looks correct, but the path coefficients for *A* → *Y* and *B* → *X* are superimposed. We can try an alternate layout for the graph. See `?lavaan::semPaths` for details on the options (under 'layout'). In this case, the `"spring"` layout works better and avoids overplotting.

``` r
plot(result_1a, type="model", layout="spring")
```

![](readme_files/figure-markdown_github/unnamed-chunk-13-1.png)

Having verified that the data-generating model was set up correctly, let's examine the `regsim()` output.

``` r
result_1a
```

    ## $targetval
    ## [1] 0.25
    ## 
    ## $expected.b
    ## [1] 0.5843389
    ## 
    ## $bias
    ## [1] 0.3343389
    ## 
    ## $analytic.SE
    ## [1] 0.10401
    ## 
    ## $empirical.SE
    ## [1] 0.09922978
    ## 
    ## $analytic.CI
    ##      2.5%     97.5% 
    ## 0.3779346 0.7907433 
    ## 
    ## $empirical.CI
    ##      2.5%     97.5% 
    ## 0.3788186 0.7775613 
    ## 
    ## $coverage
    ## [1] 0.104
    ## 
    ## $RMSE
    ## [1] 0.3487395
    ## 
    ## $adjustment.sets
    ##  { A, B }

-   `$targetval` The true value of the target parameter, *β*, which was specified when the function was called. It is used to calculate bias, RMSE, and coverage.

-   `$expected.b` The mean value of the target parameter estimate across all the simulated repetitions, *E*(*b*). This value should closely approximate the corresponding true value.

-   `$bias` Bias is calculated as the mean of the difference between the true value of the target parameter and its estimate across repetitions, *E*(*b* − *β*). This value should be close to zero when the model is correctly specified.

-   `analytical.SE` The average estimated standard error for the target parameter across the repetitions. The analytical standard error should closely approximate the empirical standard error if the model is correctly specified. A violation of certain regression assumptions can cause it to diverge from the empirical standard error.

-   `empirical.SE` The calculated standard deviation of the parameter estimates, *b*. The empirical standard error is an estimated of what the sampling variability of the parameter actually is.

-   `$analytic.CI` The mean lower and upper analytic confidence interval boundaries for the target parameter across repetitions based on the limits specified by the `interval=` argument. These should closely approximate the empirical confidence interval boundaries; when this does not occur, an assumption violation has likely occurred.

-   `$empirical.CI` The empirical percentiles of the estimated values of the target parameter at the limits specified by the `interval=` argument. By default these will be the 2.5th and 97.5th percentiles, corresponding to the 95% empirical confidence interval.

-   `$coverage` The proportion of repetitions in which the true value of the target parameter, *β*, is contained in the analytic confidence interval. The coverage rate should should approximate the confidence level for the interval if model assumptions are met.

-   `$RMSE` The root mean squared error, calculated as $\\sqrt{E\[(b - \\beta)^2\]}$. The RMSE is includes contributions from both sampling error and bias and is a single-value summary of how close, on average, estimates come to the true value.

-   `$adjustment.sets` The set(s) of covariates that, if included in the fitted regression model (argument `fit.model=`), would allow for unbiased estimation of the target parameter. The adjustment sets are calculated using the `adjustmentSets()` function of the `dagitty` package. When the set is empty, no covariates are required.

The `regsim()` output includes a few other components which are not printed but are nonetheless available. The full contents can be inspected by running `str()` on the `regsim()` output. The other components of output include:

-   `$true.model` The data generating (true) model.

-   `$fit.model` The fitted model.

-   `$b` The vector of parameter estimates across repetitions.

-   `$data` A simulated dataset from the first repetition. (The data for the other repetitions is not available).

-   `$true.DAG` The data generating model expressed as a DAG of class `dagitty`.

### Interpreting the results

The results indicate that the fitted model is producing a biased estimated of the true effect of *X* → *Y*. This is can clearly be observed in the discrepancy between the mean estimated value (`$expected.b`) and the true value, which is also represented in the `$bias` value. The empirical and analytic standard error are quite close to one another, as are the confidence intervals, indicating that the problem with this model is in its fixed (structural) component. The bias results in a coverage rate that is far below its 95% nominal value. The RMSE is much larger than the standard error as a result of the bias.

Because the fitted model is misspecified, the estimated effect of *X* → *Y* is biased. `regsim()` returns the parameter estimates across all of the repetitions in a list component `$b`. Below is a histogram with a superimposed density plot of the parameter estimates. The true value of the relationship is indicated with a dashed vertical red line. The locations of the 2.5th and 97.5th percentiles (representing the empirical 95% confidence interval) are indicated with dotted vertical lines.

This plot can be generated by running `plot()` on the `regsim()` output with `type="perf"`. (The optional `breaks=30` argument increases the number of bars in the plot).

``` r
plot(result_1a, type="perf", breaks=30)
```

![](readme_files/figure-markdown_github/unnamed-chunk-15-1.png)

We can visualize the distribution of the simulated data in lots of ways. The easiest is to `plot()` the `regsim()` output using `type="data"`. This plots the data from the first repetition of the simulation with the `pairs.panels()` function from the `psych` package.

``` r
plot(result_1a, type="data")
```

![](readme_files/figure-markdown_github/unnamed-chunk-16-1.png)

### Obtaining an Unbiased Estimate of the Effect of X on Y

The `$adjustment.sets` output indicates the source of the problem: variables *A* and *B* must be included in this model as covariates in order to correctly estimate the *X* → *Y* effect. Let's correctly specify this model and re-run the simulation.

The `fit.model=` argument of the code is altered to include *A* and *B* in the model formula.

``` r
result_1b <- regsim(reps=1000, n=100, true.model=true.model, 
                 fit.model="Y~X+A+B", targetparm="X",
                 targetval=.25)

result_1b
```

    ## $targetval
    ## [1] 0.25
    ## 
    ## $expected.b
    ## [1] 0.2435049
    ## 
    ## $bias
    ## [1] -0.006495071
    ## 
    ## $analytic.SE
    ## [1] 0.1133802
    ## 
    ## $empirical.SE
    ## [1] 0.1167926
    ## 
    ## $analytic.CI
    ##       2.5%      97.5% 
    ## 0.01844698 0.46856288 
    ## 
    ## $empirical.CI
    ##       2.5%      97.5% 
    ## 0.02176829 0.47546433 
    ## 
    ## $coverage
    ## [1] 0.942
    ## 
    ## $RMSE
    ## [1] 0.1169148
    ## 
    ## $adjustment.sets
    ##  { A, B }

The results indicate good performance. Bias is very low, the coverage rate is nearly 95%, there is close correspondence between the analytic and empirical standard error and confidence intervals, and the RMSE is dominated by sampling error.

An updated histogram / density plot of the parameter estimates illustrate that they are now centered on the true value of the parameter.

``` r
plot(result_1b, type="perf", breaks=30)
```

![](readme_files/figure-markdown_github/unnamed-chunk-18-1.png)

### Misspecifying the Model by Conditioning on the Mediator

In general, downstream descendents of the target parameter should not be included as covariates in an analysis. The only exception would be when explicit mediation modeling is intended. What happens if *M* is added to the model as a covariate, in addition to *A* and *B*? Simulation can answer that question.

``` r
result_1c <- regsim(reps=1000, n=100, true.model=true.model, 
                 fit.model="Y~X+A+B+M", targetparm="X",
                 targetval=.25)

result_1c
```

    ## $targetval
    ## [1] 0.25
    ## 
    ## $expected.b
    ## [1] 0.0003883205
    ## 
    ## $bias
    ## [1] -0.2496117
    ## 
    ## $analytic.SE
    ## [1] 0.114455
    ## 
    ## $empirical.SE
    ## [1] 0.1122119
    ## 
    ## $analytic.CI
    ##       2.5%      97.5% 
    ## -0.2268336  0.2276103 
    ## 
    ## $empirical.CI
    ##       2.5%      97.5% 
    ## -0.2231423  0.2207332 
    ## 
    ## $coverage
    ## [1] 0.417
    ## 
    ## $RMSE
    ## [1] 0.2736511
    ## 
    ## $adjustment.sets
    ##  { A, B }

As the results indicate, conditioning on the mediator is a bad idea. It results in an intensely biased estimate of the effect of *X* on *Y*. In fact, the average estimated effect across all the repetitions is now zero. Note that the variable *M* is not included in the `$adjustment.sets` component of the output.

``` r
plot(result_1c, type="perf", breaks=30)
```

![](readme_files/figure-markdown_github/unnamed-chunk-20-1.png)

Example 2: Violating distributional assumptions
===============================================

The `regsim()` function can generate non-normal data via the optional `skewness=` and `kurtosis=` arguments. (These are passed to `lavaan`'s `simulateData` function via `...`).

In this case the model will be a simple regression of *Y* on *X*, but the variables will be generated tuch they they are skewed and kurtotic. The `skewness=` and `kurtosis=` arguments expect a vector of these values for each variable in the data-generating model. However, I do not understand how the underlying `simulateData()` function orders the variables. Thus, I find it necessary to experiment with the order, each time plotting the the first generated dataset (via the `$data` output component) until the desired result is obtained.

I intend to generate data such that *Y* is severely skewed and *X* is not. I will try specifying `skewness=c(4,0)`, hoping that the skewness will be applied to *Y*. Also, it is important to note that skewness and kurtosis are not independent; skewness implies kurtosis. If *κ* is the full kurtosis of a distribution, and *γ* is the skewness, then the following inequality must be satisfied.

*κ* ≥ 1 + *γ*<sup>2</sup>

Non-normal data generation in `regsim()` is based on the algorithm presented by Vale and Maurelli (1983, as implemented in `lavaan`'s `simulateData()` function) and is quite sensitive to the kurtosis value. Successful convergence appears to require kurtosis values well above those that would satistify the inequality given above. Through trial and error I discovered that a kurtosis value of 26 or higher is needed to simulate data with a skewness of ±4. This is an extreme level of non-normality.

``` r
result_2 <- regsim(reps=1000, n=500, true.model="Y~.5*X", 
                   fit.model="Y~X", targetparm="X", targetval=.5, 
                   kurtosis=c(26, 0), skewness=c(4, 0))
```

If the provided kurtosis value is too small for the skewness, you will see the following error message.

    ## Warning in fleishman1978_abcd(skewness = SK[i], kurtosis = KU[i]): lavaan
    ## WARNING: ValeMaurelli1983 method did not convergence, or it did not find
    ## the roots

In this case, you should not trust `regsim()`'s output!

Next, I examine the `$data` with `psych`'s `pairs.panels()` function.

``` r
plot(result_2, type="data")
```

![](readme_files/figure-markdown_github/unnamed-chunk-23-1.png)

The plot show that the *Y* variable was skewed as intended.

Now let's examine the simulation results.

``` r
result_2
```

    ## $targetval
    ## [1] 0.5
    ## 
    ## $expected.b
    ## [1] 0.4997811
    ## 
    ## $bias
    ## [1] -0.0002188914
    ## 
    ## $analytic.SE
    ## [1] 0.04450613
    ## 
    ## $empirical.SE
    ## [1] 0.08573579
    ## 
    ## $analytic.CI
    ##      2.5%     97.5% 
    ## 0.4123382 0.5872240 
    ## 
    ## $empirical.CI
    ##      2.5%     97.5% 
    ## 0.3393903 0.6701404 
    ## 
    ## $coverage
    ## [1] 0.696
    ## 
    ## $RMSE
    ## [1] 0.08569319
    ## 
    ## $adjustment.sets
    ##  {}

The results are unbiased, but the analytic standard errors are wrong. In fact, it is about half the magnitude that it should be. Consequently, the analytic 95% confidence interval is too narrow, and the coverage rate is only about 70%. A hypothesis test for this parameter would be *liberal*, meaning that its long-run false positive rate would be much higher than the alpha criterion should allow.

The regression assumption that has been violated in this case is that the distribution of the residuals does not follow a normal distribution. The strutural part of the model is fine. The incorrect residual distribution manifests as an improperly-shaped shaped sampling distribution for the parameter estimates.

Below is a density plot of the empirical sampling distribution. The analytic sampling distribution is shown as the blue dotted curve. This plot can be obtained by running `plot()` on the `regsim()` output with `type="compareSE"`. When the regreesion assumptions are met, these two distributions are identical in expectation. But they diverge when assumptions are violated. The true value of the parameter is denoted by the red vertical line.

``` r
plot(result_2, type="compareSE")
```

![](readme_files/figure-markdown_github/unnamed-chunk-25-1.png)

Not only is the empirical sampling distribution platykurtotic, it is also slightly skewed to the right.

This is a rather extreme example of nonnormality. Linear regression works reasonbly well for minor departures from the normality assumption.

Example 3: Measurement error
============================

Measurement error is a ubiquitous feature of real data. Linear regression models assume that the predictor variables are measured without error. No such assumption is made about the response variable. We can explore this issue with simulation.

The **reliability** of a variable describes the proportion of its variance that is true score rather than measurement error. Many assessments in psychology have reliability coefficients of 0.7 to 0.9. Single item indicators often have reliabilitiy coefficients in the 0.5 range.

Measurement error in predictor variables
----------------------------------------

Let us imagine that we can measure variable *X* with a reliability of 0.8. According to classical test theory, the observed score *X*<sub>*O*</sub> is composed of a true score component *X*<sub>*T*</sub> plus measurement error component *X*<sub>*E*</sub>.

This situation can be represented by the following figure.

![](readme_files/figure-markdown_github/unnamed-chunk-26-1.png)

Two key equations in classical test theory are

*V**a**r*(*X*)=*V**a**r*(*T*)+*V**a**r*(*E*)\]

and

*ρ*<sub>*XX*</sub> = *V**a**r*(*T*)/\[*V**a**r*(*T*)+*V**a**r*(*E*)\]

where *ρ*<sub>*XX*</sub> is the reliability coefficient, *X* is the observed score, *T* is the true score, and *E* is measurement error. These equations provide guidance for setting the parameters of the data generating model in order to produce an observed *X* variable with the desired reliability. Let us assume that the exogenous latent true score *X*<sub>*T*</sub> exists on a standardized metric with unit variance. Let us set the loading relating the observed score *X*<sub>*O*</sub> to the true score to be one. Then the true score contribution to *X*<sub>*O*</sub>'s variance is 1.0 (because of the rule *V**a**r*(*bX*)=*b*<sup>2</sup>*V**a**r*(*X*)). The variance of the observed score must exceed the true score variance because the measurement error variance is additive. A bit of algebra allows us to solve for the residual variance required to produce *X*<sub>*O*</sub> with the desired reliability.

*V**a**r*(*E*)=(1/*ρ*<sub>*XX*</sub>)−1

In this case, to produce a reliability of 0.8, the error variance must be set to 0.25.

The following diagram represents the data-generating model for this situation.

![](readme_files/figure-markdown_github/unnamed-chunk-27-1.png)

Note that the residual variance for *X*<sub>*O*</sub> has been set to 0.25, which represents the error variance. Its its true score variance is 1, and therefore its total variance is therefore 1.25. Finally, 1/1.25 = .80, the desired reliability coefficient.

Note that there is no path from *X*<sub>*O*</sub> to *Y*. *Y* is caused by the true value, *X*<sub>*T*</sub>, and not by the observed *X*<sub>*O*</sub>. In fact, the true score is a confounder. The true value of *X*<sub>*O*</sub> → *Y* is actually zero!

The code to define the data generating model is as follows.

``` r
true.model <- "X_T =~ 1*X_O  
               X_T ~~ 1*X_T   
               X_O ~~ .25*X_O
               Y ~ .5*X_T"
```

In `lavaan` syntax, the symbol `=~` is read "measured by" and is used to define latent variables in terms of observed indicators. In this case, the true score is latent and is measured by the observed score *X*<sub>*O*</sub> with its loading fixed to one. The symbol `~~` is used to specify variances or residual variances. In the syntax above, *X*<sub>*T*</sub>'s total variance is set to one, and *X*<sub>*O*</sub>'s residual variance is set to 0.25. Finally, `~` represents path coefficients (read "regressed on" or "caused by") as before.

Next, we run the simulation. Also note that the target value I have set is equal to the effect of *X*<sub>*T*</sub> → *Y*, not the zero effect of *X*<sub>*O*</sub> → *Y*. Because *X*<sub>*T*</sub> is defined as a latent variable, it will not appear in the simulated datasets produced by `regsim()`.

``` r
result_3a <- regsim(reps=1000, n=300, true.model=true.model, 
                    fit.model="Y~X_O", targetparm="X_O", targetval=.5)
result_3a
```

    ## $targetval
    ## [1] 0.5
    ## 
    ## $expected.b
    ## [1] 0.4023442
    ## 
    ## $bias
    ## [1] -0.09765579
    ## 
    ## $analytic.SE
    ## [1] 0.05314453
    ## 
    ## $empirical.SE
    ## [1] 0.05083078
    ## 
    ## $analytic.CI
    ##      2.5%     97.5% 
    ## 0.2977581 0.5069303 
    ## 
    ## $empirical.CI
    ##      2.5%     97.5% 
    ## 0.3049644 0.5061742 
    ## 
    ## $coverage
    ## [1] 0.555
    ## 
    ## $RMSE
    ## [1] 0.1100811
    ## 
    ## $adjustment.sets

The results indicate that the estimated effect of *X*<sub>*O*</sub> → *Y* is biased. Linear regression models assume that the predictor variables are measured without error. Thus, we cannot obtain an unbiased estimate of the effect of *X* → *Y* when *X* is measured with error. In general, predictor variables cannot be replaced by noisy proxies without harming inferences. The bias can be visualized by plotting the `regsim()` output with `type="perf"`.

``` r
plot(result_3a, type="perf", breaks=30)
```

![](readme_files/figure-markdown_github/unnamed-chunk-30-1.png)

Noisy proxies of predictor variables are not ideal, but they are generally best that we can achieve, because all acts of measurement involve measurement error. The magnitude of the bias is directly related to the measurement reliability of the predictor variable. Thus, it is crucial to measure predictor variables as precisely as possible. Measurement error bias is no different than confounding bias. Both forms are caused by omitted variables.

Measurement error in the response variable
------------------------------------------

Linear regression models do not assume that the response variable is measured without error. We can modify our previous analysis to illustrate what happens when *Y*, rather than *X*, is measured with error. The new data generating model looks like this:

![](readme_files/figure-markdown_github/unnamed-chunk-31-1.png)

It can be defined as follows.

``` r
true.model <- "Y_T =~ 1*Y_O  
               Y_T ~~ 1*Y_T   
               Y_O ~~ .25*Y_O
               Y_T ~ .5*X"
```

The `regsim()` output confirms that the results are unbiased even though *Y* is measured with error.

``` r
result_3b <- regsim(reps=1000, n=300, true.model=true.model, 
                    fit.model="Y_O~X", targetparm="X", targetval=.5)
result_3b
```

    ## $targetval
    ## [1] 0.5
    ## 
    ## $expected.b
    ## [1] 0.4999318
    ## 
    ## $bias
    ## [1] -6.81728e-05
    ## 
    ## $analytic.SE
    ## [1] 0.0646636
    ## 
    ## $empirical.SE
    ## [1] 0.06623344
    ## 
    ## $analytic.CI
    ##      2.5%     97.5% 
    ## 0.3726767 0.6271870 
    ## 
    ## $empirical.CI
    ##      2.5%     97.5% 
    ## 0.3755090 0.6346035 
    ## 
    ## $coverage
    ## [1] 0.946
    ## 
    ## $RMSE
    ## [1] 0.06620036
    ## 
    ## $adjustment.sets
    ##  {}

Plotting the results illustrates the unbiased estimation.

``` r
plot(result_3b, type="perf", breaks=30)
```

![](readme_files/figure-markdown_github/unnamed-chunk-34-1.png)
