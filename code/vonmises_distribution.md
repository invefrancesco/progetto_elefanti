MLE Von Mises
================
Francesco Invernizzi
2025-06-26

# Von Mises distribution

Il post di [Whittenbury
Daniel](https://dlwhittenbury.github.io/ds-1-sampling-and-visualising-the-von-mises-distribution.html)
contiene alcune note riguardanti la distribuzione von Mises.

La “statistica direzionale” ha numerosi usi nelle scienze

- la statistica classica, come sviluppata da Gauss, assume che gli
  errori di misura siano “piccoli” e distribuiti in uno **spazio lineare
  infinito** (come la retta dei numeri reali o uno spazio euclideo)

- questo funziona bene quando le misure sono **molto precise** come nel
  caso di astronomia o dei geometria che usano strumenti accurati

- **le misure angolari** (come direzioni su una bussola, ore del giorno,
  orientamenti spaziali, ecc.) **vivono naturalmente su una struttura
  topologica diversa** (un cerchio o una sfera).

  - Se gli errori diventano grandi (cioè, se i dati sono molto
    dispersi), **non ha più senso ignorare la curvatura** dello spazio
    su cui si trovano i dati.

  - Per esempio, dire che la media tra 10° e 350° è 180° è assurdo in un
    contesto circolare, ma potrebbe sembrare corretto se si ignorasse la
    topologia.

La distribuzione di von Mises è stata introdotta da Richard von Mises
nel 1918 e la sua PDF è
$$ f(x | \mu, \kappa) =  \frac{\exp(k \cos(x - \mu))}{2 \pi I_0(k)}$$

- x è la variabile aleatoria angolare

- $\mu$ è la **direzione media** (con $\mu \in [0, 2 \pi)$)

- $\kappa$ è un parametro di concentrazione della distribuzione, che può
  assumere valori non negativi ($\kappa \in [0, \infty)$)

  - più è grande $\kappa$, più la distribuzione è concentrata attorno
    alla direzione media

- $I_0(\kappa)$ è la funzione di Bessel modificata del primo tipo e
  ordine zero

Le distribuzioni circolari hanno diversi vantaggi

- sono definite su uno spazio campionario corretto

  - **la normale può produrre stime errate vicino ai limiti (es. 0 e
    $2\pi$)**

- la distribuzione di von Mises ha un comportamento limite corretto
  quando $\kappa \to 0$: tende a una distribuzione uniforme sul cerchio
  dove tutte le direzioni sono equiprobabili

- **usare la normale in problemi circolari può portare a errori seri
  come fondere in modo scorretto una distribuzione a priori con una
  funzione di verosimiglianza**

  - **nelle misurazioni dell’angolo di arrivo (in applicazioni di
    tracciamento), si rischia di aggiornare male la stima bayesiana se
    si ignora la natura ciclica della variabile**

## Optim function

A common task in statistics is maximizing (or minimizing) complex
univariate or multivariate functions. This is typically done in the
context of maximizing likelihood functions with respect to a vector of
unknown parameters, given the observed data.

- More generally, this task pertains to an “optimization” problem, where
  we seek to finding a set of parameter values that minimizes or
  maximizes some pre-defined objective function.

- When this objective function is a log likelihood function, this
  reduces to the problem of maximum likelihood estimation.

- There are many potential approaches to use for optimization in this
  context.

- How best to optimize a function of interest depends on the nature of
  the function to be optimized, as well as practical concerns regarding
  a candidate procedure to be used.

By default `optim` performs minimization, but it will maximize if
control\$fnscale is negative. optimHess is an auxiliary function to
compute the Hessian at a later stage if hessian = TRUE was forgotten.

- The default method is an implementation of that of Nelder and Mead
  (1965), that uses only function values and is robust but relatively
  slow. It will work reasonably well for non-differentiable functions.

- Method “BFGS” is a quasi-Newton method (also known as a variable
  metric algorithm), specifically that published simultaneously in 1970
  by Broyden, Fletcher, Goldfarb and Shanno. This uses function values
  and gradients to build up a picture of the surface to be optimized.

- Method “CG” is a conjugate gradients method based on that by Fletcher
  and Reeves (1964) (but with the option of Polak–Ribiere or
  Beale–Sorenson updates). Conjugate gradient methods will generally be
  more fragile than the BFGS method, but as they do not store a matrix
  they may be successful in much larger optimization problems.

- Method “L-BFGS-B” is that of Byrd et al. (1995) which allows box
  constraints, that is each variable can be given a lower and/or upper
  bound. The initial value must satisfy the constraints. This uses a
  limited-memory modification of the BFGS quasi-Newton method. If
  non-trivial bounds are supplied, this method will be selected, with a
  warning.

- Nocedal and Wright (1999) is a comprehensive reference for the
  previous three methods.

- Method “SANN” is by default a variant of simulated annealing given in
  Belisle (1992).

  - Simulated-annealing belongs to the class of stochastic global
    optimization methods. It uses only function values but is relatively
    slow.

  - It will also work for non-differentiable functions. This
    implementation uses the Metropolis function for the acceptance
    probability.

  - By default the next candidate point is generated from a Gaussian
    Markov kernel with scale proportional to the actual temperature. If
    a function to generate a new candidate point is given, method “SANN”
    can also be used to solve combinatorial optimization problems.

  - Temperatures are decreased according to the logarithmic cooling
    schedule as given in Belisle (1992, p. 890); specifically, the
    temperature is set to
    `temp / log(((t-1) %/% tmax) * tmax + exp(1))`, where `t` is the
    current iteration step and `temp` and `tmax` are specifiable via
    control, see below.

  - Note that the “SANN” method depends critically on the settings of
    the control parameters. It is not a general-purpose method but can
    be very useful in getting to a good value on a very rough surface.

- Method “Brent” is for one-dimensional problems only, using
  `optimize(<ff>, lower, upper, tol = control$reltol)` where `<ff>` is
  `function(par) fn(par, ...)/control$fnscale`. It can be useful in
  cases when `optim()` is used inside other functions where only method
  can be specified, such as in `mle` from package stats4.

Function `fn` can return `NA` or `Inf` if the function cannot be
evaluated at the supplied value, but the initial value must have a
computable finite value of fn. (Except for method “L-BFGS-B” where the
values should always be finite.)

`optim` can be used recursively, and for a single parameter as well as
many. It also accepts a zero-length par, and just evaluates the function
with that argument \## Stima di una Von Mises con una covariata

1.  Creo due vettori $\mathbf{x}_n$ e $\mathbf{y}_n$ in cui

- $X_i \sim N(0, 1)$

- $Y_i \sim VonMises(\mu = 2atan(\beta_0 +\beta_1 \cdot x_i), k = (10))$

``` r
rndm_x <- rnorm(100)
rndm_y <- sapply(rndm_x, function(xi) circular::rvonmises(n = 1, mu = 2*atan(-1 + 3.5 * xi), kappa = 10))
```

2.  Definisco la funzione di log-verosimiglianza

``` r
log.lik.vm <- function(par, x, y){
  beta.0 <- par[1]
  beta.1 <- par[2]
  kappa <- exp(par[3])
  
  l <- array(dim = length(x))
  
  for(i in seq_along(x)){
    l[i] = kappa * cos(y[i] - 2 * atan(beta.0 + beta.1*x[i])) - log((besselI(kappa, nu = 0)))
  }
  
  return(-sum(l))
}
```

3.  Minimizzo la funione
    - `hessian = T` include nei risultati la matrice Hessiana, cioè la
      matrice di informazione osservata
    - Se questa matrice è invertibile, corrisponde asintoticamente alla
      matrice di varianza e covarianza dei parametri
    - La radice quadrata dei valori sulla diagonale principale,
      corrispondo agli standard error asintotici dei parametri

``` r
result <- optim(
  par = c(beta.0 = 0, beta.1 = 0, kappa = log(1)),
  fn = log.lik.vm,
  x = rndm_x,
  y = rndm_y, 
  hessian = T)

result$par
```

    ##    beta.0    beta.1     kappa 
    ## -1.031485  3.367014  2.459555

``` r
det(result$hessian) != 0 # condizione di invertibilità
```

    ## [1] TRUE

``` r
cov.matrix <- solve(result$hessian)
se <- sqrt(diag(cov.matrix))
se
```

    ##     beta.0     beta.1      kappa 
    ## 0.05104836 0.14179574 0.13804482

### Simulazione di $K$ repliche per lo stesso set di parametri veri

Definisco una funzione che dati

- `n` = la numerosità del campione casuale

- `K` = il numero di repliche

- `beta.0` e `beta.1` = il valore reale dei parametri

- `kappa` = il valore reale del parametro di concentrazione della Von
  Mises

Genera $K$ set di random $x$ e random $y$ e stima i parametri
ottimizzando la funzione di verosimiglianza per ogni set La seguente
funzione serve per generare i set di dati

``` r
gen.data <- function(n, beta.0, beta.1, kappa){
  x <- rnorm(n)
  y <- sapply(x, function(xi) circular::rvonmises(n= 1, mu=(beta.0 + beta.1 * xi), kappa = kappa))
  data <- tibble(x=x, y=y)
}
```

Per riproducibilità genero una sola volta la lista di df e la salvo in
un file locale. Il codice per generare la funzione non è incluso nel
file per poterlo compilare facilmente. - Inserisco come parametri
`n = 1000`, `beta.0 = -1`, `beta.1 = 3.5`, `kappa = 10` - Creo una lista
di 100 dataset

``` r
load(paste0(dir_data, "/simulazioni.RData"))

results <- map_dfr(
    .x = datasets,
    .f = ~ {
      mle <- optim(
        par = c(beta.0 = 0, beta.1 = 0, kappa = log(2)),
        fn = log.lik.vm,
        x = .x$x,
        y = .x$y
      )
      
      tibble(
        beta.0 = mle$par[1],
        beta.1 = mle$par[2],
        kappa  = exp(mle$par[3]),
        convergence = mle$convergence
      )
    }
  )

results
```

    ## # A tibble: 100 × 4
    ##    beta.0 beta.1    kappa convergence
    ##     <dbl>  <dbl>    <dbl>       <int>
    ##  1  5.55    3.73 3.91e- 1           0
    ##  2  5.50    3.68 3.36e- 1           0
    ##  3 13.3    46.0  1.99e-36           0
    ##  4 13.3    46.0  1.99e-36           0
    ##  5  7.03   19.2  2.14e-25           0
    ##  6 -0.907   3.21 1.18e+ 0           0
    ##  7 29.5     8.57 7.96e-24           0
    ##  8  5.50    3.65 3.17e- 1           0
    ##  9 13.3    46.0  1.99e-36           0
    ## 10 11.9    48.9  4.03e-40           0
    ## # ℹ 90 more rows

I risultati della simulazione sono decisamente da quelli che ci si
potrebbe aspettare. A titolo di esempio, seleziono le righe in cui il
parametro stimato per `beta.0` esce dall’intervallo $\pm 3 \sigma$ dove
$\sigma$ corrisponde al valore stimato tramite l’informazione di Fisher

``` r
print(
results %>% 
  filter(!between(beta.0, -1 - 3*se[1], -1 + 3*se[1])),
n = 50)
```

    ## # A tibble: 79 × 4
    ##    beta.0  beta.1    kappa convergence
    ##     <dbl>   <dbl>    <dbl>       <int>
    ##  1  5.55    3.73  3.91e- 1           0
    ##  2  5.50    3.68  3.36e- 1           0
    ##  3 13.3    46.0   1.99e-36           0
    ##  4 13.3    46.0   1.99e-36           0
    ##  5  7.03   19.2   2.14e-25           0
    ##  6 29.5     8.57  7.96e-24           0
    ##  7  5.50    3.65  3.17e- 1           0
    ##  8 13.3    46.0   1.99e-36           0
    ##  9 11.9    48.9   4.03e-40           0
    ## 10  7.92   31.8   4.30e-26           0
    ## 11 30.1   -10.4   6.84e- 2           0
    ## 12 13.3    46.0   1.99e-36           0
    ## 13  5.69    3.85  3.50e- 1           0
    ## 14  9.25    2.79  3.37e- 2           0
    ## 15 13.3    46.0   1.99e-36           0
    ## 16 13.3    46.0   1.99e-36           0
    ## 17  6.53    4.50  2.62e- 1           0
    ## 18  8.73   44.2   1.04e-30           0
    ## 19 10.8    16.2   8.93e-29           0
    ## 20 11.9    48.9   4.03e-40           0
    ## 21  6.42   32.1   1.60e-22           0
    ## 22 13.3    46.0   1.99e-36           0
    ## 23 13.3    46.0   1.99e-36           0
    ## 24 41.6    10.4   3.88e-33           0
    ## 25  8.73   44.2   1.04e-30           0
    ## 26  5.50    3.64  3.62e- 1           0
    ## 27 41.6    10.4   3.88e-33           0
    ## 28 29.5     8.57  7.96e-24           0
    ## 29 13.3    46.0   1.99e-36           0
    ## 30 11.9    48.9   4.03e-40           0
    ## 31 13.3    46.0   1.99e-36           0
    ## 32 10.8    16.2   8.93e-29           0
    ## 33 10.8    16.2   8.93e-29           0
    ## 34  7.03   19.2   2.14e-25           0
    ## 35  5.28    3.53  2.82e- 1           0
    ## 36  5.72    3.85  1.94e- 1           0
    ## 37  8.43   28.8   8.84e-25           0
    ## 38 10.8    16.2   8.93e-29           0
    ## 39 11.9    48.9   4.03e-40           0
    ## 40  8.73   44.2   1.04e-30           0
    ## 41 13.3    46.0   1.99e-36           0
    ## 42 13.3    46.0   1.99e-36           0
    ## 43 13.3    46.0   1.99e-36           0
    ## 44  5.63    3.78  2.43e- 1           0
    ## 45  5.81    3.93  3.75e- 1           0
    ## 46 29.5     8.57  7.96e-24           0
    ## 47 13.3    46.0   1.99e-36           0
    ## 48  5.83    3.93  3.44e- 1           0
    ## 49  0.181  -0.444 9.33e- 2           0
    ## 50 13.3    46.0   1.99e-36           0
    ## # ℹ 29 more rows

- Provo a visualizzare i campioni in cui i parametri sono anomali

``` r
ggplot(data = datasets[[1]]) + 
  geom_histogram(aes(x)) +
  theme_classic()
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](vonmises_distribution_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
circular::plot.circular(datasets[[1]]$y)
```

![](vonmises_distribution_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

``` r
ggplot(data = datasets[[1]]) +
  geom_point(aes(x, y)) + 
  theme_classic()
```

![](vonmises_distribution_files/figure-gfm/unnamed-chunk-7-3.png)<!-- -->

- Profile likelihood dei parametri

``` r
x <- datasets[[1]]$x
y <- datasets[[1]]$y

tibble(
  beta_0 = seq(from = -15, to = 15, length.out = 1000),
  neg_log_lik = map_dbl(beta_0, ~ log.lik.vm(par = c(.x, results$beta.1[1], results$kappa[1]), 
                                         x = x, 
                                         y = y
  ))) %>% 
  ggplot() +
  geom_line(aes(x = beta_0, y = neg_log_lik)) +
  theme_classic() +
  labs(
    x = expression(beta[0]),
    y = expression(-log * L(beta[0])),
    title = expression("Profilo di verosimiglianza per " * beta[0])
  ) 
```

![](vonmises_distribution_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
tibble(
  beta_1 = seq(from = -15, to = 15, length.out = 1000),
  neg_log_lik = map_dbl(beta_1, ~ log.lik.vm(par = c(results$beta.0[1], .x, results$kappa[1]), 
                                             x = x, 
                                             y = y
  ))) %>% 
  ggplot() +
  geom_line(aes(x = beta_1, y = neg_log_lik)) +
  theme_classic() +
  labs(
    x = expression(beta[1]),
    y = expression(-log * L(beta[1])),
    title = expression("Profilo di verosimiglianza per " * beta[1])
  ) 
```

![](vonmises_distribution_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
tibble(
  kappa = seq(from = -15, to = 35, length.out = 1000),
  neg_log_lik = map_dbl(kappa, ~ log.lik.vm(par = c(results$beta.0[1], results$beta.1[1], .x), 
                                             x = x, 
                                             y = y
  ))) %>% 
  ggplot() +
  geom_line(aes(x = kappa, y = neg_log_lik)) +
  theme_classic() +
  labs(
    x = expression(kappa),
    y = expression(-log * L(kappa)),
    title = expression("Profilo di verosimiglianza per " * kappa)
  ) 
```

![](vonmises_distribution_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->
