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
with that argument.

### Stima dei parametri della Von Mises

Siano $X_1, \dots, X_n$ `n` v.a. iid, con $X_i \sim$ VonMises($\mu, k$)

1.  Definisco la funzione di log-verosimiglianza delle $n$ v.a. come
    sommatoria delle singole log-verosimiglianze in virtù dell’ipotesi
    i.i.d.

``` r
log_lik_vm <- function(data, par){
  log_lik <- sum(circular::dvonmises(data, mu = par[1], kappa = par[2], log = T))
  return(-log_lik)
}
```

2.  Estraggo un campione casuale dalla popolazione VM($0, 10$)

``` r
data_rndm <- circular::rvonmises(n = 100, mu = 0, kappa = 10)
circular::plot.circular(data_rndm)
```

![](vonmises_distribution_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

3.  Minimizzo la funzione di -log-verosimiglianza (e controllo i
    risultati con mle.vonmises)

``` r
optim(
  par = c(0,0),
  fn = log_lik_vm,
  data = data_rndm,
)
```

    ## $par
    ## [1] 0.0191494 9.4489001
    ## 
    ## $value
    ## [1] 32.51089
    ## 
    ## $counts
    ## function gradient 
    ##       73       NA 
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## NULL

``` r
circular::mle.vonmises(data_rndm, mu=NULL, kappa=NULL, bias=FALSE, control.circular=list())
```

    ## 
    ## Call:
    ## circular::mle.vonmises(x = data_rndm, mu = NULL, kappa = NULL,     bias = FALSE, control.circular = list())
    ## 
    ## mu: 0.01917  ( 0.03347 )
    ## 
    ## kappa: 9.443  ( 1.295 )

### Stima di una Von Mises con una covariata

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
  kappa.log <- exp(par[3])
  
  l <- array(dim = length(x))
  
  for(i in seq_along(x)){
    l[i] = kappa.log * cos(y[i] - 2 * atan(beta.0 + beta.1*x[i])) - log((besselI(kappa.log, nu = 0)))
  }
  
  log.lik <- sum(l)
  return(-log.lik)
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

det(result$hessian) != 0 # condizione di invertibilità
```

    ## [1] TRUE

``` r
cov.matrix <- solve(result$hessian)

se <- sqrt(diag(cov.matrix))
se
```

    ##     beta.0     beta.1      kappa 
    ## 0.04995093 0.13589698 0.13843907

### Simulazione di $K$ repliche per lo stesso set di parametri veri

Definisco una funzione che dati

- `n` = la numerosità del campione casuale

- `K` = il numero di repliche

- `beta.0` e `beta.1` = il valore reale dei parametri

- `kappa` = il valore reale del parametro di concentrazione della Von
  Mises

Genera $K$ set di random $x$ e random $y$ e stima i parametri
ottimizzando la funzione di verosimiglianza per ogni set. La funzione
restituisce una tabella in cui le colonne rappresentano i valori stimati
dei parametri e le righe sono le $K$ simulazioni

``` r
vonmises_mle_sim <- function(n, K, beta.0, beta.1, kappa){
  
  # funzione per generare i set di dati 
  gen.data <- function(n, beta.0, beta.1, kappa){
    x <- rnorm(n)
    y <- sapply(x, function(xi) circular::rvonmises(n= 1, mu=(beta.0 + beta.1 * xi), kappa = kappa))
    data <- tibble(x=x, y=y)
  }
  
  # definizione della funzione `log.lik`
  log.lik.vm <- function(par, x, y){
    
    beta.0 <- par[1]
    beta.1 <- par[2]
    kappa <- exp(par[3])
    
    l <- array(dim = length(x))
    
    for(i in seq_along(x)){
      l[i] = kappa * cos(y[i] - 2 * atan(beta.0 + beta.1*x[i])) - log((besselI(kappa, nu = 0)))
    }
    
    log.lik <- sum(l)
    return(-log.lik)
  }
  
  datasets <- map(1:K, ~ gen.data(n, beta.0, beta.1, kappa))
  
  # stima dei parametri 
  results <- map_dfr(
    .x = datasets,
    .f = ~ {
      mle <- optim(
        par = c(beta.0 = 0, beta.1 = 0, kappa = log(1)),
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
  
  return(list(results = results, datasets = datasets))
}

sim_1 <- vonmises_mle_sim(100, 100, -1, 3.5, 10)
sim_1$results
```

    ## # A tibble: 100 × 4
    ##    beta.0 beta.1 kappa convergence
    ##     <dbl>  <dbl> <dbl>       <int>
    ##  1 -0.874  2.88  1.88            0
    ##  2 -1.07   3.45  1.26            0
    ##  3  2.00  -1.13  0.211           0
    ##  4 -0.820  3.14  1.22            0
    ##  5  4.27   0.684 0.139           0
    ##  6 -0.856  3.05  1.29            0
    ##  7 -0.849  3.37  1.40            0
    ##  8  6.71   4.67  0.411           0
    ##  9 -0.884  3.00  1.15            0
    ## 10  2.40  -2.14  0.328           0
    ## # ℹ 90 more rows

I risultati della simulazione sono decisamente da quelli che ci si
potrebbe aspettare. A titolo di esempio, seleziono le righe in cui il
parametro stimato per `beta.0` esce dall’intervallo $\pm 3 \sigma$ dove
$\sigma$ corrisponde al valore stimato tramite l’informazione di Fisher

``` r
print(
sim_1$results %>% 
  filter(!between(beta.0, -1 - 3*se[1], -1 + 3*se[1])),
n = 50)
```

    ## # A tibble: 68 × 4
    ##         beta.0     beta.1    kappa convergence
    ##          <dbl>      <dbl>    <dbl>       <int>
    ##  1      2.00      -1.13   2.11e- 1           0
    ##  2     -0.820      3.14   1.22e+ 0           0
    ##  3      4.27       0.684  1.39e- 1           0
    ##  4     -0.849      3.37   1.40e+ 0           0
    ##  5      6.71       4.67   4.11e- 1           0
    ##  6      2.40      -2.14   3.28e- 1           0
    ##  7      3.61      -1.17   2.07e- 1           0
    ##  8      5.26       3.54   3.53e- 1           0
    ##  9      5.15       3.46   3.78e- 1           0
    ## 10      0.749     -0.405  3.61e- 1           0
    ## 11     20.4       -6.80   7.94e-28           0
    ## 12      3.12       0.811  1.58e- 1           0
    ## 13      7.43       5.04   3.26e- 1           0
    ## 14     -0.828      3.05   9.17e- 1           0
    ## 15     -0.760      2.76   1.18e+ 0           0
    ## 16      5.43       3.43   4.94e- 1           0
    ## 17     -0.705      2.64   1.28e+ 0           0
    ## 18     42.5       12.4    1.74e-34           0
    ## 19      0.231      0.305  2.80e- 1           0
    ## 20      5.44       3.66   3.37e- 1           0
    ## 21      4.43       0.934  1.40e- 1           0
    ## 22      6.06       4.18   5.37e- 1           0
    ## 23     -0.848      3.04   1.33e+ 0           0
    ## 24     -0.837      2.99   1.26e+ 0           0
    ## 25     23          8.03   1.81e-37           0
    ## 26     -0.192     -0.114  2.23e- 1           0
    ## 27     42.5       12.4    1.74e-34           0
    ## 28     -0.838      2.93   1.10e+ 0           0
    ## 29     -0.790      2.88   1.11e+ 0           0
    ## 30 367635.     78007.     3.00e- 1           0
    ## 31      1.50      -0.787  2.92e- 1           0
    ## 32      0.174      0.102  2.88e- 1           0
    ## 33      0.337      0.159  1.77e- 1           0
    ## 34     -0.816      3.26   1.01e+ 0           0
    ## 35     -0.831      3.30   1.17e+ 0           0
    ## 36     12.4       22.9    1.74e-34           0
    ## 37     -0.0701    -0.222  2.31e- 1           0
    ## 38     -0.759      3.19   1.35e+ 0           0
    ## 39     -0.798      2.77   1.35e+ 0           0
    ## 40     39.4       10.7    6.57e-32           0
    ## 41     -0.827      2.90   1.20e+ 0           0
    ## 42      5.63       3.59   2.93e- 1           0
    ## 43     -0.842      3.09   1.33e+ 0           0
    ## 44     -0.839      3.11   1.20e+ 0           0
    ## 45     -0.828      3.02   1.28e+ 0           0
    ## 46      4.21      -0.0148 2.01e- 3          10
    ## 47     -0.841      2.71   1.45e+ 0           0
    ## 48     11.7       19.6    2.16e-40           0
    ## 49     42.5       12.4    1.74e-34           0
    ## 50      2.37      -0.949  1.88e- 1           0
    ## # ℹ 18 more rows
