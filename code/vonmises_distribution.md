Informazioni iniziali
================
Francesco Invernizzi
2025-06-25

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
    ## [1] -0.001591691 13.262618582
    ## 
    ## $value
    ## [1] 14.65037
    ## 
    ## $counts
    ## function gradient 
    ##      101       NA 
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
    ## mu: -0.001635  ( 0.028 )
    ## 
    ## kappa: 13.26  ( 1.837 )

### Stima di una Von Mises con una covariata

1.  Creo due vettori $\mathbf{x}_n$ e $\mathbf{y}_n$ in cui

- $X_i \sim N(0, 1)$

- $Y_i \sim VonMises(\mu = 2atan(\beta_0 +\beta_1 \cdot x_i), k = (10))$

``` r
set.seed(143)
rndm_x <- rnorm(100)
rndm_y <- array()

for(i in seq_along(rndm_x)) {
  rndm_y[i] = circular::rvonmises(n = 1, mu = 2*atan(-1 + 3.5 * rndm_x[i]), kappa = 10)
}
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
    ## 0.05243765 0.15812123 0.13827019

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
  gen.data <- function(n){
    x <- rnorm(n)
    y <- c()
    for (i in 1:n) {
      y[i] <- circular::rvonmises(n= 1, mu=(beta.0 + beta.1 * x[i]), kappa = kappa)
    }
    
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
  
  # stima dei parametri 
  results <- map_dfr(
    .x = map(1:K, ~ gen.data(n)),
    .f = ~ {
      mle <- optim(
        par = c(beta.0 = 0, beta.1 = 0, log.kappa = log(1)),
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
  
  return(results)
}

sim_1 <- vonmises_mle_sim(100, 100, -1, 3.5, 10)
sim_1
```

    ## # A tibble: 100 × 4
    ##     beta.0  beta.1    kappa convergence
    ##      <dbl>   <dbl>    <dbl>       <int>
    ##  1 -1.33    -2.00  2.44e- 1           0
    ##  2 -0.808    2.80  1.56e+ 0           0
    ##  3 32.5    -12.4   5.44e-36           0
    ##  4 -0.939    3.23  1.28e+ 0           0
    ##  5 10.6     17.2   2.57e-30           0
    ##  6 -0.985    3.17  1.35e+ 0           0
    ##  7 -0.734    2.84  1.19e+ 0           0
    ##  8 -0.0404  -0.156 1.82e- 1           0
    ##  9 -0.769    2.94  1.54e+ 0           0
    ## 10  0.411   -2.11  3.02e- 1           0
    ## # ℹ 90 more rows

Verifico se la media dei parametri corrisponde al valore reale del
parametro e se la sd ai valori stimati tramite l’inversa
dell’Informazione di Fisher

``` r
map(.x = sim_1,
    .f = ~ mean(.x))
```

    ## $beta.0
    ## [1] 174.5847
    ## 
    ## $beta.1
    ## [1] 54.62559
    ## 
    ## $kappa
    ## [1] 0.7351773
    ## 
    ## $convergence
    ## [1] 0

``` r
map(.x = sim_1, 
    .f = ~ sd(.x))
```

    ## $beta.0
    ## [1] 1722.113
    ## 
    ## $beta.1
    ## [1] 511.4374
    ## 
    ## $kappa
    ## [1] 0.5481125
    ## 
    ## $convergence
    ## [1] 0
