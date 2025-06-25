#' ---
#' title: "Analisi del primo elefante (id: GON09)"
#' author: "Francesco Invernizzi"
#' date: "`r format(Sys.Date())`"
#' output:
#'  html_document:
#'    code_folding: hide
#' ---
#+ include = F
rm(list = ls())
file.remove(".Rhistory")
setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
source("code/00_master_code.R")
knitr::opts_chunk$set(fig.align = "center", echo = T, warning = F, message = F)
#' Usando il pacchetto `amt` e seguendo i passaggi svolti nel [paper di riferimento del pacchetto](https://doi.org/10.1002/ece3.4823) svolgo le analisi sul primo elefante del `df`. 
#' 
#' ### Preparazione dei dati
#' 
#' Le principali funzioni per la preparazione dei dati sono le seguenti: 
#' 
#' - `amt::make_track()` crea un "track" cioè una tibble contenente le coordinate x e y, il tempo (t) e un sistema di riferimento di coordinate (CRS). Questa tibble verrà utilizzata dal pacchetto per le operazioni successive. Nel nostro caso il CRS è identificato con EPSG:21036, Arc 1960, UTM zone 36S. 
load(paste0(dir_data, "/data_elephants_clean_envcovs.RData"))
GON09 <- data_elephants %>% #ordino i dati per id e data/ora
  arrange(id, timepoint) %>% 
  filter(id == "GON09")

crs <- capture.output(st_crs(GON09))
cat(crs[1:2], sep = "\n")
#' 
#' - `amt::summarize_sampling_rate()` fornisce la distribuzione della frequenza di campionamento e consente di vedere se le registrazione della posizione sono state costanti e con che intervallo. 
GON09_trk <- amt::make_track(GON09,
                             .x = "XUTM",
                             .y = "YUTM", 
                             .t = "timepoint",
                             crs = 21036) 

amt::summarize_sampling_rate(GON09_trk)
#'
#' - `amt::track_resample()` unito a `filter_min_n_burst(3)` permettono di conservare le osservazioni solo nel caso in cui in cui le registrazioni della posizione avvengono ad un intervallo costante (nel nostro caso 4h $pm$ 15 min) per almeno tre registrazioni consecutive 
#' 
#' -  `steps_by_burst()` passa da un `df` in cui ogni riga corrisponde ad una posizione ad uno in cui ogni riga corrisponde ad un passo. Per ogni passo, sono incluse le variabili coordinate di inizio (x1_, y1_ ), cordinate finali (x2_, y2_ ), tempo di inizio e di fine (t1_, t2_ ), step length (sl_; in CRS units), turning angles (ta_; in degrees), la time difference (dt_ ) e la burst (sequenza di passi con intervallo di registrazione di 4h) (burst_ ) a cui appartiene il passo
GON09_trk <-  track_resample(GON09_trk,
                             rate = hours(4), 
                             tolerance = minutes(20)) %>% 
  filter_min_n_burst(3) %>% 
  steps_by_burst() 

knitr::kable(head(GON09_trk, 15), format = "markdown", digits = 2)
#'
#' - `amt::extract_covariates` consente di estrarre delle covariate ambientali da un raster. Dato che il dataset originale associa già alla posizione la distanza dalla riva (distriv), l'indive NDVI (ndvi), l'elevazione (elev) e la stagione (seas) ad ogni posizione, uso `left_join` per associare le covariate alla posizione di arrivo. Come specificato dal paer questa è una prima scelta:
#' 
#' > Depending on the target process under investigation (habitat selection or movement), covariates might be extracted at the end of the step (habitat selection process) or at the start of the step (movement process). If covariates are extracted at the end of the step, they are typically included in the model as main effects, to answer questions of the type: How do covariates influence where the animal moves? In contrary, if covariates are extracted at the beginning of the step, they are typically included in the model as an interaction with movement characteristics (step length, log of the step length, or the cosine of the turn angle), to test hypotheses of the type: Do animals move faster/more directed, if they start in a given habitat? Finally, covariate values at the start and the end of a step can also be included in the model as an interaction with each other, to test hypotheses of the type: Are animals more likely to stay in a given habitat, if they are already in that habitat?
GON09_fnl <- GON09_trk %>% 
  left_join(
    dplyr::select(GON09, XUTM, YUTM, distriv, ndvi, elev, seas),
    join_by("x2_" == "XUTM", "y2_" == "YUTM")
  ) 
#' ### Analisi esplorativa dei dati 
#' Gli istogrammi rappresentano le distribuzioni di step length, turn angle e direction nelle tre stagioni climatiche 
#' 
#' >Three climatic seasons can be recognised inthe Gonarezhou landscape, the hot wet (HW) season (November to March, when 90% of annual rain falls); the cool dry (CD) season (April to August); and the hot dry (HD) season (September to October) (Gandiwa 2014; Republic of Zimbabwe 2016).
ggplot(data = GON09_fnl, aes(x = sl_, y = ..density..)) +
  geom_histogram(color = "black", fill = NA) +
  geom_density() +
  facet_wrap(~ seas) +
  theme_test() +
  xlab("Step length [km]") +
  ylab("Density")

ggplot(data = GON09_fnl, aes(x = ta_, y = ..density..)) +
  geom_histogram(color = "black", fill = NA) +
  geom_density() +
  facet_wrap(~ seas) +
  theme_test() +
  xlab("Turn angle") +
  ylab("Density")

ggplot(data = GON09_fnl, aes(x = direction_p, y = ..density..)) +
  geom_histogram(color = "black", fill = NA) +
  geom_density() +
  facet_wrap(~ seas) +
  theme_test() +
  xlab("Direction") +
  ylab("Density")

leaflet(data = GON09) %>%
  addProviderTiles(providers$Esri.NatGeoWorldMap) %>%  
  addPolylines(lng = ~long, lat = ~lat, color = "blue", weight = 0.5) 
#TODO 
#1. mappa
#2. scatter plots 
# vorrei usare il mio df + tmap mode view

#' ### Fit
#' Nel paper, come riporto qua sotto, viene utilizzato il metodo iSSFs. La variabile dipendente è la step length mentre il turn angle viene usato per stimare le step length random di controllo.
#' 
#' > To fit SSFs, the observed covariates associated with observed steps are compared to covariates associated with random (or control) steps. Random steps can be generated by either (a) sampling from the observed turn step‐length and turn angle distribution (resulting in a traditional SSF), or (b) by fitting a parametric distribution to the observed step lengths (either a negative‐exponential, a halfnormal, a log‐normal, or a gamma; see Avgar et al., 2016, Appendix 2) and turn angles (a von Mises; Duchesne et al., 2015). As mentioned above, an iSSF is arguably less biased and also provides the user with a mechanistic movement model that can be used to simulate space use, and hence utilization distributions (Avgar et al., 2016; Signer et al., 2017). Currently, amt only implements the iSSFs with gamma and von Mises distributions. Thus, we proceed by fitting a gamma distribution to the step lengths and a von Mises distribution to the turn angles using maximum likelihood (Agostinelli & Lund, 2017; Delignette‐Muller & Dutang, 2015), and use these distributions to generate and pair nine random steps with each observed step. The number of random steps effects the estimation error; the more the steps, the lower the error, but the higher the computational burden (Avgar et al., 2016). We then extract the covariates at the end point of each step (observed and random) using the function amt::extract_ covariates, and fit a conditional logistic regression model to the resulting data including movement‐related covariates with the function amt::fit_issf (a wrapper to survival::clogit; Therneau & Grambsch, 2000; see Supporting information Data S2 for details of the model). The target variable case_ is one for observed steps and zero for random (or control) steps. Each step is paired with several (here 9) control steps that form together a stratum (indicated by strat(step _id_) in the model formula). The function amt::random _steps automatically creates a new column, step _id_ , that identifies different strata
#'
#' In questo caso fitto un modello in cui la variabile dipendente è la step length e sono incluse tutte le covariate ambientali a disposizione.
mod1 <- lm(log(sl_) ~ distriv + ndvi + elev + seas, data = GON09_fnl)
summary(mod1)