#' ---
#' title: Informazioni iniziali
#' author: Francesco Invernizzi
#' date: "`r format(Sys.Date())`"
#' output: github_document
#' ---
#+ setup, include=FALSE
setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
source("code/00_master_code.R")
options(warn = -1)
#' 
#' Il documento contiene informazioni di base sul pacchetto  `amt` usato nell'analisi
#' 
#' # Misure usate per il tracciamento degli animali e calcolate dal pacchetto `amt`
#' 
#' Creo un `data.frame` con quattro punti e calcolo sia analiticamente che con il pacchetto `amt` le tre misure di base: step length, turn angle e direction.
#' 
#' - la "step length" è misurata come [distanza euclidea](https://en.wikipedia.org/wiki/Euclidean_distance) tra due rilevazioni successive, quindi conserva l'unità di misura dei punti.
#' 
#' - direction e turn angle sono misurati rispettivamente come [`atan2`](https://en.wikipedia.org/wiki/Atan2) e differenza tra `atan2` consecutivi.
df <- data.frame(
  x = c(0, 1, 1, 0),
  y = c(0, 0, 1, 1),
  t = as.POSIXct(c("2025-01-01 00:00:00", "2025-01-01 01:00:00",
                 "2025-01-01 02:00:00", "2025-01-01 03:00:00"))
)
# plot:
ggplot(df, aes(x, y)) +
  geom_path() +
  geom_point() +
  coord_equal()
#' ### Calcoli analitici: 
#'  1. **Step length** (distanza euclidea)
sqrt((1 - 0)^2 + (0 - 0)^2) 
sqrt((1 - 1)^2 + (1 - 0)^2) 
sqrt((0 - 1)^2 + (1 - 1)^2) 
#'  2. **Direction** (`atan2`)
dir1 = atan2(0, 1)           
dir2 = atan2(1, 0) 
dir3 =atan2(0, -1) 
#'
dir1
dir2
dir3
#'  3. **Turn angle** (differenza tra `atan2` successivi)
ta_2 = dir2 - dir1 
ta_3 = dir3 - dir2
#'
ta_2
ta_3
#' ### Usando il pacchetto `amt`
trk <- make_track(df, x, y, t) # sistema proiettato (UTM zona 32N)

steps <- trk %>% 
  amt::track_resample(rate = hours(1)) %>% 
  steps_by_burst()

steps %>% dplyr::select(x1_, y1_, x2_, y2_, sl_, direction_p, ta_)