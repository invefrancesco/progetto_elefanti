#' ---
#' title: Movimento elefanti
#' author: Francesco Invernizzi
#' date: "`r format(Sys.Date())`"
#' output: github_document
#' ---

#' Prima di iniziare è necessario caricare i pacchetti e impostare le directories
setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
source("code/00_master_code.R")

#' ## Introduzione 
#' In questo file l'obiettivo è verificare come vengono calcolati turn angle e step length dal pacchetto `samt` 
#' 
#'  - Nella prima parte utilizzo le prime dieci righe del dataset per calcolare sia analiticamente che con il pacchetto entrambe le misure 
#'  
#'  - Nella seconda parte invece analizzo i dati osservati per il primo elefante a disposizione (GON09)
#'  
#' ## Operazioni basilari sui movimenti degli animali 
#' 
#' The basic building blocks of amt are tracks. Tracks are tibbles with at least two columns that contain the coordinates: x_ and y_. A track behaves exactly like a tibble (the only difference being that we added an other S3 class). Below is an example of creating a track with some dummy locations.
#' 
#' Importo il df
load(paste0(dir_data, "/data_elephants_clean_envcovs.RData"))
data_elephants <- data_elephants %>% #ordino i dati per id e data/ora
  arrange(id, timepoint) %>% 
  st_drop_geometry()
str(data_elephants, width = 80, strict.width = "no", nchar.max = 80, give.attr = FALSE)

#' ### Distanza euclidea 
#' 
#' Il dataframe contiene latitudine e longitudine espresse sia in gradi che in UTM per la zona di riferimento, cioè la zona 36S (codice EPSG: 21036)
#' 
#' Seleziono solamente le prime 10 righe per svolgere le operazioni basilari e confrontare i risultati con quelli del pacchetto `amt` 
#' 
#' Svolgo la trasformazione basilare del `df` in track `amt` che crea per ogni riga una x e una y di partenza e una x e una y di arrivo. Calcolo la distanza euclidea a mano e la confronto con quella calcolata dal pacchetto
#' 
#' The function `amt::make_track` creates a track (the basic building block of the amt package), given the names of the columns containing x and y coordinates, time (t), and we can set a coordinate referencesystem (CRS) 
#' stps is a regular data_frame with 11 attributes of steps (e.g., start, end, and step length; columns). For each step, the start (x1_, y1_ ) and end (x2_, y2_ ) coordinates, as well as the start and end time (t1_, t2_ ) are given. In addition, the following derived quantities are calculated: step length (sl_; *in CRS units*), turning angles (ta_; in degrees; notice that it cannot be calculated for steps that are not preceded by a valid step), the time difference (dt_ ), and the burst (burst_ ) to which the step belongs.

GON09 <- data_elephants %>% #
  filter(row_number() <= 10) 

GON09_trk <- amt::make_track(GON09,
                             .x = "XUTM",
                             .y = "YUTM", 
                             .t = "timepoint",
                             crs = 21036) 

amt::summarize_sampling_rate(GON09_trk)

GON09_trk <-  track_resample(GON09_trk,
                             rate = hours(4), 
                             tolerance = minutes(20)) %>% 
  filter_min_n_burst(3) %>% 
  steps_by_burst()

