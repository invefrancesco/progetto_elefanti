Movimento elefanti
================
Francesco Invernizzi
2025-04-17

Prima di iniziare è necessario caricare i pacchetti e impostare le
directories

``` r
setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
source("code/00_master_code.R")
```

    ## Date: 17 April 2025
    ## By: Francesco Invernizzi
    ## Description: code - master code to upload packages and set directories
    ## Version of R used: R version 4.4.3 (2025-02-28)

    ## Loading required package: pacman

## Introduzione

In questo file l’obiettivo è verificare come vengono calcolati turn
angle e step length dal pacchetto `samt`

- Nella prima parte utilizzo le prime dieci righe del dataset per
  calcolare sia analiticamente che con il pacchetto entrambe le misure

- Nella seconda parte invece analizzo i dati osservati per il primo
  elefante a disposizione (GON09)

## Operazioni basilari sui movimenti degli animali

The basic building blocks of amt are tracks. Tracks are tibbles with at
least two columns that contain the coordinates: x\_ and y\_. A track
behaves exactly like a tibble (the only difference being that we added
an other S3 class). Below is an example of creating a track with some
dummy locations.

Importo il df

``` r
load(paste0(dir_data, "/data_elephants_clean_envcovs.RData"))
data_elephants <- data_elephants %>% #ordino i dati per id e data/ora
  arrange(id, timepoint)
str(data_elephants, width = 80, strict.width = "no", nchar.max = 80, give.attr = FALSE)
```

    ## sf [67,002 × 20] (S3: sf/tbl_df/tbl/data.frame)
    ##  $ id       : chr [1:67002] "GON09" "GON09" "GON09" "GON09" ...
    ##  $ long     : num [1:67002] 32.2 32.2 32.2 32.2 32.2 ...
    ##  $ lat      : num [1:67002] -21.2 -21.2 -21.2 -21.1 -21.1 ...
    ##  $ tag      : chr [1:67002] "African Elephant: 1099301" "African Elephant: 1099301" "African Elephant: 1099301" "African Elephant: 1099301" ...
    ##  $ sex      : chr [1:67002] "Female" "Female" "Female" "Female" ...
    ##  $ date     : Date[1:67002], format: "2016-02-01" "2016-02-01" ...
    ##  $ time     : 'hms' num [1:67002] 05:57:45 09:57:49 17:57:52 01:58:49 ...
    ##  $ hr       : num [1:67002] 5 9 17 1 5 9 13 19 13 17 ...
    ##  $ month    : num [1:67002] 2 2 2 2 2 2 2 2 2 2 ...
    ##  $ yr       : num [1:67002] 2016 2016 2016 2016 2016 ...
    ##  $ timepoint: POSIXct[1:67002], format: "2016-02-01 05:57:45" "2016-02-01 09:57:49" ...
    ##  $ id2      : chr [1:67002] "ele.1099301" "ele.1099301" "ele.1099301" "ele.1099301" ...
    ##  $ year2    : chr [1:67002] "2016" "2016" "2016" "2016" ...
    ##  $ seas     : Factor w/ 3 levels "HW","CD","HD": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ geometry :sfc_POINT of length 67002; first list element:  'XY' num [1:2] 413 -2341
    ##  $ XUTM     : num [1:67002] 413 414 414 414 415 ...
    ##  $ YUTM     : num [1:67002] -2341 -2341 -2340 -2338 -2338 ...
    ##  $ distriv  : num [1:67002] 1.173 0.99 0.162 1.551 1.214 ...
    ##  $ ndvi     : num [1:67002] 3634 3364 2912 2822 2564 ...
    ##  $ elev     : num [1:67002] 499 502 504 516 520 ...
