#' ---
#' title: "Rmarkdown"
#' author: "Kirill"
#' date: "2020-11-16"
#' output:
#'   html_document:
#'     code_folding: hide
#' ---


#+ echo=F
knitr::opts_chunk$set(fig.align = "center", echo = T, warning = F, message = F)

#+ 
rm(list = ls())
setwd("~/pCloud Drive/francesco_invernizzi/Phd - Roma Tre/06_elefanti")
source("code/00_master_code.R")


# Import dataset ----
load(paste0(dir_data, "/data_elephants_clean_envcovs.RData"))
data_elephants <- data_elephants %>% #ordino i dati per id e data/ora
  arrange(id, timepoint) %>% 
  st_drop_geometry() # elimino la sf perché non serve
str(data_elephants, width = 80, strict.width = "no", nchar.max = 80, give.attr = FALSE)

#' We proceed using nearly all the same steps as in the first example, but with a different data structure: data _ frames with list columns (Müller & Wickham, 2018). List columns are best thought of as regular columns of a data _ frame that are R lists and can contain any objects (in our case tracks and fitted models). The purrr::nest command can be used to nest data into a list column (Henry & Wickham, 2017).

#+ 
data_elephants_tmp <- data_elephants %>% 
  nest(data = -id) 

#We can now apply the steps as before for all animals. We first create a track for each animal and transform the coordinate reference system using the function amt::transform _ coords. 
data_elephants_tmp <- data_elephants_tmp %>% 
  mutate(trk = map(data, function(d) {
    amt::make_track(d,
                    .x = "long",
                    .y = "lat", 
                    .t = "timepoint",
                    crs = 4326) 
  }))

data_elephants_tmp %>%
  mutate(sr = map(trk, summarize_sampling_rate)) %>%
  dplyr::select(id, sr) %>% 
  unnest(cols = c(sr))

# Anche in questo caso nella maggior parte dei casi la registrazione avviene ogni quattro ore. A parte gli ultimi due elefanti (in cui avviene ogni ora), primo secondo e terzo quartile si discostano molto poco dalle quattro ore. Per questo come nel caso precedente faccio il resample a quattro ore. 

data_elephants_tmp <- data_elephants_tmp %>% 
  mutate(steps = map(trk, function(x) {x %>% 
      amt::track_resample(rate = hours(4), tolerance = minutes(20)) %>% 
      amt::filter_min_n_burst(min_n = 3) %>% 
      amt::steps_by_burst()
  }))

data_elephants_fnl <- data_elephants_tmp %>%
  unnest(steps) %>%
  dplyr::select(!c(data, trk, burst_)) %>% 
  left_join(
    dplyr::select(data_elephants, long, lat, id, sex, distriv, ndvi, elev, seas),
    join_by("x2_" == "long", "y2_" == "lat", "id" == "id")
  ) 

str(data_elephants_fnl, width = 80, strict.width = "no", nchar.max = 80, give.attr = FALSE)

# A questo punto come prima visualizzo i grafici delle densità e faccio il modello lineare
ggplot(data = data_elephants_fnl, aes(x = sl_, y = ..density..)) +
  geom_histogram(color = "black", fill = NA) +
  geom_density() +
  facet_wrap(~ sex) +
  theme_test() +
  xlab("Step length [m]") +
  ylab("Density")

ggplot(data = data_elephants_fnl, aes(x = ta_, y = ..density..)) +
  geom_histogram(color = "black", fill = NA) +
  geom_density() +
  facet_wrap(~ sex) +
  theme_test() +
  xlab("Turn angle") +
  ylab("Density")

ggplot(data = data_elephants_fnl, aes(x = direction_p, y = ..density..)) +
  geom_histogram(color = "black", fill = NA) +
  geom_density() +
  facet_wrap(~ sex) +
  theme_test() +
  xlab("Direction") +
  ylab("Density")

mod2 <- lm(sl_ ~ distriv + ndvi + elev + seas + sex, data = data_elephants_fnl)
summary(mod2)
