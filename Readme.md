# Progetto elefanti

## Paper 1: Animal movement tools (amt): R package for managing tracking data and conducting habitat selection analyses

Nathan et al. (2008) introduced the movement ecology paradigm that conceptually connects different factors shaping the realized movement path of animals (e.g., the internal state of an animal, interaction with intra and conspecifics, and varying environmental conditions)

-   The movement ecology paradigm can serve as a framework for generating new hypotheses about animal movements

-   We had two primary objectives in developing the amt R package, namely to provide:

    1.  a set of functions for exploratory analyses of movement data in R

    2.  functions for inferring habitat selection from movement data in R

> **Da capire questa parte**: Methods that quantify habitat selection by linking environmen‐ tal covariates to location data of animals have been around for a long time. Traditionally, resource‐selection functions (RSF; Boyce & McDonald, 1999; Manly, McDonald, Thomas, McDonald, & Erickson, 2007) were used to study habitat selection of animals. RSFs compare covariates associated with locations where the animal was observed with covariates associated with random locations within the avail‐ ability domain, a spatial domain within which any location is assumed available for the animal to use at any given time. Despite the sensi‐ tivity of the resulting inference to habitat availability (Beyer et al., 2010), no consensus exists as to the most suitable approach to delin‐ eate the spatial domain of availability (Northrup, Hooten, Anderson, & Wittemyer, 2013; Paton & Matthiopoulos, 2016; Prokopenko, Boyce, & Avgar, 2017b). Moreover, the assumption that the availabil‐ ity domain can be considered temporally static might have been jus‐ tifiable for very coarse sampling rates (e.g., daily or weekly positions of the animal), but is challenging for modern GPS data with sampling rates \<1 hr. SSFs (Fortin et al., 2005; Thurfjell, Ciuti, & Boyce, 2014) resolve these issues by pairing each observed location with a set of random locations deemed accessible from the previously observed location. SSFs estimate conditional selection coefficients using a likelihood equivalent of a Cox proportional hazards model (Gail, Lubin, & Rubinstein, 1981).

A typical workflow to analyze animal tracking data can be divided into four main steps (described in detail below):

1.  Data preparation, inspection, and management: Load and inspect gaps in the data, resample tracks if needed, and adjust coordinate reference systems.

2.  Exploratory data analysis and descriptive analyses: Explore patterns in the data graphically, consider multiple movement characteristics (e.g., step‐length distribution, net square displacement, or home‐range size) across several animals and/or time periods.

3.  Modeling: Fit models to answer questions or test hypothesis related to movement and space use of animals.

4.  Simulation: Use fitted models to simulate derived quantities (e.g., space use) and assess model fit.

## Cose che non conosco

-   sl e ta analiticamente

-   modello logit condizionato e distribuzione von Mises

-   come stimare la distribuzione dei dati / verificare se la distribuzione è corretta

-   ragioni per usare il log della distanza come covariata
