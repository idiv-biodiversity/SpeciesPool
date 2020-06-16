# SpeciesPool
R Package to calculate plot-specific species pools using large vegetation databases

For installing the package in a local Rstudio, use   
library(devtools)  
install_github("idiv-biodiversity/SpeciesPool", ref="master")  
  
as Published in 
Bruelheide H, Jiménez-Alfaro B, Jandt U, Sabatini FM (2020) Deriving site-specific species pools from large databases. Ecography, 43, 1–14.
https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.05172


# Abstract

Defining the species pool of a community is crucial for many types of ecological analyses, providing a foundation to metacommunity, null modelling or dark diversity frameworks. It is a challenge to derive the species pool empirically from large and heterogeneous databases. Here, we propose a method to define a site‐specific species pool (SSSP), i.e. the probabilistic set of species that may co‐occur with the species of a target community. Using large databases with geo‐referenced records that comprise full plant community surveys, our approach characterizes each site by its own species pool without requiring a pre‐defined habitat classification. We calculate the probabilities of each species in the database to occur in the target community using Beals’ index of sociological favourability, and then build sample‐based rarefaction curves from neighbouring records with similar species composition to estimate the asymptotic species pool size. A corresponding number of species is then selected among the species having the highest occurrence probability, thus defining both size and composition of the species pool. We tested the robustness of this approach by comparing SSSPs obtained with different spatial extents and dissimilarity thresholds, fitting different models to the rarefaction curves, and comparing the results obtained when using Beals co‐occurrence probabilities or presence/absence data. As an example application, we calculated the SSSPs for all calcareous grassland records in the German Vegetation Reference Database, and show how our method could be used to 1) produce grain‐dependent estimations of species richness across plots, 2) derive scalable maps of species richness and 3) define the full list of species composing the SSSP of each target site. By deriving the species pool exclusively from community characteristics, the SSSP framework presented here provides a robust approach to bridge biodiversity estimations across spatial scales.
