# Physiological_stresses_drive_range_shifts

R code and data associated with the paper "_Physiological Stresses Drive Range Shifts in Ectotherms_" by Anderson & Chapple, accepted in Journal of Biogeography.

## Abstract

Climate change is expected to impose severe physiological stress on ectotherms, potentially reshaping their geographic distributions. Many species are expected to move to higher latitudes to avoid detrimental effects of climate change on their homeostatic balance, yet little is known about the ecophysiological mechanisms driving their redistribution and the extent their geographic range will alter. Here, we assess the impacts of future climate change scenarios on Lampropholis skinks, a group of small lizards from eastern Australia, by integrating physiological data with biophysical and species distribution models. We compared the current physiological stress they face currently to those predicted in future climatic scenarios (+2 °C and +4 °C), and we evaluated whether their geographic range will shrink or expand with climate change using mechanistic species distribution modelling. We found that all species will experience increased dehydration, higher metabolic costs, and prolonged exposure to critical thermal limits. These physiological constraints will reduce activity time and drive range contractions, particularly under +4 °C scenarios, where suitable habitats could shrink by half. Most species are predicted to shift to higher latitudes in search of more suitable habitats. Even widespread generalist species, often considered more resilient, are projected to face significant physiological stress, challenging the assumption that broad climatic tolerance ensures future persistence. Our findings underscore the vulnerability of both specialists and generalists ectotherms to climate change, with implications for biodiversity conservation.


## Contents #

### Scripts 

_ensembleSDM_lampropholis.R_: This script runs Ensemble Species Distribution Models (SDMs) for Lampropholis skinks using biophysical modeling-derived data. It also includes code for estimating range size, calculating centroid shifts, and generating Figures 2, 3, and 4.

_analyses_lampro.R_: This script performs phylogenetic data analysis to assess changes in physiological stress (dehydration, thermal stress, energetic demands, and activity hours), range size, and range shifts in Lampropholis skinks across different climatic scenarios.

_skinks_model.R_: This script performs biophysical modeling for Lampropholis skinks across Australia using NicheMapR (Kearney & Porter, 2020), parameterized with data from Anderson et al. (2023).

### Data 

_Lampro_data.zip_: Data generated through biophysical modeling for Lampropholis skinks in Australia. Each dataset contains physiological, behavioral, and stress information during activity and total periods across 11,005 coordinates in Australia.

_aus_coord_0.25_.csv: Coordinates for biophysical modeling, generated in QGIS at 0.25-degree resolution across Australia.

Physiological, behavioural and morphological data to parametrise the biophysical models were extracted from Table S2 in Anderson et al. (2023).

###  References 

Anderson, R. O., Tingley, R., Hoskin, C. J., White, C. R., & Chapple, D. G. (2023). Linking physiology and climate to infer species distributions in Australian skinks. Journal of Animal Ecology, 92(10), 2094-2108.

Kearney, M. R., & Porter, W. P. (2020). NicheMapR–an R package for biophysical modelling: the ectotherm and dynamic energy budget models. Ecography, 43(1), 85-96.

