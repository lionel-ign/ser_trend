# Companion repository to the manuscript: "Turning point in the productivity of western European forests associated with a climate change footprint"

[![DOI](https://zenodo.org/badge/773863581.svg)](https://zenodo.org/doi/10.5281/zenodo.10836876)

This companion repository allow to reproduce the results shown in the manuscript (soon to be published in STOTEN) and available has a [preprint](https://www.biorxiv.org/content/10.1101/2024.04.12.589202v1.abstract).

## Structure

This repository is divided into the following subfolders:

- data: the necessary data to run the analysis and reproduce the figures. Note that the raw data are extracted from the french inventory computation engine via an internal R package (inventR) that is currently not distributed outside of IGN. The preprocessed data are available in this subfolder
- script: the following three R scripts
    - 01_create_data.R: create the dataset, as noted above this script is for reference only as the package inventR that derives the estimation is not available for external use. The meteorological data are available from Christian Piedallu (christian.piedallu@agroparistech.fr) upon reasonable request
    - 02_fitting_script.R: from the pre-processed data run the different models using the package brms
    - 03_plotting_script.R: from the fitted models and the preprocessed data, make the figures shown in the manuscript

## Contact

Inquiries should be sent to lionel.hertzog@ign.fr


## Licence

The code developed in this repository is licensed under [CC BY-NC-SE 4.0 DEED](https://creativecommons.org/licenses/by-nc-sa/4.0/). The data originating from the National Forest Inventory are under the Open Licence Etalab Version 2.0.
