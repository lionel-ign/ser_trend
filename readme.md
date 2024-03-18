# Companion repository to the manuscript: "Turning point in forest productivity"

This companion repository allow to reproduce the results shown in the manuscript (under review) and available has a preprint in [add link to arxiv].

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