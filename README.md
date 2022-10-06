## Ecological Niche Models Applied to Fossil Data

Repository maintainer: Marlon E. Cobos [E-mail](mailto:manubio13@gmail.com)

License: <a href="https://github.com/marlonecobos/ENM_paleo/blob/main/LICENSE" target="_blank">GPL3</a> 

### Description

This repository contains R scripts and materials to run examples for the course. Example data can be obtained using the scripts included in this repository.

<hr>

### Course mentors by section 

- Introduction: **Alycia Stigall** (University of Tennessee)
- ENM theory: **Erin Saupe** (University of Oxford)
- Data collection: **Cori Myers** (University of New Mexico)
- Algorithms, modeling, comparisons: **Hannah Owens** (Copenhagen University)
- Algorithms, modeling, comparisons: **Marlon Cobos** (University of Kansas)
- Model interpretation, MegaSDM: **Jenny McGuire** (Geogria Tech)

<hr>

### Information in this repository

Information is organized in folders that contain data and commented R scripts, as follows:

- <a href="https://github.com/marlonecobos/ENM_paleo/tree/main/Data" target="_blank">Occurrence data</a>
    - <a href="https://github.com/marlonecobos/ENM_paleo/blob/main/Data/mammut_americanum_all.csv" target="_blank">Mammoth data</a> downloaded from <a href="https://paleobiodb.org" target="_blank">paleobiodb.org</a>
    - <a href="https://github.com/marlonecobos/ENM_paleo/blob/main/Data/mammut_americanum_0.00623mya.csv" target="_blank">Mammoth data (period ~0.00623 mya)</a>
    - <a href="https://github.com/marlonecobos/ENM_paleo/blob/main/Data/mammut_americanum_0.01585mya.csv" target="_blank">Mammoth data (period ~0.01585 mya)</a>
    - <a href="https://github.com/marlonecobos/ENM_paleo/blob/main/Data/mammut_americanum_0.13mya.csv" target="_blank">Mammoth data (period ~0.13 mya)</a>
    - <a href="https://github.com/marlonecobos/ENM_paleo/blob/main/Data/mammut_americanum_0.787mya.csv" target="_blank">Mammoth data (period ~0.787 mya)</a>
    - <a href="https://github.com/marlonecobos/ENM_paleo/blob/main/Data/mammut_americanum_3.205mya.csv" target="_blank">Mammoth data (period ~3.205 mya)</a>
    - <a href="https://github.com/marlonecobos/ENM_paleo/blob/main/Data/mammut_americanum_3.3mya.csv" target="_blank">Mammoth data (period ~3.3 mya)</a>
    
- <a href="https://github.com/marlonecobos/ENM_paleo/blob/main/Scripts/data_preparation.R" target="_blank">Data preparation</a>
    - Download occurrence data
    - Download environmental data
    - Occurrence data cleaning and spatial thinning
    - Defining areas for model calibration
    - Preparing environmental data for modeling
    
- <a href="https://github.com/marlonecobos/ENM_paleo/blob/main/Scripts/model_calibration_projections.R" target="_blank">Model calibration, projections, and uncertainty</a>
    - Preparing data for model calibration
    - Model calibration (model selection) 
    - Model projections
    - Uncertainty analysis (MESS, MOP)

- <a href="https://github.com/marlonecobos/ENM_paleo/blob/main/Scripts/model_comparisons.R" target="_blank">Model comparisons</a>
    - Suitable areas in projections
    - Predicting records with projections 
    - Niche overlap analysis
