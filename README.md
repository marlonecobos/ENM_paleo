## Ecological Niche Models Applied to Fossil Data

Repository maintainer: Marlon E. Cobos [E-mail](mailto:manubio13@gmail.com)

License: <a href="https://github.com/marlonecobos/ENM_paleo/blob/main/LICENSE" target="_blank">GPL3</a> 


-   [Description](#description)
-   [Course mentors by section](#course-mentors-by-section)
-   [Information in this repository](#information-in-this-repository)
-   [Installation instructions](#installation-instructions)
    -   [Download software](#download-software)
    -   [Detailed instructions](#detailed-instructions)

<hr>

### Description

This repository contains R scripts and materials to run examples for the 2022 Paleontological Society short course "Ecological Niche Models Applied to Fossil Data". Example data can be obtained using the scripts included in this repository.

<hr>

### Course mentors by section 

- Introduction: **Alycia Stigall** (University of Tennessee)
- ENM theory: **Erin Saupe** (University of Oxford)
- Data collection: **Cori Myers** (University of New Mexico)
- Algorithms, modeling, comparisons: **Hannah Owens** (Copenhagen University)
- Algorithms, modeling, comparisons: **Marlon E. Cobos** (University of Kansas)
- Model interpretation, MegaSDM: **Jenny McGuire** (Geogria Tech)

<hr>

### Information in this repository

Information is organized in folders that contain data and commented R scripts, as follows:

- <a href="https://github.com/marlonecobos/ENM_paleo/tree/main/Data" target="_blank">Occurrence data</a>
    - <a href="https://github.com/marlonecobos/ENM_paleo/blob/main/Data/https://github.com/marlonecobos/ENM_paleo/blob/main/Data/mammut_americanum_raw.csv" target="_blank">Mastodon data (raw)</a> downloaded from <a href="https://paleobiodb.org" target="_blank">paleobiodb.org</a>
    - <a href="https://github.com/marlonecobos/ENM_paleo/blob/main/Data/mammut_americanum_all.csv" target="_blank">Complete Mastodon data (initial subset)</a>
    - <a href="https://github.com/marlonecobos/ENM_paleo/blob/main/Data/mammut_americanum_0.00623mya.csv" target="_blank">Mastodon data (period ~0.00623 mya)</a>
    - <a href="https://github.com/marlonecobos/ENM_paleo/blob/main/Data/mammut_americanum_0.01585mya.csv" target="_blank">Mastodon data (period ~0.01585 mya)</a>
    - <a href="https://github.com/marlonecobos/ENM_paleo/blob/main/Data/mammut_americanum_0.13mya.csv" target="_blank">Mastodon data (period ~0.13 mya)</a>
    - <a href="https://github.com/marlonecobos/ENM_paleo/blob/main/Data/mammut_americanum_0.787mya.csv" target="_blank">Mastodon data (period ~0.787 mya)</a>
    - <a href="https://github.com/marlonecobos/ENM_paleo/blob/main/Data/mammut_americanum_3.205mya.csv" target="_blank">Mastodon data (period ~3.205 mya)</a>
    - <a href="https://github.com/marlonecobos/ENM_paleo/blob/main/Data/mammut_americanum_3.3mya.csv" target="_blank">Mastodon data (period ~3.3 mya)</a>
    
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
    
- <a href="https://github.com/marlonecobos/ENM_paleo/blob/main/Scripts/model_evaluation.R" target="_blank">Model evaluation</a>
    - Variable importance (contribution to models)
    - Thresholding results
    - Statistical significance, fitting, and complexity 

- <a href="https://github.com/marlonecobos/ENM_paleo/blob/main/Scripts/model_comparisons.R" target="_blank">Model comparisons</a>
    - Suitable areas in projections
    - Predicting records with projections 
    - Niche overlap analysis

<hr>

### Installation instructions

All the analysis that can be performed using the scripts and data in this repository were done in R using functions from different packages. See the list of programs and packages you will need and how to download them below.

#### Download software

This is the list of programs and packages to be used and where to get them. See detailed instructions to install them as needed. R packages can be installed following instructions in the next section.

- R available at <a href="https://cran.r-project.org/bin/windows/base/">link for Windows</a>, <a href="https://cran.r-project.org/bin/macosx/">link for Mac</a>, and <a href="https://cran.r-project.org/">link for Linux</a>.
- Rstudio (optional) available at <a href="https://www.rstudio.com/products/rstudio/download/#download" target="_blank">link</a>.
- Rtools available at <a href="https://cran.r-project.org/bin/windows/Rtools/" target="_blank">link</a>.
- Xcode available at <a href="https://apps.apple.com/us/app/xcode/id497799835?mt=12" target="_blank">link</a>
- JDK available at <a href="https://www.oracle.com/java/technologies/downloads/#java19" target="_blank">link</a>, also at this <a href="https://aws.amazon.com/corretto/?filtered-posts.sort-by=item.additionalFields.createdDate&filtered-posts.sort-order=desc" target="_blank">link</a>.
- Maxent available at <a href="https://biodiversityinformatics.amnh.org/open_source/maxent/">link</a>.
- R packages from CRAN (remotes, raster, maps, paleobioDB, and spThin).
- R packages from GitHub (kuenm, ellipsenm).

#### Detailed instructions

1. Install R and Rstudio in your computer. 
2. Install the compilations tools needed:
    - In Mac, if this is the first time doing so, try running in the terminal `xcode-select --install`. More instructions <a href="https://www.freecodecamp.org/news/how-to-download-and-install-xcode/">here</a>.
    - In Windows, for RTools 4.0 follow this <a href="https://cran.r-project.org/bin/windows/Rtools/rtools40.html/">instructions</a>.
    - In Linux, install build essential tools for compilation and other pieces of software may be required. Searches of how to install each of them will be required.
3. Install JDK or Amazon corretto.
4. Download Maxent, unzip it and place it in an appropriate fixed location.
5. Install packages from CRAN in R using the function `install.packages("package_name")`.
6. Install packages from GitHub in R as follows:
    - `remotes::install_github("marlonecobos/ellipsenm")`
    - `remotes::install_github("marlonecobos/kuenm")`
