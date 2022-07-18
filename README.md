🔴❗UNDER DEVELOPMENT❗🔴

# Principal Component Analysis Automation (PCAauto)

PCAauto is a script developed to automate the PCA analysis of Molecular Dynamics.

Available methods are:

* PCA analysis
* Calculation of variance, cumulated varianc and explained variance (this last one only for 3 components for now)
* Calculation of PCA through simulation time
* Visualization o the movements o the components on model (needs bugfix)
* Also, there is the possibility to plot graphs (under development)

## Instalation

Download the code and unzip it on the desirable directory. To prepare the environment use the following command:

```
conda env create -f environment.yml
``` 

Be aware to uncomment the sections on the environment.yml file depending on which OS you are using.

## How to use

Activate the environment using:

```
conda activate pca-mdaanalysis
```

You may modify the main.py file depending on which calculation you are interested.

## Observations

This script was developed following the MDA Analysis documentation, specifically the [MDA Analysis examples](http://minium.com.au/UserGuide/stable/examples/analysis/reduced_dimensions/pca.html) page.

## Authorship

* Author: **Brenda Ferrari** ([brendaferrari](https://github.com/brendaferrari))
