# GEn Analysis
SBS GEn-II analysis framework. This was initially copied over from Provakar Datta's framework for GMn but has since diverged significantly.

# Configuration files

The configuration files are located in the **config** directory. These tells the analysis scripts all the relevant information about a specific run, i.e. beam energy, kinematic angles, raw data location.

# Database files

The database files are located in **DB** directory. These files hold information about the following information.

**F1/2/u/d.csv**: These files all hold form factor data for the RCQM theoretical model

**World_Data.dat**: Holds the world data measurements for GEN.

**Beam_pol.csv**: Holds the beam polarization results from the Moller measurements.

**He3_pol.csv**: Holds the He3 polarization results from the NMR analysis.

**Field_Meas.csv**: Holds the He3 polarization direction.

**Helicity_quality.csv**: Lists the runs where the helicity readback was inaccurate and cannot be used.

**Moller_quality.csv**: Lists the runs where the Moller measurements are unknown and therefore the beam polarization is unknown.

**corrections**: This directory holds files that list all the GEN analysis results with the file name based on the cuts used in the analysis.

# Framework

Several classes are defined to help organize the analysis process. See the **include** and **src** files for all the details. The **Analysis.h** files holds most of the important functions which are used for the form factor extraction. 

# Analyze Raw Data

The main analysis scripts are located in the directory **scripts/replay**. The main script for production data is:
```shell
QuasiElastic_ana.C
```
The main script for simulation data is:
```shell
QuasiElastic_sim_ana.C
```

These files will output a condensed analyzed root file that is much less data than the raw root file. This scripts can take hours to run if you are analyzing hundreds of millions of events.

# Asymmetry Analysis

The asymmetry analysis scripts are located in the directory **scripts/analysis**. These scripts are for GEN physics analysis and focus on making elastic cuts and obtaining asymmetry yields. There are many scripts and all of them have header comments explaining what they do. The main script calculating all the asymmetries is
```shell
Full_GEN_analysis.C
```

This will run all of the asymmetry corrections, however the data must all be checked to ensure that the settings were all appropriate. The list of the most important scripts are as follows:

**accidental_contamination.C**: Calculates the accidental correction.

**Nitrogen_contamination.C**: Calculates the nitrogen correciton.

**pion_contamination.C**: Calculates the pion correction.

**Inelastic_contamination.C**: Calculates the inelastic correction.

**GEN_Extraction.C**: Calculates the proton correction and also finalizes all corrections and writes the results in a text file.