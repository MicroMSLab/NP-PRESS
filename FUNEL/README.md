# FUNEL

## Init ##

To run all scripts in this part, you need to install `R`. We recommend using conda to manage the version of `R` and its packages. For your convenience, we provide the `environment.yml` file, which can be used to create the required environment. This file specifies the dependencies needed for running the R code in this part of the repository.

```bash
conda env create -f environment.yml
```

Even if renv is installed automatically the first time you launch R, you still need to install `renv` manually within `R`.

```R
install.packages("renv")
```

You need to install all packages in `renv.lock` file. 

```R
renv::restore()
``` 

You also need to provide your own SQLite database file to search against. This File should contain a Table whose name is the same as the `dbtype` in `FUNEL_sample.json`. The table should contain `ID` and `pubchem_mass`.

## Config file ##

`FUNEL_control.json` and `FUNEL_sample.json` are the config files for `FUNEL_control.r` and `FUNEL_sample.r` file respectively. 

`MergeControlPeaklistFile` is the path of the merged control file. 

`MergeControlPeaklist2File` is also the path of the merged control file, but it is used for the second round peak detection.

`SampleFileExtension` is the file type of the sample file, file type can be `mzXML` or `mzML`. Recommend to use mzXML.

`ControlFileExtension` is the file type of the control file, file type can be `mzXML` or `mzML`. 

`DatabaseFile` is the path of your database file.

`dbtype` is the database Tabel to search. 

`rulesfile` is the rules file for adduct calculation, used by CAMERA. 

`IsotopeMatch` is a logical value, if it is TRUE, a signal is only filtered if it is identified as the same isotope peak in both control and experimental samples (e.g. both as M+1). Default is `FALSE`. 

`SampleThreshold` is the Intensity threshold for signals in experimental samples. Default is 10000. 

`ControlThreshold` is the Intensity threshold for signals in control samples. Default is 1000. 

`method` is the method to detect peaks, it can be `centWave` or `massifquant`. Default is `centWave`. 

`ppm` is the ppm value used for matching signals between samples. It should be a integer value between 0 and 1000, depending on the MS instrument conditions. Default is `20`.

`deltaMZ` is the Mass difference tolerance for matching signals between samples. If set as 0, then only ppm parameter will be used. It should be a float value between 0 and 0.2. This parameter is complementary to ppm. Default is `0.01`.

`lowPeakwidth` is the Low boundary for EIC peak width (in seconds.). Normally set as 90% of the average peak width.It should be an integer value between 1 and 100, default is `9`. 

`highPeakwidth` is the Up boundary for EIC peak width (in seconds.). Normally set as 150% of the average peak width.It should be an integer value between 1 and 100 and large than `lowPeakwidth`, Default is 20. 

`rtCheck` is the Retention time difference tolerance for maching signals within a sample. It should be an integer value between 1 and 300. Default is `30` (seconds.). 

`ControlNoise` is the Noise level for controls, signals with intensities lower than this value will be considered as noise. Default is `100`. 

`SampleNoise` is the Noise level for experimental samples, signals with intensities lower than this value will be considered as noise. Default is `1000`. 

`calAdductS` is a logical value, if it is `TRUE`, adducts (e.g. [M+Na]+, [2M+H]+) will be annotated for experimental samples. Default is `TRUE`. 

`snthresh` is the Signal to noise ratio threshold for a signal to be considered valid. It should be an integer value. Default is 2. 

`perfwhm` is the Percentage of the width of the Full Width at Half Maximum, a parameter for centWave method. It should be a float  value between 0.5 and 1.5. Default is 1. 

`consecMissedLimit` is the Parameter for massifquant method, ignored if centWave method is selected. It can be `1`, `2`or`3`. Default is 1. 

`integrate` is the Parameter for massifquant method, ignored if centWave method is selected. It can be `1` or `2`. Default is 2.

`looseCtrRemove` is a logical value, if it is TRUE, control signals will be filtered from experimental samples with more stringent criteria. Default is `TRUE`.

`rtThrFilter` is the Retention time difference tolerance for maching signals between samples. It should be an integer value between 1 and 600. Default is `60` (seconds.).

`foldchangeFilter` is the parameter used to determine how to filter signals in experimental samples. When comparing the intensity of a common signal in experimental samples and reference/control samples, if the fold change of intensity is higher than a threshold, the signal will not be filtered. When set as 0, all the signals from experimental samples that can be found in reference/control samples will be filtered. This parameter should be an integer value between 1 and 100. Default is `0`.

`cutoff` is the parameter used to determine what signals can be used for adduct/fragment calculation. Only signals with intensities higher than this threshold will be used for calulation. Default is 5000.

`mzdiffCalc` is a logical value, if it is TRUE, context-based adduct/fragment calculation will be performed. Default is `TRUE`.

## Script ##

Place the reference/control sample spectrum files (.mzXML) in the path set in `FUNEL_control.json` by `ControlFolder`, and the experimental sample spectrum files (.mzXML) in the path set in `FUNEL_sample.json` by `SampleFolder`. 

Then, run the `run.sh` file to perform FUNEL.

**Notice** : This Scrip  can run on Linux. To run on Windows, you can convert `run.sh` to a PowerShell script.

## Known Issues ##

1. It is known that installing the `Amber/AmberTools` by compiling the source code results in the failure of the `ncdf4` installation in `renv::restore()` when using the same conda (whether it's the same or different environment). The reason for this issue is still unclear. Deleting the section related to Amber in the .bashrc file and then restarting the shell can resolve this issue.