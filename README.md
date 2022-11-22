#This branch specific requirements 
The changes in this branch are made for decentralized mancova(this tool) to be compatible with COINSTAC vaults. 
This repo usually requires two csv files, one with fixed covariates and the other specifying those covariate names and types. 
COINSTAC Vault currently allows only one csv file for communication, so both these covariates.csv and covariate_keys.csv 
files are merged into one as follows: 
```
         filename   -  This column lists the GICA output files that are required for MANCOVA compuation
     niftifilename  -  This column lists the nifti file names of the subjects; if GICA output files are provided as inputs then
                       these nifti files need not be actually present in the nifti file, however this column still needs to be 
                       populated for computation. 
covariatename:type  -  The column name contains covariate name and its type separated by a ':' .
                       Example age:continuous or gender:categorical 
                       
An example of this format can be found in test/input/local0/simulatorRun/covariates.csv file.
```

```
For groupica file, please download the GIFT standalone from: https://drive.google.com/file/d/1U-IoiOJQYPikj7uwxS9n4wUdwb66LOum/view 
Extract GroupICATv4.0b_standalone_Aug_05_2021.zip and copy these file into groupicatv4.0b/GroupICATv4.0b_standalone/ directory.  
```

# Coinstac Group ICA/dMANCOVA

This repository compiles submodules utilized for Group ICA and MANCOVA

The dMANCOVA pipeline consists of two distinct parts - Group or Spatially Constrained ICA first, followed by dMANCOVA.

## Running in the Simulator

This pipeline has been tested with version 5.1.0 of the COINSTAC simulator


Install the simulator:

```
npm i -g coinstac-simulator
```

Download this repository

```
git clone https://github.com/trendscenter/coinstac-mancova.git
```

Initialize submodules

```
git submodule update --init --recursive
```

Copy the mask and template into the local input folders, using the bash script

```
bash copy_data.sh
```

*or* the following commands

```
mkdir -p test/input/remote/simulatorRun ;
cp local_data/* test/input/remote/simulatorRun
cp local_data/* test/input/local0/simulatorRun
cp local_data/* test/input/local1/simulatorRun
```

Copy your NIfTI files and corresponding cavariates CSV files into the test site folders.
The `covariates.csv` file should only have data for the NIfTI files in the corresponding site.
Below is an example `covariates.csv` file.

| age | diagnosis | filename | 
| --- | --------- | -------- | 
| 25.0 | 1.0 | vsdwa_000123.nii.gz |
| 42.0 | 0.0 | vsdwa_000456.nii.gz |

After those steps, the `test` folder should look similar to this (probably with different .nii.gz files).

```
test
├── input
│   ├── local0
│   │   └── simulatorRun
│   │       ├── NeuroMark.nii
│   │       ├── covariate_keys.csv
│   │       ├── covariates.csv
│   │       ├── mask.nii
│   │       ├── vsdwa_000123.nii.gz
│   │       ├── vsdwa_000456.nii.gz
│   ├── local1
│   │   └── simulatorRun
│   │       ├── NeuroMark.nii
│   │       ├── covariate_keys.csv
│   │       ├── covariates.csv
│   │       ├── mask.nii
│   │       ├── vsdwa_000789.nii.gz
│   │       ├── vsdwa_000800.nii.gz
│   └── remote
│       └── simulatorRun
│           ├── NeuroMark.nii
│           └── mask.nii
└── inputspec.json
```



Finally, run using the bash script (will require **sudo** unless your user is in the `docker` group)

```
bash run.sh
```

*or*

Run using the following commands

```
docker build  -t dmancova .
coinstac-simulator --loglevel silly
```

## Spatially Constrained ICA

The stages of Spatially Constrained ICA are

 - local spatially constrained ICA (performed with GIFT)

## Group ICA

The stages of Group ICA are:

 - decentralized row means
 - decentralized PCA
 - local ICA (either with Infomax ICA, spatially-constrained ICA, or other)

## ddFNC

The stages of dMANCOVA are:
 - ...
