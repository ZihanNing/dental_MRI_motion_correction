# 3D-moco-siemens
This repository contains the reconstruction code to perform data-driven motion correction on Cartesian acquisitions. For optimal performance, a self-navigated trajectory should be deployed (e.g. DISORDER). 
It is built on the methods developed under the [DISORDER framework](https://onlinelibrary.wiley.com/doi/10.1002/mrm.28157). The current version includes the following extensions:
- Dealing with Siemens data structures.
- Extension to non-self-navigated trajectories as well.
- The option to correct for pose-dependent B0 correction, as presented in [Brackenier et al., MRM, 2022](https://onlinelibrary.wiley.com/doi/10.1002/mrm.29255)
- The option to correct for temporally varying low-order B0 variations, as presented in [Brackenier et al., ISMRM, 2024]
- Using Pilot Tone signals to guide the motion correction [Brackenier et al., MRM, 2024](https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.30161)

## Table of Contents
- [Dependencies](#Dependencies)
- [Installation](#Installation)
- [Usage](#Usage)
- [Contributors](#Contributors)
- [License](#License)

## Dependencies
The code runs in MATLAB and does not require a specific version.
This repository depends on 2 other repositories built in-house:
- The code base developed for the [Developing Human Connectome Project](https://www.developingconnectome.org/). 
- The code base developed during Yannick Brackenier's thesis (brackenier-tools). 

## Installation
To install all necessary scripts, this and other repositories need to be cloned for this repository to run. In order to do so, run the following commands in the shell:
```sh
git clone git@github.com:ybrackenier/3D-moco-siemens.git
cd 3D-moco-siemens
git clone git@github.com:ybrackenier/brackenier-tools.git
git clone git@github.com:ybrackenier/dhcp-repo-release07.git
```

## Usage
First, the names of the files that need to be reconstructed must be provided in the script /Studies/studies.m.
Next, the script deployRecon.m needs to be run. The output will automatically be written into the next folders:
- *.mat: These are the raw data files after converting the raw Siemens files (.dat files).
- *.json: These are the JSON files that contain most of the raw data and sequence parameters. This can be useful to skim the sequence and acquisition setup.
  
- An-Ve/: The NIFTI files of the reconstructions. Files with the suffix **_Aq** are the uncorrected files whereas the suffix **_Di** indicates motion-corrected files.
- An-Ve_Sn/: The snapshots of (intermediate) results.
- An-Ve_Log/: The log files of the reconstruction. This can be useful to assess convergence, resolution levels, etc.
  
- Re-Se/: The sensitivity maps stored as NIFTI images.
- Re-Se_Sn/: Snapshots of the sensitivity estimation process.
- Re-Se_Log/: The log files of the sensitivity estimation.

- Parsing_Sn: Snapshots of the raw data parsing.
- Parsin_Log/: The log files of the data parsing step.

## Contributors
Yannick Brackenier & Lucilio Cordero-Grande

## License



