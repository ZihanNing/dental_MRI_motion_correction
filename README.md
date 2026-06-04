# Motion-Robust Dental MRI for Paediatric Dental Trauma

This repository contains the main reconstruction and workflow code accompanying our paper:

**Motion-Robust Dental MRI for Imaging of Paediatric Dental Trauma**

In the paper, we describe the proposed method as a **region-adaptive motion correction reconstruction tailored to dental MRI**. 

The proposed automated region-adaptive motion-correction reconstruction identifies maxillary and mandibular regions, performs region-specific motion estimation and correction, and fuses the corrected regions into a single image for clinical review. 

The core workflow is implemented in [batch_dental_multiple.m](https://github.com/ZihanNing/dental_MRI_motion_correction/blob/main/batch_dental_multiple.m), while the main reconstruction routines remain in this Git repository.

## Repository Contents

- `batch_dental_multiple.m`: end-to-end batch workflow for case processing
- `deployRecon_dental_SENSE.m`: baseline SENSE reconstruction
- `deployRecon_dental_MoCo.m`: region-adaptive motion-correction reconstruction
- `zihan_tools/image_fusion.m`: fusion of upper-jaw and lower-jaw MoCo reconstructions and other dependent scripts
- `Python/`: preprocessing, intensity normalization, and mask restoration scripts for nnUNet-based segmentation
- `brackenier-tools/` and `dhcp-repo-release07/`: code dependencies included with this repository for motion-correction reconstruction

## Release Plan

The full open-source release is split across two locations:

- **GitHub repository**: MATLAB workflow, reconstruction code, and utility scripts
- **Zenodo**: demo raw data and trained nnUNetv2 models for head segmentation and teeth segmentation

Data and trained model download:

- Zenodo record: [https://zenodo.org/records/20544237](https://zenodo.org/records/20544237)

## Requirements

### MATLAB

The workflow has been run on **MATLAB R2019a** on our server.

At minimum, the current workflow expects MATLAB with NIfTI and image-processing functionality used by:

- `niftiinfo`, `niftiread`, `niftiwrite`
- `imgaussfilt3`
- `imregtform`
- `imwarp`
- `imref3d`
- `bwdist`

In practice, **Image Processing Toolbox** is required for the segmentation post-processing and image-fusion steps.

Some reconstruction code paths in this repository also use GPU-aware MATLAB arrays (`gpuArray`). A GPU is not strictly required for every script, but it may be important for practical reconstruction speed depending on your setup.

### Python

The segmentation preprocessing and mask restoration scripts in `Python/` require Python packages including:

- `numpy`
- `scipy`
- `nibabel`

These scripts are called from MATLAB in `batch_dental_multiple.m`.

### nnUNetv2

To run the segmentation-enabled workflow, you will need **nnUNetv2** installed in a Python environment.

Essential links:

- nnUNet repository: [MIC-DKFZ/nnUNet](https://github.com/MIC-DKFZ/nnUNet)
- nnUNetv2 installation and setup: [official installation guide](https://github.com/MIC-DKFZ/nnUNet/blob/master/documentation/getting-started/installation-and-setup.md)
- PyTorch installation: [PyTorch local installation guide](https://pytorch.org/get-started/locally/)

The current batch script assumes:

- a Conda environment name such as `nnunetv2`
- `nnUNetv2_predict` is available in that environment
- nnUNet paths are configured via `nnUNet_raw`, `nnUNet_preprocessed`, and `nnUNet_results`

In [batch_dental_multiple.m](./batch_dental_multiple.m), these are currently set through the local variables `CONDA`, `ENVNAME`, and `NNUNET_BASE`. You will likely need to edit these paths for your system before running the workflow.

### Trained Model Weights

The repository does **not** store the trained nnUNet model weights directly. Instead, the trained models for:

- teeth segmentation
- head segmentation

are distributed together with the demo raw data on Zenodo:

- [https://zenodo.org/records/20544237](https://zenodo.org/records/20544237)

After downloading them, place them into the nnUNet results structure expected by your local nnUNetv2 installation.

## Installation

Clone the repository and add it to your MATLAB path. The workflow already uses:

```matlab
addpath(genpath(pwd))
```

so the included helper code under this repository is made available automatically when launched from the project root.

If you are preparing a fresh environment, make sure the following are ready before running the batch workflow:

1. MATLAB R2019a or a compatible MATLAB release
2. Image Processing Toolbox
3. A Python environment with `numpy`, `scipy`, and `nibabel`
4. nnUNetv2 installed and callable from that environment
5. The released demo raw data and head/teeth nnUNet model weights downloaded from Zenodo

## Data Layout

The batch workflow expects case folders under `Studies-deploy/`, for example:

```text
Studies-deploy/
  1/
    *.dat
  2/
    *.dat
```

During processing, outputs are written into folders such as:

- `An-Aq/`: baseline reconstruction and segmentation-related intermediates
- `An-Ve/`: motion-corrected reconstructions and fused images
- `An-Aq/seg/`: restored head and teeth masks used by the workflow

## Quick Test with Demo Cases

Demo raw data and trained nnUNetv2 models are provided through Zenodo so users can test the released workflow without preparing their own dataset first:

- [https://zenodo.org/records/20544237](https://zenodo.org/records/20544237)

Recommended quick-start:

1. Clone this GitHub repository and start MATLAB from the repository root.
2. Download and unpack the demo raw data and trained nnUNetv2 model folders from Zenodo.
3. Place the unpacked demo case folders inside `Studies-deploy/`. The expected layout is:

```text
Studies-deploy/
  1/
    *.dat
  2/
    *.dat
```

4. Install nnUNetv2 locally by following the official nnUNetv2 installation instructions.
5. Locate the trained nnUNetv2 model folders downloaded from Zenodo:

```text
Dataset003_Dental_autolandmark_PDwMPRAGET2w
Dataset004_Dental_headseg_PDT1T2
```

6. Place the downloaded model folders under your local `nnUNet_results` directory. For example, if your nnUNet workspace is `/path/to/nnUNet`, the expected layout is:

```text
/path/to/nnUNet/
  nnUNet_raw/
  nnUNet_preprocessed/
  nnUNet_results/
    Dataset003_Dental_autolandmark_PDwMPRAGET2w/
    Dataset004_Dental_headseg_PDT1T2/
```

7. Open [batch_dental_multiple.m](https://github.com/ZihanNing/dental_MRI_motion_correction/blob/main/batch_dental_multiple.m) in MATLAB.
8. Update the demo case settings near the top of the script if needed:

```matlab
rootFolder  = './Studies-deploy';
caseList    = [1];
seqSelect   = {'MPRAGE','T2wSPACE','PDwSPACE'};
```

9. Update the local Python and nnUNet paths in the script:

```matlab
CONDA = '/path/to/anaconda3/bin/conda';
ENVNAME = 'nnunetv2';
NNUNET_BASE = '/path/to/nnUNet';
```

`NNUNET_BASE` should point to the folder that contains `nnUNet_raw`, `nnUNet_preprocessed`, and `nnUNet_results`. The batch script sets these nnUNet environment variables automatically from `NNUNET_BASE` before calling `nnUNetv2_predict`.

10. Run the batch workflow from MATLAB:

```matlab
batch_dental_multiple
```

For a minimal test, start with a single case:

```matlab
rootFolder  = './Studies-deploy';
caseList    = [1];
seqSelect   = {'MPRAGE','T2wSPACE','PDwSPACE'};
```

## What the Batch Workflow Does

For each selected case and sequence, `batch_dental_multiple.m` performs:

1. baseline SENSE reconstruction
2. preprocessing for nnUNet inference
3. teeth segmentation
4. head segmentation
5. restoration of masks to the original resolution and field of view
6. landmark extraction from the teeth mask
7. region-specific motion-corrected reconstruction
8. fusion of upper-jaw and lower-jaw motion-corrected images

The fused image is saved in `An-Ve/` together with the other reconstruction outputs.

## Notes on Current Configuration

This repository currently includes some environment-specific paths in the MATLAB scripts, for example Conda and nnUNet installation locations. These should be treated as templates and adapted locally before use.

The present workflow has been developed and tested in our server environment. If you are porting it to a new machine, the main things to check first are:

- MATLAB version and toolbox availability
- Python environment activation
- nnUNetv2 installation and path variables
- location of downloaded trained weights
- case folder layout under `Studies-deploy/`

## Acknowledgements

This repository builds on included code from:

- `brackenier-tools`
- `dhcp-repo-release07`

We thank the original contributors to those components.

## Citation

If you use this code, please cite the accompanying paper:

**Motion-Robust Dental MRI for Imaging of Paediatric Dental Trauma**

Citation details will be added here upon publication.

## Contributors

Zihan Ning and collaborators.

## License

This repository is released under the **Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0)**, unless otherwise stated.

This means that public use, sharing, and modification are allowed with appropriate attribution, but commercial use is not permitted.

See [LICENSE](./LICENSE) for details. Third-party code and dependencies included in this repository, including `brackenier-tools` and `dhcp-repo-release07`, remain subject to their own original license terms.
