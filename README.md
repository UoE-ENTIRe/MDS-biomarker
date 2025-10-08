# Mean Diffusion Slope (MDS) Imaging Biomarker in Glioblastoma

This repository contains the computational code and analysis pipeline supporting the validation of a cancer imaging biomarker for glioblastoma using the Imaging Biomarker Roadmap framework proposed by O'Connor et al.

## Overview

This work implements a quantitative imaging approach for characterizing tumor boundaries in glioblastoma using apparent diffusion coefficient (ADC) maps derived from diffusion-weighted magnetic resonance imaging. The analysis pipeline includes:

- Gaussian infiltration front simulation models
- 2D vector analysis of tumor boundary profiles
- Multi-spectral image data preprocessing
- Quantitative biomarker validation against histological reference standards

## Repository Structure

- `ADC_Transition_IB_code/`: Primary analysis code including simulation models and vector analysis
  - `GaussianfrontSim.m`: Gaussian infiltration front simulation 
  - `Vector_Analysis_2D.m`: 2D tumor boundary vector analysis
  - `Utilities_Code/`: Supporting functions including NIfTI image processing tools
- `Data matrix/`: Data preprocessing utilities and processed datasets
- `Vector_Analysis_Results/`: Analysis output directory
- `Vector_Analysis_Results_Sim/`: Simulation analysis output directory

## Requirements

### MATLAB Toolboxes
- Image Processing Toolbox
- Fuzzy Logic Toolbox

### External Dependencies
- NIfTI Tools: Jimmy Shen (2023). Tools for NIfTI and ANALYZE image
- Error bars: Jean-Yves Tinevez (2023). errorbarxy  
- Image overlay: Matt Smith (2023). imoverlay
- Subplot layout: Felipe G. Nievinski (2023). subtightplot

## Usage

The analysis pipeline consists of three main components:

1. **Data Preprocessing** (`Data matrix/Data_preprocessing.m`): Harmonizes multi-spectral imaging data and applies normalization procedures.

2. **Simulation Analysis** (`ADC_Transition_IB_code/GaussianfrontSim.m`): Generates synthetic infiltration models with configurable parameters for validation studies.

3. **Vector Analysis** (`ADC_Transition_IB_code/Vector_Analysis_2D.m`): Performs quantitative analysis of tumor boundary profiles using 2D sampling vectors.

## Data Format

The pipeline processes co-registered multi-spectral imaging data including:
- T1-weighted and T2-weighted MRI
- Diffusion tensor imaging (FA, ADC, AD, RD)
- Perfusion imaging (mbASL)
- Histological validation data (H&E, immunohistochemistry)

## Biomarker Validation Framework

This implementation follows the technical performance assessment criteria outlined in the imaging biomarker roadmap, focusing on:
- Measurement repeatability and reproducibility
- Biological correlation with histological markers
- Technical validation through simulation studies

## Citation

[Publication details to be added upon acceptance]

## License

This work is released under the MIT License. See LICENSE file for details.

## Contact

Dr. Gerard Thompson  
University of Edinburgh  
gerard.thompson@ed.ac.uk

Dr. Antoine Vallatos  
University of Glasgow  
antoine.vallatos@glasgow.ac.uk
