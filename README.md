# A data-driven framework for structure-property correlation in ordered and disordered cellular metamaterials
This repository provides the data and code for a data-driven framework focusing on structure-property correlation in cellular metamaterials. The proposed framework enables the prediction of macroscopic properties, but more importantly, also reveals their connection to key morphological characteristics, as identified by the integration of machine-learning models and interpretability algorithms.

### This paper is currently under revision with PNAS, and the link of the paper is coming soon.

## Table of Contents
- [Data](#data)
- [Code Overview](#code-overview)
- [Prerequisites](#prerequisites)
- [Contributors](#contributors)
- [Contact](#contact)
  
## Data
### Microstructure Data

The microstructure data contains 1646 different tessellations, including both ordered and disordered. For visualization, each tessellation is represented by the corresponding node and connection (in the *Tessellation Dataset* folder), and two demos are provided to display the tessellation and/or microstructure.

> `Tessellation_Demo.m` Display of tessellation for a certain sample in the dataset.

> `Microstructure_Demo.m` Display of microstructure for a certain sample at a specific relative density in the dataset.

### Structure-Property Data

The structure-property data contains 42 microstructural features and also the corresponding macroscopic property for 1646 different tessellations at 5 different relative densities, which in total consists of 8230 microstructures with 43 properties.

> `Structure-Property-Data.csv` Each column represents a property with name listed in the header. Each five rows denote a microstructure at different relative densities in the same order with sample index in *Tessellation Dataset*.

## Code Overview
1. `Virtual microstructure generation` Generation of cellular metamaterial microstructure dataset.
2. `Feature and property measurement` Measurement of feature and property for cellular metamaterials in the dataset.
3. `Data Preparation` Converting the original mechanical data and generation of formal data files for future data-driven techniques.
4. `Random Forest and SHAP` Model of random forest for property prediction and utilization of SHAP for structure-property correlation.
5. `Generalized Additive Model` Model of GAM for property prediction and structure-property correlation.

## Prerequisites
- Matlab (R2020a or later, full toolbox installation recommended)
- Simulia Abaqus (2021)
- Python (3.8 or later)
- R (3.6.1 or later)


## Contributors
Shengzhi Luan, Enze Chen, Joel John, Stavros Gaitanaros

## Contact
For any further information, please feel free to contact us through email: stavrosg@jhu.edu.

![](./Figure/Framework.png)
