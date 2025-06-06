# SquyresNewmanBiofilmAnalysis
Biofilm image processing and modeling code used in Squyres and Newman, “Real-time high-resolution microscopy reveals how single-cell lysis shapes biofilm matrix morphogenesis”. 

[Preprint Available Here](https://www.biorxiv.org/content/10.1101/2024.10.13.618105v1)

Each processing step is briefly described below; see manuscript for method details. Relevant scripts and functions are listed, with the main function indicated for each. Dependencies indicates necessary software from others. 

### Image preprocessing: 

Image file pre-processing for confocal Z stack time lapses. Configured for .ims files from an Andor Dragonfly confocal. Merges multiple files with different names and potentially different Z stack heights, combines, downsamples to 8-bit, optionally resizes, and registers. Memory efficient: only loads pairs of images at a time at the cost of cropping registered images. Portions of this function are adapted from the MATLAB bioformats reader and from Guizar-Sicairos et al. Opt. Lett. 2008.

**Main method**:  imagePreprocessing.m

**Additional functions**: dftreg3D.m

**Dependencies**: [MATLAB Bioformats](https://www.openmicroscopy.org/bio-formats/downloads/), [Natural-Order Filename Sort](https://www.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort)

### Biofilm segmentation:

Negative stained biofilms were segmented using [Cellpose-SAM](https://www.biorxiv.org/content/10.1101/2025.04.28.651001v1) following the [documentation](https://cellpose.readthedocs.io/en/latest/). Training was performed using a set of manually annotated images, which are included here. Biofilm segmentation was also performed by thresholding, as described in the next section. 

**Training files**: TrainingData.zip

### Lysis depth analysis:

Identifies lysis events and measures their depth in the biofilm based on whole-biofilm segmentation. Biofilm segmentation can either be performed using the threshold-based method provided, or masks can be loaded that were generated through a separate method (e.g. Cellpose). Lysis events can either be detected using automatic peak finding or, as a verification method, can be manually identified in the provided interactive interface. Lysis depths are compared to a random simulation and results are plotted. 

**Main method**:  getLysisDepths.m

**Additional functions**: biofilmSeg_Threshold.m, lysisCaller_Automated.m, lysisCaller_Interactive.m

**Dependencies**: [MATLAB Bioformats](https://www.openmicroscopy.org/bio-formats/downloads/)

### Modeling:

Models comparing the consequences of patterned cell lysis (cell lysis that follows the experimentally observed probability distribution) with uniform cell lysis (unregulated, probability of cell lysis is equal everywhere). The lysisCoupling models determine how the cumulative amount of cell lysis is coupled to biofilm growth in each case, and the matrixDistribution models determime the final eDNA matrix distribution that results in each case. In both cases, we have included both both a simple 1D model, an equation for each case that captures the underlying principle of these relationships, and a 3D model for comparison to real data. We have also included a visualization of simulated cell lysis events that result from both models. 

**Lysis coupling models**: lysisCoupling1D.m, lysisCoupling3D.m

**eDNA matrix distribution models**: matrixDistribution1D.m, matrixDistribution3D.m

**Model visualization**: viewModels.m


