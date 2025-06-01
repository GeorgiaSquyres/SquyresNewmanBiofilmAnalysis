# SquyresNewmanBiofilmAnalysis
Biofilm image processing and modeling code used in Squyres and Newman, “Real-time high-resolution microscopy reveals how single-cell lysis shapes biofilm matrix morphogenesis”. 

[Preprint Available Here](https://www.biorxiv.org/content/10.1101/2024.10.13.618105v1)

Each processing step is briefly described below; see manuscript for method details. Relevant scripts and functions are listed, with the main function indicated for each. Dependencies indicates necessary software from others. 


### Image preprocessing: 

Image file pre-processing for confocal Z stack time lapses. Configured for .ims files from an Andor Dragonfly confocal. Merges multiple files with different names and potentially different Z stack heights, combines, downsamples to 8-bit, optionally resizes, and registers. Memory efficient: only loads pairs of images at a time at the cost of cropping registered images. Portions of this function are adapted from the MATLAB bioformats reader and from Guizar-Sicairos et al. Opt. Lett. 2008.

**Main method**:  imagePreprocessing.m

**Additional functions**: dftreg3D.m

**Dependencies**: [MATLAB Bioformats](https://www.openmicroscopy.org/bio-formats/downloads/), [Natural-Order Filename Sort](https://www.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort)
