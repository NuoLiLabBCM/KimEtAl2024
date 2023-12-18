# KimEtAl2024
Data analysis code for Kim, Daie, Li et al A combinatorial neural code for long-term motor memory

These are scripts to analyze key results in the paper. You should download 'Source_codes' folder with codes inside it to run the scripts.
Brief descriptions about each script below:

STEP1: Identification of co-registered neurons across imaging sessions and saving mat files (deconvolved dF/F0 PSTH of individual trial types) 
We used Suite2P package (Pachitariu et al) to perform motion correction and automatic ROI identification.
We used CellReg package (Sheintuch et al 2017) to find co-registered neurons from the same imaging field of view across sessions.
We used OASIS package (Friedrich et al 2017) to calculated deconvolved dF/F0.

STEP2: Sorting out reliable co-registered neurons

STEP3: Estimating Sample CD, Delay CD, Response CD.
We estimated a coding direction (CD) which is maximally separating different trials types (anterior pole vs. posterior pole: Sample CD; lick left versus lick right: Delay CD and Response CD) in an epoch-specific manner. We used 
