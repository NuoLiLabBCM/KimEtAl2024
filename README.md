# KimEtAl2024
Data analysis code for Kim, Daie, Li et al A combinatorial neural code for long-term motor memory
Updated 18 Dec 2023

These are scripts to analyze key results in the paper. You should download 'Source_codes' folder with codes inside it to run the scripts.
Brief descriptions about each script below:

STEP1: Identification of co-registered neurons across imaging sessions and saving mat files (deconvolved dF/F0 PSTH of individual trial types) 
We used Suite2P package (Pachitariu et al) to perform motion correction and automatic ROI identification.
We used CellReg package (Sheintuch et al 2017) to find co-registered neurons from the same imaging field of view across sessions.
We used OASIS package (Friedrich et al 2017) to calculated deconvolved dF/F0.

STEP2: Sorting out reliable co-registered neurons
We further sorted out 'reliable' co-registered neurons for further analyses. We computed Pearson's correlation between concatenated (lick left and lick right) mean PSTHs from 1st/2nd halves. We only used reliable (R>=0.5) from at least one imaging session. 

STEP3: Estimating Sample CD, Delay CD, Response CD.
We estimated a coding direction (CD) which is maximally separating different trials types (anterior pole stimulus vs. posterior pole stimulus: Sample CD; lick left versus lick right: Delay CD and Response CD) in an epoch-specific manner. We used Li/Daie et al 2016 method to determine CD.
Based on the estimated CD, we built an epoch-specific linear classifier to decode different trial types. We used Chen/Kang et al 2021 method to build decoder.
