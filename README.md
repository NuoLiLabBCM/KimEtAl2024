# KimEtAl2024
Data analysis coder for "A combinatorial neural code for long-term motor memory." Jae-Hyun Kim, Kayvon Daie, Nuo Li. Nature 2024
Updated 17 Sep 2024 by Jae-Hyun Kim

All pre-processed imaging dataset (233 imaging sessions from 14 mice) were deposited at DANDI (https://doi.org/10.48324/dandi.001188/0.240912.1925).
These dataset is output of Suite2P package. Further running of CellReg and OASIS are required to find co-registered neurons across sessions and calculate deconvolved dF/F0, respectively.

These are MATLAB scripts to analyze key results in the paper. 
Before running codes, make sure:
1. Download 'Source_codes' folder with codes inside it to run the scripts.
2. Download the large data associated with these codes from the Zenodo links.
- Imaging dataset 1 (2023_07_07) https://zenodo.org/uploads/10402559
- Imaging dataset 2 (2023_07_25) https://zenodo.org/uploads/10402610
3. Make '2023_07_07' folder and put imaging dataset 1
4. Make '2023_07_25' folder and put imaging dataset 2
5. Run STEP1, STEP2, STEP3 scripts one after another.

Brief descriptions about each script below:

STEP1: Identification of co-registered neurons across imaging sessions and saving mat files (deconvolved dF/F0 PSTH of individual trial types) 
We used Suite2P package (Pachitariu et al 2016, v0.1.2) to perform motion correction and automatic ROI identification.
We used CellReg package (Sheintuch et al 2017, v1.4.5) to find co-registered neurons from the same imaging field of view across sessions.
We used OASIS package (Friedrich et al 2017, version of this paper directly) to calculate deconvolved dF/F0.

STEP2: Sorting out reliable co-registered neurons
We further sorted out 'reliable' co-registered neurons for further analyses. We computed Pearson's correlation between concatenated (lick left and lick right) mean PSTHs from 1st/2nd halves. We only used reliable (R>=0.5) from at least one imaging session. 

STEP3: Estimating Sample CD, Delay CD, Response CD.
We estimated a coding direction (CD) which is maximally separating different trials types (anterior pole stimulus vs. posterior pole stimulus: Sample CD; lick left versus lick right: Delay CD and Response CD) in an epoch-specific manner. We used Li/Daie et al 2016 method to determine CD.
Based on the estimated CD, we built an epoch-specific linear classifier to decode different trial types. We used Chen/Kang et al 2021 method to build decoder.

