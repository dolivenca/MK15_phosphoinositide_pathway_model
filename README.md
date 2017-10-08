To run this code you will need to install R. ( https://www.r-project.org/ )

The code was created in windows 10. It should run on Mac, but the command windows() must be changed to quartz(). The code was not tested in Linux. 

The content of each folder should be fully downloaded to the same disk folder for the code execution to perform without errors.

In the 2_Model folder, the model is presented in two configurations Apical and Basolateral. Files : MK15_apical.R and MK15_basolateral.R. The files map_matrix_MK6_perfect.txt, node_coords.txt, vertex_coords_MK6_3D.txt and vertex_coords_MK6_perfect.txt are required to build the model map with igraph.

In the 4_Global_sensitive_analysis folder, the .RData files are where the data from Monte Carlo search is stored. In the Global_sens_analysis_for_MK15_5_10_20_50_100.R file the initial chunk of code was used to perform the Monte Carlo search and should not be executed. 

You will find the same situation in the 5_Monte_Carlo_search folder. The .RData files are where the data from Monte Carlo search is stored. In the MK15_Monte_Carlo_score.R file the initial chunk of code was used to perform the Monte Carlo search and should not be executed. 

The fluxes and their parameters have a different notation in the code and in the paper.
Notation in the paper:
Fluxes
V->0 = Y->0
V->4 = Y->4
V->3 = Y->3
V0->3 = Y0->3 * PI^f0->3 * (PI3KII+PI3III)
V3->0 = Y3->0 * PI(3)P^f3->0 * (SYNJ+SAC1+MTMR)
V0->4 = Y0->4 * PI^f0->4 * PI4K 
V4->0 = Y4->0 * PI(4)P^f4->0 * (SYNJ+SAC1) 
V0->5 = Y0->5  * PI^f0->5 * PIKfyve 
V5->0 = Y5->0 * PI(5)P^f5->0 * (SYNJ+SAC1) 
V3->35 = Y3->35 * PI(3)P^f3->35 * PIKfyve 
V35->3 = Y35->3 * PI(3,5)P2^f35->3 * (SYNJ+SAC1+SAC3) 
V4->45 = Y4->45 * PI(4)P^f4->45 * PIP5KI 
V45->4 = Y45->4 * PI(4,5)P2^f45->4 * (SIOSS) 
V5->45 = Y5->45 * PI(5)P^f5->45 * PIP5KII 
V45->5 = Y45->5 * PI(4,5)P2^f45->5 * (SYNJ+TMEM55) 
V45->345 = Y45->345 * PI(4,5)P2^f45->345 * PI3KI 	
V345->45 = Y345->45 * PI(3,4,5)P3^f345->45 * PTEN 
V35->5 = Y35->5 * PI(3,5)P2^f35->5 * MTMR 
V34->3 = Y34->3 * PI(3,4)P2^f34->3 * INPP4 
V345->34 = Y345->34 * PI(3,4,5)P3^f345->34 * (SIOSS+SHIP2) 
V45-> = Y45-> * PI(4,5)P2
V0-> = Y0-> * PI
V4-> = Y4-> * PI(4)P
V345-> = Y345-> * PI(3,4,5)P3
V3-> = Y3-> * PI(3)P
V35-> = Y35-> * PI(3,5)P2
V5-> = Y5-> * PI(5)P
V34-> = Y34-> * PI(3,4)P2
V0->45 = Y0->45 * PI^f0->45 * (PI4K+PIP5KI) 
V4->34 = Y4->34 * PI(4)P^f4->34 * PI3KII 
V34->4 = Y34->4 * PI(3,4)P2^f34->4 * PTEN 
V45->0 = Y45->0 * PI(4,5)P2^f45->0 * SYNJ
Differential equations
dPI = V->0 + V3->0 + V4->0 + V5->0 + V45->0 - V0->3 - V0->4 - V0->5 - V0->45 - V0-> 
dPI3P = V->3  + V0->3 + V35->3 + V34->3 - V3->0 - V3->35 - V3->
dPI4P = V->4 + V0->4 + V45->4+ V32 - V4->0 - V4->45  - V24 - V4->
dPI5P = V0->5 + V35->5 + V45->5 - V5->0 - V5->45 - V5->
dPI35P2 = V3->35 - V35->5 - V35->3 - V35->
dPI45P2 = V4->45 + V5->45 + V345->45 + V0->45 - V45->4 - V45->5 - V45->345 - V45->0 - V45->
dPI34P2 = V345->34 + V4->34 - V34->4 - V34->3 - V34->
dPI345P3 = V45->345 - V345->45 - V345->34 - V345->

Notation in the code:
Fluxes		
V1 = gamma1
V16 = gamma16
V18 = gamma18
V2 = gamma2 * PI^f2_PI * pi_3KII_III
V3 = gamma3 * PI3P^f3_PI3P * (SYNJ_SAC1_MTMR)
V4 = gamma4 * PI^f4_PI * pi_4K
V5 = gamma5 * PI4P^f5_PI4P * (SYNJ_SAC1)
V6 = gamma6 * PI^f6_PI * pi_Kfyve
V7 = gamma7 * PI5P^f7_PI5P * (SYNJ_SAC1)
V8 = gamma8 * PI3P^f8_PI3P * pi_Kfyve
V9 = gamma9 * PI35P2^f9_PI35P2 * (SYNJ_SAC1_SAC3)
V10 = gamma10 * PI4P^f10_PI4P * pip_5KI
V11 = gamma11 * PI45P2^f11_PI45P2 * (SIOSS)
V12 = gamma12 * PI5P^f12_PI5P * pip_5KII
V13 = gamma13 * PI45P2^f13_PI45P2 * (SYNJ_TMEM55)
V14 = gamma14 * PI45P2^f14_PI45P2 * pi_3KI
V15 = gamma15 * PI345P3^f15_PI345P3 * PTEN
V17 = gamma17 * PI35P2^f17_PI35P2 * MTMR
V19 = gamma19 * PI34P2^f19_PI34P2 * INPP4
V21 = gamma21 * PI345P3^f21_PI345P3 * (SIOSS_SHIP2)
V22 = gamma_e * PI45P2
V23 = gamma_e * PI
V24 = gamma_e * PI4P
V25 = gamma_e * PI345P3
V26 = gamma_e * PI3P
V27 = gamma_e * PI35P2
V28 = gamma_e * PI5P
V29 = gamma_e * PI34P2
V30 = gamma30 * PI^f30_PI * pi_4K_pip_5KI
V31 = gamma31 * PI4P^f31_PI4P * pi_3KII
V32 = gamma32 * PI34P2^f32_PI34P2 * PTEN 
V33 = gamma33 * PI45P2^f33_PI45P2 * (SYNJ)
Differential equations
dPI = V1 + V3 + V5 + V7 + V33 - V2 - V4 - V6 - V23 - V30 
dPI3P = V18 + V2 + V9 + V19 - V3 - V8 - V26
dPI4P = V16 + V4 + V11 + V32 - V5 - V10 - V24 - V31
dPI5P = V6 + V17 + V13 - V7 - V12 - V28
dPI35P2 = V8 - V9 - V17 - V27
dPI45P2 = V10 + V12 + V15 + V30 - V11 - V13 - V14 - V22 - V33
dPI34P2 = V21 + V31 - V19 - V29 - V32
dPI345P3 = V14 - V15 - V21 - V25



Any questions, please use the mail dvolivenca@fc.ul.pt
