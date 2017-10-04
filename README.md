# MK15_phosphoinositide_pathway_model
ODE model of the phosphoinositide pathway


To run this code you will need to install R. ( https://www.r-project.org/ )

The code was created in windows 10. It should run on Mac, but the command windows() must be changed to quartz(). The code was not tested in Linux. 

The content of each folder should be fully downloaded to the same disk folder for the code execution to perform without errors.

In the 2_Model folder, the model is presented in two configurations Apical and Basolateral. Files : MK15_apical.R and MK15_basolateral.R. The files map_matrix_MK6_perfect.txt, node_coords.txt, vertex_coords_MK6_3D.txt and vertex_coords_MK6_perfect.txt are required to build the model map with igraph.

In the 4_Global_sensitive_analysis folder, the .RData files are where the data from Monte Carlo search is stored. In the Global_sens_analysis_for_MK15_5_10_20_50_100.R file the initial chunk of code was used to perform the Monte Carlo search and should not be executed. 

You will find the same situation in the 5_Monte_Carlo_search folder. The .RData files are where the data from Monte Carlo search is stored. In the MK15_Monte_Carlo_score.R file the initial chunk of code was used to perform the Monte Carlo search and should not be executed. 

Any questions, please use the mail dvolivenca@fc.ul.pt
