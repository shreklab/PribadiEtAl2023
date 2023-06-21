======================PribadeEtAl2023 README=============================

The code in this repository relates to the following manuscript:
"Dopamine signaling regulates predator-driven changes in Caenorhabditis elegans’ egg laying behavior"
Amy Pribadi, Michael A Rieger, Kaila Rosales, Kirthi C Reddy, Sreekanth H Chalasani

Note: The final figure numbering/lettering scheme has changed for the final submission
Because changing all the code would have been laborious, code has been kept intact but please use the key at the bottom for figure/code correspondencies

Dependencies:


R analyses depend on the car and multcomp packages for omnibus and post hoc comparisons.
Additionally the readxl package is used to read Excel data files.

The original data spreadsheet and planned post hoc comparisons sheets used by this code are included for reproducibility, but the code for each set of analyses has a preamble section where these sheets are defined and this can be modified as desired.


Root directory:

alwaysLoad.R : Includes a number of wrapper functions for preparing code outputs as well as initialization options and determining plot axes.
salkcolorpalette.R : Code that generates Salk approved color palette.
plotsalkcolorpalette.R : Generates a plot of the color palette for manual selection of colors in plots. 


Egg assays:
For egg assays, "_analysis.R" files compute:
		P(off), # of eggs, Mean, Lower quartile, Upper quartile, CV (Coefficient of variation) with associated omnibus and post hoc comparisons.
		
Figure Numbering Discrepancies
FIG 1\
	Figure1*
FIG 1-S1\
	Figure1-S1B*.R 						Figure 1-figure supplement 2
	Figure1-S1C_plot.R					Figure 1-figure supplement 3
	Figure1-S1D*.R						Figure 1-figure supplement 4
FIG 1-S2\			
	Figure1-S2*.R						Figure 1-figure supplement 5
FIG 1-S3\			
	Figure1-S3*.R						Figure 1-figure supplement 6
FIG 2\			
	Figure2*.R							No labeling discrepancies with presented figure.
FIG 3\			
	Figure3*.R							No labeling discrepancies with presented figure.
FIG 4\			
	Figure4B-C*.R						No labeling discrepancies with presented figure.
	Figure4D-G*.R						Figure 4d (original "d-g" combined as a single panel label).
FIG 5\			
	Figure5*.R							No labeling discrepancies with presented figure.
FIG 6\			
	Figure6A-C_analysis.R				Figure 6a & 6b
	Figure6A_plot.R						Figure 6a plot
	Figure6B_plot.R						Unused, log2 fold change in baseline P(off), now represented below 6a
	Figure6C_plot.R						Figure 6b plot
	Figure6D-F_analysis.R				Figure 6c & 6d
	Figure6D_plot.R						Figure 6c
	Figure6E_plot.R						Unused, see Figure6B_plot.R
	Figure6F_plot.R						Figure 6d
	Figure6G-I_analysis.R				Figure 6e & 6f
	Figure6G_plot.R						Figure 6e
	Figure6H_plot.R						Unused, see Figure6B_plot.R
	Figure6I_plot.R						Figure 6f
FIG	6-S4\
	Figure6-S4*.R						Figure 6-figure supplement 1
FIG 7\
	Figure7B-D_analysis.R				Figure 7b and 7c
	Figure7B_plot.R						Figure 7b
	Figure7C_plot.R						Unused as in Figure6B_plot.R above
	Figure7D_plot.R						Figure 7c
FIG 8\
	Figure8A_dopreceptors_alignment.m	Note: this is a MATLAB code file to generate Figure 8a
	Figure8B-D_analysis.R				Figure 8b & 8c
	Figure8B_plot.R						Figure 8b
	Figure8C_plot.R						Unused, as in Figure6B_plot.R above
	Figure8D_plot.R						Figure 8c
FIG 9\
	Figure9A-C_analysis.R				Figure 9a & 9b
	Figure9A_plot.R						Figure 9a
	Figure9B_plot.R						Unused, as in Figure6B_plot.R above
	Figure9C_plot.R						Figure 9b
	Figure9D-F_analysis.R				Figure 9c & 9d
	Figure9D_plot.R						Figure 9d
	Figure9E_plot.R						Unused, as in Figure6B_plot.R above
	Figure9F_plot.R						Figure 9f