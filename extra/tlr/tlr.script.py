###
### This script investigates if delta pt can be explained by ribosomal footprints.
###


# take all pt that have a consistent n > 2 signficant FC

# filter data. Remove pt/mRNA/RF whose standard error of the mean is larger than 30%.

# x.1. Scatter plot of FC_pt vs FC_mRNA. This figure would reveal how well mRNA explains pt changes.

# x.2. Scatter plot of FC_pt vs FC_RF. This figure would reveal how well RF explains pt changes.

# x.3 Scatter plot of FC_pt vs FC_mRNA + FC_RF. This figure would reveal how well an equal linear model would explain pt changes.

# x.4. Scatter plot of model: FC_pt = w_1 FC_mRNA + w_2 (FC_RF - FC_RF_0).
# For this section I need to obtain FC_RF_0 from linear regression on FC_mRNA vs FC_RF, from all genes.
# w_1 is defined within the range [0,1].
# w_2 is defined within the range [-1,1].

# x.5. visualization of w_1 vs w_2 to distinguish TCR from TLR.

# x. A sanity check for the script is the capture of RF causing pausing on rps10-like.
