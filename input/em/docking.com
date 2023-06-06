# Intended to find the best docking of 1XM9 PKP structure to three of the PKP-layer density blobs
# This docked orientation is used for G_PKP

open maps/pkp_separate.mrc
volume #0 level 0.006
open maps/sub_pkp_segments/pkp_individual_1.seg
open maps/sub_pkp_segments/pkp_individual_2.seg
open maps/sub_pkp_segments/pkp_individual_4.seg
open maps/sub_pkp_segments/pkp_individual_6.seg
close #1

# generate a separate map for all the blobs by successively masking them with the segments
mask #0 #2 modelId 5
mask #0 #3 modelId 6
mask #0 #4 modelId 7

# close all the un-required models
close #0 #2 #3 #4

# load the initial-position manually-docked 1xm9 PDBs (the order is to just match the blobs)
open docking/docked_start_2_1xm9.pdb
open docking/docked_start_3_1xm9.pdb
open docking/docked_start_1_1xm9.pdb

# fit (map-map) for each of the molecules
fitmap #0 #5 maxSteps 6000 resolution 32 metric correlation
fitmap #1 #6 maxSteps 6000 resolution 32 metric correlation
fitmap #2 #7 maxSteps 6000 resolution 32 metric correlation

# write the new PDBs (the order here does not match the numbers in the start-PDBs
write #0 docking/docked_2_gpkp.pdb
write #1 docking/docked_3_gpkp.pdb
write #2 docking/docked_1_gpkp.pdb
