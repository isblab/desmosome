# Intended to find "blob" wise centers (for randomization) and the "blob"-wise fit of the 1XM9 structure
# Opening a segmentation map also loads an extra density map some times which is not desirable
# hence opening the first segmentation map twice to later delete the first opened seg

open maps/pkp_separate.mrc
volume #0 level 0.006
open segments/sub_pkp_segments/pkp_individual_1.seg
open segments/sub_pkp_segments/pkp_individual_1.seg
open segments/sub_pkp_segments/pkp_individual_2.seg
open segments/sub_pkp_segments/pkp_individual_3.seg
open segments/sub_pkp_segments/pkp_individual_4.seg
open segments/sub_pkp_segments/pkp_individual_5.seg
open segments/sub_pkp_segments/pkp_individual_6.seg
open segments/sub_pkp_segments/pkp_individual_7.seg
close #1

# generate a separate map for all the blobs by successively masking them with the segments
mask #0 #2 modelId 9
mask #0 #3 modelId 10
mask #0 #4 modelId 11
mask #0 #5 modelId 12
mask #0 #6 modelId 13
mask #0 #7 modelId 14
mask #0 #8 modelId 15

# close all the un-required maps
close #0 #2 #3 #4 #5 #6 #7 #8

# measure the centers of all the blobs -> mark with a marker and get the coordinates of the marker
# the measure command returns the center in grid coordinates, the marker measurement is in the required coordinates
measure center #9 mark true modelId 16
getcrd #16
measure center #10 mark true modelId 17
getcrd #17
measure center #11 mark true modelId 18
getcrd #18
measure center #12 mark true modelId 19
getcrd #19
measure center #13 mark true modelId 20
getcrd #20
measure center #14 mark true modelId 21
getcrd #21
measure center #15 mark true modelId 22
getcrd #22

# close unnecessary parts
close #16-22

# load the 1xm9 pdb (check the path) and generate a molmap
open ../homology/input_data/pdb_files/1xm9.pdb
molmap #0 32 gridSpacing 6 sigmaFactor 0.187 modelId 1
volume #1 sdLevel 1
close #0

# get the center of the molmap
measure center #1

# Sequentially, bring each blob to the origin (i.e. center of mass centered at the origin)
# the i, j, k values are taken from the reply log as returned by measure center
volume #1 originIndex 23.10,22.04,19.70
volume #15 originIndex 6.40,5.49,7.91
volume #14 originIndex 5.40,6.43,7.90
volume #13 originIndex 6.90,6.23,5.89
volume #12 originIndex 6.21,4.80,7.47
volume #11 originIndex 6.34,5.04,6.98
volume #10 originIndex 6.98,7.10,4.96
volume #9 originIndex 7.16,6.13,5.80

# fit the molmap to each of the blob densities
fitmap #9 #1 maxSteps 6000
fitmap #10 #1 maxSteps 6000
fitmap #11 #1 maxSteps 6000
fitmap #12 #1 maxSteps 6000
fitmap #13 #1 maxSteps 6000
fitmap #14 #1 maxSteps 6000
fitmap #15 #1 maxSteps 6000