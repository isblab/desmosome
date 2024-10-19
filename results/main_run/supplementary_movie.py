import os
from chimerax.core.commands import run

# Names of proteins/domains for which we have created densities
prots = [
    "DP-N", "DP-S", "DSC", "DSG", "GPKP", "PG-C", "PG-N", "PG-S", "PKP-C", "PKP-N", "PKP-S"
]

if f"LPD_{prots[0]}.mrc" not in os.listdir():
    raise Exception(
        "Density maps not found."
        " Please change the directory to the"
        " one containing the density maps by"
        " running the following command in ChimeraX:"
        " cd /path/to/density/maps"
    )

prot_names = [
    "DP1-N", "DP1-S", "DSC1a", "DSG1a", "PKP1aG", "PG-C", "PG-N", "PG-S", "PKP1a-C", "PKP1a-N", "PKP1a-S"
]

threshold = [
    0.0167, 0.036, 0.008, 0.0117, 0.071, 0.0084, 0.0193, 0.063, 0.0035, 0.0151, 0.09
]

colors = [
    "#e31a1c", "#ef7678", "#1f78b4", "#6a3d9a", "#4daf4a",
    "#ffb366", "#994d00", "#ff7f00", "#95d293", "#377e35",
    "#4daf4a"
]

pd_transparency = 50 # transparency of localization probability density maps
total_models = 9000 # Upper limit of models to load out of 24016
step_model = 20 # Load every 20th model upto 9000

num_frames = int(total_models/step_model)
per_frame_angle = 360/num_frames
wait_for_movie = num_frames
framerate = 30
per_frame_move = 1
num_frames_move = 80

def dothis(command_list):
    for command in command_list:
        run(session, command)

# load localization probability density maps
for i, p in enumerate(prots):
    dothis(
        [
            f"open ./LPD_{p}.mrc",
            f"volume #{str(i+1)} step 1",
            f"volume #{str(i+1)} level {str(threshold[i])}",
            f"color #{str(i+1)} {colors[i]}"
        ]
    )
    i += 1

# set view
dothis(
    [
        f"move X 300 models #5-11",
        f"move X -300 models #4",
        f"move Y 300 models #1-2,6-8",
        "view all",
        "view name all_view",
        f"move X -300 models #5-11",
        f"move X 300 models #4",
        f"move Y -300 models #1-2,6-8"
    ]
)

# set scence
dothis(
    [
        "set bgcolor #2B2B2B",
        "view all",
        "lighting soft",
        f"transparency all {pd_transparency} target s",
    ]
)

# load RMF file
dothis(
    [
        "open ./sampcon_0_extracted_main_run.rmf3",
        f"rmf readtraj #12.1 first 1 last {total_models} step {step_model}",
        f"coordset slider #12.1",
        "show all; view all",
        "view name one_view",
        "hide all target m"
    ]
)

# set labels
for i, p in enumerate(prot_names):
    dothis(
        [
            f"label #{str(i+1)} color white text {p} font arial height 20 offset -10,-10,0"
        ]
    )

####### start movie recording #######
dothis(
    [
        "movie record"
    ]
)

# Pkp models
dothis(
    [
        "show #5,9-11,12 target m",
        "hide #12",
        "show #12.1/A-G",
        "hide #12.1::resolution==10.0",
        "close #5.2", # remove ghost Pkp label
        "move X 150 models #9-11.2",
        f"coordset #12.1 1,{total_models}",
        f"roll y {per_frame_angle} {num_frames} models #5,9-11,12.1",
        f"wait {wait_for_movie}",
        "coordset stop",
        "hide all target m"
    ]
)

# Pg models
dothis(
    [
        "show #6,7,8,12 target m",
        "hide #12",
        "show #12.1/H-K",
        "hide #12.1::resolution==10.0",
        "move X 150 models #6,7,8,12.2",
        f"coordset #12.1 1,{total_models}",
        f"roll y {per_frame_angle} {num_frames} models #6,7,8,12.1",
        f"wait {wait_for_movie}",
        "coordset stop",
        "hide all target m"
    ]
)

# Dp models
dothis(
    [
        "show #1,2,12 target m",
        "hide #12",
        "show #12.1/L-O",
        "hide #12.1::resolution==10.0",
        "move X 150 models #1,2,12.2",
        f"coordset #12.1 1,{total_models}",
        f"roll y {per_frame_angle} {num_frames} models #1,2,12.1",
        f"wait {wait_for_movie}",
        "coordset stop",
        "hide all target m"
    ]
)

# Dsg models
dothis(
    [
        "show #4,12 target m",
        "hide #12",
        "show #12.1/R-S",
        "hide #12.1::resolution==10.0",
        "move X 150 models #4,12.2",
        f"coordset #12.1 1,{total_models}",
        f"roll y {per_frame_angle} {num_frames} models #4,12.1",
        f"wait {wait_for_movie}",
        "coordset stop",
        "hide all target m"
    ]
)

# Dsc models
dothis(
    [
        "show #3,12 target m",
        "hide #12",
        "show #12.1/P-Q",
        "hide #12.1::resolution==10.0",
        "move X 150 models #3,12.2",
        f"coordset #12.1 1,{total_models}",
        f"roll y {per_frame_angle} {num_frames} models #3,12.1",
        f"wait {wait_for_movie}",
        "coordset stop",
        "hide all target m"
    ]
)

# Show all maps
dothis(
    [
        "show #1-11.1 target m",
        f"transparency all {pd_transparency} target s",
        f"move X  300 models #5-11",
        f"move X -300 models #4",
        f"move Y 300 models #1-2,6-8",
        f"view all_view {num_frames_move}",
        f"wait 120",
        f"show #1-11.2 target m",
        f"move X -50 models #3.2",
        f"move X -50 models #4.2",
        f"move X -25 models #1,2.2",
        f"wait 80",
        f"hide #1-11.2",
        f"move X -{300/num_frames_move} {num_frames_move} models #5-11",
        f"move X {300/num_frames_move} {num_frames_move} models #4",
        f"move Y -{300/num_frames_move} {num_frames_move} models #1-2,6-8",
        f"view one_view {num_frames_move}",
        f"wait 120",
        f"roll y {per_frame_angle} {num_frames} models #1-11.1",
        f"wait {num_frames}"
    ]
)

####### stop movie recording #######
dothis(
    [
        f"movie encode ./Desmosome_viz.mp4 framerate {framerate} quality highest format h264 qscale 1",
        "coordset stop"
    ]
)