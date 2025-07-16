# MITF_MD
Molecular dynamics experiments for MITF

![image](https://github.com/user-attachments/assets/0fd97db3-638b-4e26-b3b8-2508334c49e2)

Compound **7** shown in blue (in orientation 2)
Compound **8** shown in orange (in orientation 1)

## Contents

`Source_Data/`: Helix trace files

`figs/`: All figures for the paper as well as some supplementary figures

`input_structures/`: Input structures for MD simulations

`representative_structures/`: Extracted structures from last frame of each MD simulation

`mitf_analyses.R`: All computational analyses (incl statistical analyses) and figures for the paper

`utils.R`: Any utility functions required for `mitf_analyses.R`, mostly for averaging over the helix traces

`calc_distances.cpptraj`: Sample script to calculate helix distances throughout MD trajectory

`calc_mindist.trajin`: Sample script to calculate the minimum distance between the protein complex and periodic images

`image.trajin`: Sample trajectory imaging script


