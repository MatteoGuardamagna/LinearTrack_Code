# LinearTrack_Code
 Repository of analysis code for linear probe + single cell recordings in CA1


LINEAR TRACK WORKFLOW: \
\
3 Directories are assumed: \
- Raw data\
- GLM Location \
- Processed Results (Dropbox, GitHub or whatever)\
\
\
0. COMPUTE GAMMA COUPLINGS 
In code_gammas 
- Basins and layer-resolved theta-gamma coupling AND layer-resolved gamma power 
- Gamma coupling to animal speed and locations on the track 


1. COMPUTE GLM\
In code_glm: \
- Compute_PhasePosition_GLM_1D\
\
2 COMPUTE CELL PROPERTIES\
In code_properties:\
- Move_Files\
- Compute_ThetaScores\
- Compute_PlaceFields\
- Compute_CellProperties \
- Define_ThetaLimits \
\
3 SPIKE TIMING MEASURES \
In code_spikes: \
- BurstAnalysis %Computes phase-spiking properties \
- CrossAnalysis % Computes pairwise-spike coordination \
- SingleTrial_PP % Computes measure for single trial phase precession \
- SpikeGamma_Analysis % Probability of spikes conditioned on dominating gamma \
- GammaPhase_Analysis % Probability of gamma dominating conditioned on theta phase \
\
4 PROCESS RESULTS AND PLOTS \
in code_plotting (should depend only on data available in Processed Results folder)\
\
- At the moment follows Guardamagna 2022 structure \
\
\
\
\
}


## References

```
@Article{Guardamagna2021,
  author    = {Guardamagna, Matteo and Stella, Federico and Battaglia, Francesco P.},
  title     = {Distinct hippocampal network states support theta phase precession and theta sequences in CA1},
  journal   = {bioRxiv},
  year      = {2021}
}
```
