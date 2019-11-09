# Matlab Programs

This repository contains a variety of *Matlab* functions and tools to aid in oceaographic or limnologic research.  These may include tools for data analysis and processing, simple modeling, and data visualization.  This repository will be continually evolving as I learn more and continue to build new tools and refine the existing ones.

I am in the process of adding usage examples and syntax for each function.  For some tools, there are included Matlab help-style comments (so you can use Matlab's **help** function).  If you have questions about a function use where I have not included help documentation or examples, please feel free to contact me with questions.

## Contents

Air_sea_interaction
├── Heat_Flux
│   ├── BCSM.m
│   ├── COARE
│   │   └── ...
│   ├── albedo.m
│   └── insolation.m
├── WDir_conversions
│   ├── cart2ocean.m
│   └── ocean2cart.m
└── Wind_Stress
Constants
└── g.m
General
├── binAvg.m
├── binCenters.m
├── convSmooth.m
├── driftCoord2map.m
├── findbetween.m
├── map2driftCoord.m
└── sectionAPE.m
Instrument_processing
├── Nortek_Signature
│   ├── beam2xyz_enu.m
│   ├── sigBeamMapping.m
│   ├── sigBeamMappingInv.m
│   ├── sigDeclinationCorrection.m
│   ├── sigEnsembleAvg.m
│   ├── sigMagCorrection2D.m
│   └── sigWavesProcess.m
└── Underway_CTD
Plotting
├── Colormaps
│   ├── PNWColors.m
│   └── ucar.m
├── Data_Cursor
│   ├── custom_data_cursor.m
│   └── data_cursor_datest.m
├── TS_plot.m
├── bg_patch.m
├── cinterp.m
├── make_gif.m
└── trackRibbon.m
Remote_Sensing
├── Polar_Pathfinder
│   ├── EASEgrid2map.m
│   ├── READ_icemotion.m
│   ├── Usage_Example
│   │   ├── Example.m
│   │   ├── coast.mat
│   │   ├── icemotion.grid.week.2014.39.n.v3.bin
│   │   └── icemotion_coords.mat
│   ├── map2EASEgrid.m
│   └── vec_EASEgrid2map.m
└── README.md
Surface_Gravity_Waves
    └── vectWavenum.m
