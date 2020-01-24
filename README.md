# spline-aplot
Experimenting with adaptive plotting of splines given in B-form. Implementation in MATLAB.

1) Start with reading 'Dahl_thesis.pdf', at least skim through sections 3.2-3.6. Read also 'chapterX.pdf' for precisions and an update on a few points relevant to the implementation of the methods.
2) The script 'introPlotting.m' shows how to use the Curve Fitting Toolbox and some of the auxiliary functions developed for this project.
3) The general adaptive algorithm (under development) is found in spline_aplot_dev.m. See the comments in this file for more details. The core subfunctions that are called are named 'C_...'.
4) Read and run test_aplot.m to experiment with the adaptive plotting methods. This file sets up the workspace with 5 specific test-curves and an array holding randomly generated splines using specified parameters. See 'testcurves.m' for specifics on the testcurves and to generate standard plots of the curves and their original control polygons.
5) The goal is to find a good combination of methods in the general algorithm and adequate threshold values, and thus provide an algorithm that performs well for all B-splines without user interaction. Initial results indicate a potential for combining the length-method with the length-ratio method. This combination also allows for some reuse of values, lowering its cost compared to other combinations.
