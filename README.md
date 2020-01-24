# spline-aplot
Experimenting with adaptive plotting of splines given in B-form. Implementation in MATLAB.

1) Start with reading 'Dahl_thesis.pdf', at least skim through sections 3.2-3.6.
2) The script 'introPlotting.m' shows how to use the Curve Fitting Toolbox and some of the auxiliary functions developed for this project.
3) The core functions are named 'C_...'.
4) The general adaptive algorithm (under development) is found in spline_aplot.m. See the comments in this file for more details.
5) Read and run test_aplot.m to experiment with the adaptive plotting methods.
6) The goal is to find a good combination of methods in the general algorithm, to automate threshold values and thus provide an algorithm that performs well for all B-splines. Initial results indicate a potential for combining the length-method with the length-ratio method. This combination also allows for some reuse of values, lowering its cost compared to other combinations.
