clear; clc;

% mex -v OPTIMFLAGS="/Ox" MakepartM1_cblas.cpp D:\Programs\Matlab2016a\extern\lib\win64\microsoft\*.lib;


mex -v OPTIMFLAGS="/Ox" MakepartM1_c.cpp;
mex -v OPTIMFLAGS="/Ox" MakepartM2_c.cpp;
mex -v OPTIMFLAGS="/Ox" MakepartM3_c.cpp;