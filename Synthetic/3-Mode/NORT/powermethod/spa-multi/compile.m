clear; clc;

% mode-1
mex -v OPTIMFLAGS="/Ox" subSpa1Mu_c.cpp;
% mex -v OPTIMFLAGS="/Ox" subSpa1Mu_cblas.cpp  D:\Programs\Matlab2016a\extern\lib\win64\microsoft\*.lib
mex -v OPTIMFLAGS="/Ox" subSpa1Mv_c.cpp;

% mode-2
mex -v OPTIMFLAGS="/Ox" subSpa2Mu_c.cpp;
mex -v OPTIMFLAGS="/Ox" subSpa2Mv_c.cpp;

% mode-3
mex -v OPTIMFLAGS="/Ox" subSpa3Mu_c.cpp;
mex -v OPTIMFLAGS="/Ox" subSpa3Mv_c.cpp;