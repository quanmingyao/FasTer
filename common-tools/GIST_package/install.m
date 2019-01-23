clear
clc


currentpath = cd;

addpath(genpath([currentpath,'/GIST']));

cd ./GIST

mex proximalRegC.cpp
mex funRegC.cpp

cd ..

