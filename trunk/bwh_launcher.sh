#!/bin/sh

matlab_exec=/Applications/MATLAB_R2013b.app/Contents/MacOS/StartMATLAB
X="bwh_registration(${1},${2},${3})"
echo Launching command ${X}
echo ${X} > matlab_command_${1}_${2}_${3}.m
cat matlab_command_${1}_${2}_${3}.m
${matlab_exec} -nodisplay -nosplash <  matlab_command_${1}_${2}_${3}.m
rm  matlab_command_${1}_${2}_${3}.m
