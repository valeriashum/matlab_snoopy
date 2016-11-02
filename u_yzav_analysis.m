clear all;
close all;
rep1='/data/novadisk/vs391/hydro/u_iii/run_1/result1/';
files=[rep1,'uz_yzav.dat'];
full_file=importdata(files);
timevar=full_file.data;