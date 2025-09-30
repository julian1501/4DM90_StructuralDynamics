close all;
clearvars;

% import data from data folder
sch1dof = importdata("data\SCH1DOF.MAT");
dt = sch1dof.dt; % seconds
fs = 1/dt; % Hz - sampling frequency
Nsamples = sch1dof.N; % number of samples

input = sch1dof.excit;
output = sch1dof.output;

nwindows = 10;
windowsize = Nsamples/nwindows;
overlap = windowsize/2;

[tf, f] = tfestimate(input,output,windowsize,overlap,Nsamples,fs);

loglog(f,abs(tf))