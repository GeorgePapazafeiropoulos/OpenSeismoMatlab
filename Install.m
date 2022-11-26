clear all
close all
clc
format long

S = mfilename('fullpath');
f = filesep;
ind=strfind(S,f);
Fname=S(ind(end)+1:end);

S1=S(1:ind(end)-1);
addpath(genpath(S1))
cd(S1)
dbstop if error






