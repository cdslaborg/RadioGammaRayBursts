close all;
clear all;
format compact; format long;
filePath = mfilename('fullpath');
[scriptPath,fileName,fileExt] = fileparts(filePath); scriptPath = string(scriptPath); cd(scriptPath);
addpath(genpath("../../../../../libmatlab/"),"-begin") % lib codes
addpath(genpath("../data/"))

%readChandraData;
ChandraTableOneClean = importdata("ChandraTableOne_clean.xlsx");

t90Data = ChandraTableOneClean.data(:,7);
radioEmission = ChandraTableOneClean.data(:,16);
t90Data = t90Data(~isnan(radioEmission));
radioEmission = radioEmission(~isnan(radioEmission));
radioEmission = radioEmission(~isnan(t90Data));
t90Data = t90Data(~isnan(t90Data));

scatter(t90Data,radioEmission);
xlabel('T90')
ylabel('E_{iso}J')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
sprintf("correlation between T90 and radio emission is; %0.5f",corr(t90Data,radioEmission,"type","spearman"))