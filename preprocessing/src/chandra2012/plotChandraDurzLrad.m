close all;
clear all;
format compact; format long;
filePath = mfilename('fullpath');
[scriptPath,fileName,fileExt] = fileparts(filePath); scriptPath = string(scriptPath); cd(scriptPath);
addpath(genpath("../../../../../libmatlab/"),"-begin") % lib codes
addpath(genpath("../data/"))

inPath = getFullPath( fullfile(scriptPath,"..","..","data","chandra2012\") );
outPath = getFullPath( fullfile(scriptPath,"..","..","out","chandra2012\"));
%readChandraData;
ChandraTableOneClean = importdata(inPath+"ChandraTableOne_clean.xlsx");
exportDurzLradPlot = 1;


t90Data = ChandraTableOneClean.data(:,7);
redshift = ChandraTableOneClean.data(:,8);
radioEmission = ChandraTableOneClean.data(:,16);
mask = ~( isnan(radioEmission) | isnan(t90Data) | isnan(redshift));
Zone = redshift(mask) + 1;
t90Data = t90Data(mask);
radioEmission = radioEmission(mask);
Durz = t90Data ./ Zone.^0.66;

scatter(Durz,radioEmission);
xlabel('T90')
ylabel('E_{iso}J')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
sprintf("correlation between T90 and radio emission is; %0.5f",corr(Durz,radioEmission,"type","spearman"))
if exportDurzLradPlot
    filePath = getFullPath( fullfile(outPath,"DurzLradPlot.png"));
    export_fig(filePath,'-m4 -transparent')
end