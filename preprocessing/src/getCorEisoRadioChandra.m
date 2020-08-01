close all;
clear all;
format compact; format long;
filePath = mfilename('fullpath');
[scriptPath,fileName,fileExt] = fileparts(filePath); cd(scriptPath);
addpath(genpath('../../../../lib/matlab/')) % lib codes

cd([scriptPath,'/../']);
dataPath = [pwd(),'\data\'];
cd(dataPath);
filePath = [dataPath,'EradEiso.png'];
%grabit();

cd(scriptPath);
figExportRequested = 0;
markerSize = 20;
fontSize = 13;
kfacType = 'OneThird';
bivarFigExportRequested = 0;

path.root = "..";
path.out = fullfile(path.root,"out");
if ~exist(path.out,'dir'), mkdir(path.out), end
path.in = fullfile(path.root,"data");
load(fullfile(path.in,"ChandraLog10EradLog10Eiso.mat")); % loads ChandraLog10EradLog10Eiso array

% plot

rlum = 10.^ChandraLog10EradLog10Eiso(:,1); % radio luminosity
eiso = 10.^ChandraLog10EradLog10Eiso(:,2) * 1.e51;
logEiso = log(eiso);
logRlum = log(rlum);
XRFPointIndex=4;

if figExportRequested
    figure('Color','none');
else
    figure("color","white");
end

hold on; box on;

    plot( rlum ...
        , eiso ...
        , "." ...
        , "markerSize", 20 ...
        , "color", "magenta" ...
        );
    
    plot( rlum(XRFPointIndex) ... %plot XRF point
        , eiso(XRFPointIndex) ...
        , "s" ...
        , "markerSize", 10 ...
        , "color", "green" ...
        ,"MarkerFaceColor","green" ...
        );
    
    set(gca,'xscale','log','yscale','log');
    xlabel("Peak Radio Luminosity at 8.5 GHz [ ergs / s / Hz ]", "fontSize", fontSize);
    ylabel("Total Isotropic Gamma-Ray Emission [ ergs ]", "fontSize", fontSize);
    
if figExportRequested
    fileName = [dataPath,'EradEiso','.png'];
    export_fig (fileName,'-m4 -transparent');
    hold off; close(gcf);
else
    hold off;
end
    

sprintf('correlation with XRF: %d',corr(logRlum,logEiso,"type","spearman")) % correlation with XRF

sprintf('correlation without XRF: %d',corr(logRlum([1:XRFPointIndex-1 XRFPointIndex+1:end]),...
logEiso([1:XRFPointIndex-1 XRFPointIndex+1:end]),"type","spearman"))


