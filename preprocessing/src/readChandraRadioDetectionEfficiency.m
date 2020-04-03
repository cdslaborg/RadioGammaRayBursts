readPoonamFluxLimitTable % This gives the Poonam object
close all;
%clear all;
format compact; format long;
filePath = mfilename('fullpath');
[scriptPath,fileName,fileExt] = fileparts(filePath); cd(scriptPath);
addpath(genpath('../../../../lib/matlab/')) % lib codes

transparencyRequested = 1;
figExportRequested = 1;
markerSize = 20;
lineWidth = 2;
fontSize = 13;
if transparencyRequested
    backColor = 'none';
else
    backColor = 'white';
end
 

dataPath = '../data/';
outPath = '../out/';
load([dataPath,'ChandraRadioDetectionEfficiency.mat']);
figure;
plot( ChandraRadioDetectionEfficiency(:,1) ...
    , ChandraRadioDetectionEfficiency(:,2) ...
    , '-.', 'markerSize', markerSize ...
    )
close(gcf);

% get the log-scale PDF
ChandraRadioDetectionEfficiency(:,3) = ChandraRadioDetectionEfficiency(:,2) ...
                                    .* ChandraRadioDetectionEfficiency(:,1);
figure;
plot( ChandraRadioDetectionEfficiency(:,1) ...
    , ChandraRadioDetectionEfficiency(:,3) ...
    , '-.', 'markerSize', markerSize ...
    )
set(gca,'xscale','log','fontsize',fontSize);
close(gcf);

CumSumDetEff.Flux = ChandraRadioDetectionEfficiency(:,1);
CumSumDetEff.LogFlux = log( CumSumDetEff.Flux );
CumSumDetEff.DetEff = cumsum(ChandraRadioDetectionEfficiency(:,2));
CumSumDetEff.DetEffLogScale = cumsum(ChandraRadioDetectionEfficiency(:,3));
figure;
plot( CumSumDetEff.Flux ...
    , CumSumDetEff.DetEffLogScale/CumSumDetEff.DetEffLogScale(end) ...
    , '-.', 'markerSize', markerSize ...
    )
set(gca,'xscale','log','fontsize',fontSize);
close(gcf);

ChandraRadioPeakFlux = importdata([dataPath,'ChandraRadioPeakFlux_clean.xlsx']);
ChandraRadioPeakFlux.count = length(ChandraRadioPeakFlux.data(:,1));
ChandraRadioPeakFlux.Freq = cell(ChandraRadioPeakFlux.count,1);
for isample = 1:ChandraRadioPeakFlux.count
    ChandraRadioPeakFlux.Freq{isample} = sprintf('%0.2f', ChandraRadioPeakFlux.data(isample,1) );
end


%ChandraRadioFluxData = ChandraRadioPeakFlux.data(:,2);
Mask = strcmp(ChandraRadioPeakFlux.Freq,'8.46');
%Mask = logical(ones(length(ChandraRadioPeakFlux.Freq),1));
ChandraRadioFluxData = ChandraRadioPeakFlux.data(Mask,2);
[~,edges] = histcounts(log10(ChandraRadioFluxData));

if figExportRequested, figure('visible','off','Color',backColor), else, figure, end
    hold on; box on; %colormap('cool');

    histogram(ChandraRadioFluxData,10.^edges,'Normalization','probability')
   %plot( CumSumDetEff.Flux ...
   %    , CumSumDetEff.DetEffLogScale/CumSumDetEff.DetEffLogScale(end) ...
     ..., CumSumDetEff.DetEffLogScale*CumSumDetEff.DetEff(end)/CumSumDetEff.DetEffLogScale(end) ...   
	plot( Poonam.SortedRMS3 ...
        , Poonam.SortedRMS3CDF ...
        , '-.', 'markerSize', markerSize ...
        , 'lineWidth', lineWidth ...
        )
    legend  ( { 'Detected Sample' , '3\sigma Non-detection Limit' } ...
                , 'location' , 'northeast' ...
                , 'fontSize' , fontSize ...
                , 'color' , backColor ...
                )
    legend boxoff
    xlim([20,1e5]);
    set(gca,'xscale','log','fontsize',fontSize);
    xlabel('Peak Radio Flux Density at 8.5 GHz: F_{radio} [ \muJy ]', 'Interpreter', 'tex', 'fontSize', fontSize);
    ylabel('Probability Density Function', 'Interpreter', 'tex', 'fontSize', fontSize);
    set(gca,'color',backColor);

if figExportRequested
    fileName = [outPath,'ChandraRadioDetectionEfficiency.png'];
    %export_fig (fileName,'-dpdf -m2');
    export_fig (fileName,'-m4 -transparent');
    %print(gcf, '-dpdf', fileName);
    hold off; close(gcf);
else
    hold off;
end


