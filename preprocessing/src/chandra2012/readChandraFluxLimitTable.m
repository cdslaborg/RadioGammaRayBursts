close all;
clear all;
format compact; format long;
addpath(genpath('../../../../../../lib/matlab/')) % lib codes
addpath(genpath('../../../../../')) % git codes
filePath = mfilename('fullpath');
[scriptPath,fileName,fileExt] = fileparts(filePath); cd(scriptPath);

Path.input = '../../data/chandra2012/';
Path.output = '../../out/chandra2012/';
Poonam = importdata([Path.input,'POONAM-MASTER-RADIO-TABLE.xlsx']);
Poonam.DaysSinceBurst = Poonam.data(:,6);
Poonam.RadioFreq = Poonam.data(:,7);
Poonam.RadioFluxDensity = Poonam.data(:,8);
Poonam.RadioFluxDensityRMS = Poonam.data(:,9);
Poonam.RadioFluxDensityRMS3 = 3*Poonam.RadioFluxDensityRMS;
Poonam.Mask.Freq8.All = (Poonam.RadioFreq<9.0 & Poonam.RadioFreq>8.0);
Poonam.Mask.Freq8.Detection.All = Poonam.Mask.Freq8.All & Poonam.RadioFluxDensity>0;
Poonam.Mask.Freq8.NonDetect.All = Poonam.Mask.Freq8.All & Poonam.RadioFluxDensity<0;
Poonam.Mask.Freq8.Detection.Day510 = Poonam.Mask.Freq8.Detection.All & (Poonam.DaysSinceBurst>=5 & Poonam.DaysSinceBurst<=10);
Poonam.Mask.Freq8.NonDetect.Day510 = Poonam.Mask.Freq8.NonDetect.All & (Poonam.DaysSinceBurst>=5 & Poonam.DaysSinceBurst<=10);

% plot peak radio flux at 8.46 GHz

Poonam.RadioFreqUniq(:).value = unique(Poonam.RadioFreq);   % unique values
Poonam.RadioFreqUniq(:).count = countmember( Poonam.RadioFreqUniq(:).value , Poonam.RadioFreq ); % unique values counts in the sample

Poonam.SortedRMS3 = sort( Poonam.RadioFluxDensityRMS3(Poonam.Mask.Freq8.NonDetect.Day510) );
[Poonam.SortedRMS3CDF,Poonam.SortedRMS3] = ecdf( Poonam.SortedRMS3 );
plot(Poonam.SortedRMS3, Poonam.SortedRMS3CDF)
set(gca,'xscale','log');
xlim([20,1e5]);

return

% detected events flux
figure; hold on; box on;
plot( Poonam.DaysSinceBurst(Poonam.Mask.Freq8.Detection.All) ...
    , Poonam.RadioFluxDensity(Poonam.Mask.Freq8.Detection.All) ...
    , '.' ...
    , 'markersize' , 10 ...
    )
%xlim([0.1 1000]);
%ylim([5 2000]);
set(gca,'xscale','log','yscale','log');
hold off;

% Nondetection upper limits
figure; hold on; box on;
plot( Poonam.DaysSinceBurst(Poonam.Mask.Freq8.NonDetect.All) ...
    , Poonam.RadioFluxDensityRMS3(Poonam.Mask.Freq8.NonDetect.All) ...
    , '.' ...
    , 'markersize' , 10 ...
    )
xlim([0.1 1000]);
ylim([5 2000]);
set(gca,'xscale','log','yscale','log');
hold off;

figure; hold on; box on;
%histogram( log10(Poonam.RadioFluxDensityRMS(Poonam.Mask.Freq8.NonDetect.All)) )
histogram( log10(Poonam.RadioFluxDensityRMS3(Poonam.Mask.Freq8.NonDetect.All)) )
histogram( log10(Poonam.RadioFluxDensity(Poonam.Mask.Freq8.Detection.All)) )
hold off;

% linear space
figure; hold on; box on;
%histogram( log10(Poonam.RadioFluxDensityRMS(Poonam.Mask.Freq8.NonDetect.All)) )
h1 = histogram( Poonam.RadioFluxDensityRMS3(Poonam.Mask.Freq8.NonDetect.Day510) );
h2 = histogram( Poonam.RadioFluxDensity(Poonam.Mask.Freq8.Detection.Day510) );
h1.BinWidth = 50;
h2.BinWidth = 50;
xlim([0 1200]);
hold off;


