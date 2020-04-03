close all;
%clear all;
format compact; format long;
addpath(genpath('../../../../../lib/matlab/')) % lib codes
addpath(genpath('../../')) % lib codes
filePath = mfilename('fullpath');
[scriptPath,fileName,fileExt] = fileparts(filePath); cd(scriptPath);

fontSize = 13;

Log10Energy = 0:0.02:4;
Energy = [10.0.^Log10Energy]';
amplitude = 1.e7;

global alpha beta epk
epk = 700;
%alpha = -0.5;
%beta = -2.5;
alpha = -1.9;
beta = -3.5;

PhotonSpectrum = zeros(length(Energy),1);
for i = 1:length(Energy)
    PhotonSpectrum(i) = getBandFunc( Energy(i), epk, alpha, beta );
end
figure; box on;
plot( Energy, PhotonSpectrum, '-' );
set(gca,'xscale','log','fontsize',fontSize);
set(gca,'yscale','log','fontsize',fontSize);
set(gca,'color','none','XMinorTick','on','YMinorTick','on','ZMinorTick','on', 'fontSize', fontSize);

% compute the energy spectrum
figure; box on; hold on;
plot( Energy, Energy.^2 .* PhotonSpectrum, '-' , 'linewidth', 5 );
EnergySpectrum = integrandEnergyUnit( Energy );
plot( Energy, Energy' .* EnergySpectrum, '-' , 'linewidth', 2 , 'color', 'green' );
set(gca,'xscale','log','fontsize',fontSize);
set(gca,'yscale','log','fontsize',fontSize);
set(gca,'color','none','XMinorTick','on','YMinorTick','on','ZMinorTick','on', 'fontSize', fontSize);
hold off;

alphaPlus2 = alpha + 2.
alphaMinusBeta = alpha - beta
ebrk = epk * alphaMinusBeta / alphaPlus2
coef = ebrk^alphaMinusBeta * exp(-alphaMinusBeta)

photonFluence = integral(@integrand,0.1,20000,'reltol',1.e-10)
energyFluence = integral(@integrandEnergyUnit,1,10000,'reltol',1.e-10)


% compute the energy fluence, for a range of epk
LogEpk = log(0.01):0.02:log(500000);
EnergyFluence = zeros(length(LogEpk),1);
EnergyFluenceBatse = zeros(length(LogEpk),1);
EnergyFluenceSwift = zeros(length(LogEpk),1);
counter = 0;
for logepk = LogEpk
    counter = counter + 1;
    epk = exp(logepk);
    EnergyFluence(counter) = integral(@integrandEnergyUnit,1,10000,'reltol',1.e-10) ...
                           / integral(@integrand,1,10000,'reltol',1.e-10);
    EnergyFluenceSwift(counter) = integral(@integrandEnergyUnit,1,10000,'reltol',1.e-10) ...
                                / integral(@integrand,50,150,'reltol',1.e-10);
    EnergyFluenceBatse(counter) = integral(@integrandEnergyUnit,1,10000,'reltol',1.e-10) ...
                                / integral(@integrand,50,300,'reltol',1.e-10);
end

% compute the energy spectrum
figure; box on; hold on;
plot( EnergyFluence, exp(LogEpk), '-' , 'linewidth', 3 );
plot( EnergyFluenceSwift, exp(LogEpk), '-' , 'linewidth', 3 );
plot( EnergyFluenceBatse, exp(LogEpk), '-' , 'linewidth', 3 );
set(gca,'xscale','log','fontsize',fontSize);
set(gca,'yscale','log','fontsize',fontSize);
set(gca,'color','none','XMinorTick','on','YMinorTick','on','ZMinorTick','on', 'fontSize', fontSize);
hold off;
%close all

