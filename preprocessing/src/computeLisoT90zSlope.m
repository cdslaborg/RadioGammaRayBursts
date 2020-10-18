close all;
clear all;
format compact; format long;
filePath = mfilename('fullpath');
[scriptPath,fileName,fileExt] = fileparts(filePath); cd(scriptPath);

addpath(genpath('D:\Dropbox\Projects\20180101_ParaMonte\git\build\libparamonte_MATLAB'),"-begin") % add ParaMonte MATLAB
addpath(genpath('../../../../libmatlab/'),"-begin") % added by josh
addpath(genpath('./lloyd2019/')) % added by josh
addpath(genpath('./chandra2012/')) % added by josh

% plot histogram of log10(T90/Tr45)

kfacType = 'OneThird';
ZoneRequested = 1;
bivarPlotsRequested = 1;
histFigExportRequested = 0;
bivarFigExportRequested = 1;
fontSize = 13;
myYellow = [1,1,0];
myGrey = [0 0 0]; % determines the gray color strength, higher means whiter
%myGrey = [0.6 0.6 0.6]; % determines the gray color strength, higher means whiter
%SynSam.Path.root = '../../../../20181213_BatseLgrbRedshift/git/SyntheticSample/';
SynSam.Path.root = '../../../../../20181213_BatseLgrbRedshift/git/___SyntheticSample___/'; % added by Josh
SynSam.Path.output = [SynSam.Path.root,'winx64/intel/release/static/serial/bin/out/kfac',kfacType,'/'];
SynSam.Path.input = [SynSam.Path.root,'in/'];

% import MCMC data

chainFilePath = "D:\Dropbox\Projects\20181213_BatseLgrbRedshift\git\cosmicRate\___winx64Old___\intel\19.0.2.190\release\static\serial\kfacOneThirdB10\bin";
pm = paramonte();
pmpd = pm.ParaDRAM();
chain = pmpd.readMarkovChain(chainFilePath); chain = chain{1};
%chain.plot.grid.make("columns",[16:22]);
chain.df.Properties.VariableNames

slopeLisoDurzRelation = tanh(chain.df.atanhCorLisoT90z) .* exp(chain.df.stdLogT90z) ./ exp(chain.df.stdLogLiso);
slopeEisoDurzRelation = tanh(chain.df.atanhCorEisoT90z) .* exp(chain.df.stdLogT90z) ./ exp(chain.df.stdLogEiso);
figure('visible','on','Color','white'); hold on; box on;
    hLiso = histogram(slopeLisoDurzRelation); hLiso.EdgeColor = "none"; 
    hEiso = histogram(slopeEisoDurzRelation); hEiso.EdgeColor = "none"; 
    legend(["L_{iso}-T_{90z}", "E_{iso}-T_{90z}"],"fontsize",12,"interpreter","tex");
    xlabel("\alpha in T_{90z} \propto Energetics^{\alpha} Relationship","fontsize",12,"interpreter","tex");
    ylabel("MCMC Count","fontsize",12,"interpreter","tex");
    export_fig ("../out/lisoEisoDurzSlopeDist.png",'-m4 -transparent');
hold off;


slopeLisoEpkzRelation = tanh(chain.df.atanhCorLisoEpkz) .* exp(chain.df.stdLogEpkz) ./ exp(chain.df.stdLogLiso);
slopeEisoEpkzRelation = tanh(chain.df.atanhCorEpkzEiso) .* exp(chain.df.stdLogEpkz) ./ exp(chain.df.stdLogEiso);

figure; hold on; box on;
histogram(slopeLisoEpkzRelation);
histogram(slopeEisoEpkzRelation);
hold off;





