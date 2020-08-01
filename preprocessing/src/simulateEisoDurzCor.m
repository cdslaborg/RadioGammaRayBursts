%function CorSim = simulateEisoDurzCor(zmodel,kfacType)

close all;
clear all;
format compact; format long;
filePath = mfilename('fullpath');
[scriptPath,fileName,fileExt] = fileparts(filePath); cd(scriptPath);
addpath(genpath('../../../../lib/matlab/')) % lib codes

zmodel = 'B10';
kfacType = 'OneThird';
if strcmp(kfacType,'OneThird')
    kfac = 0.66;
else
    error('Only kfacType=''OneThird'' is supported as input to readLloydRadioData()')
end

% import Radio data of Lloyd 2019
Radio = readLloydRadioData(kfacType);

% read the Parameter posterior sample:
ParaPost.rootPath = '../../../20181213_BatseLgrbRedshift/git/cosmicRate/___winx64___/intel/release/static/serial/';
dum = [ParaPost.rootPath,'kfac',kfacType,zmodel,'/bin/out/',zmodel,'_image_1_sample.txt'];
ParaPost.Sample = importdata(dum);
ParaPost.Sample.path = dum;
ParaPost.Sample.count = length(ParaPost.Sample.data(:,1));

% construct the conditional distribution of Durz given Eiso
AvgLogEiso  = ParaPost.Sample.data(:,4);
AvgLogDurz  = ParaPost.Sample.data(:,5);
StdLogEiso  = exp(ParaPost.Sample.data(:,8)); % / log(10.0)
StdLogDurz  = exp(ParaPost.Sample.data(:,9)); % / log(10.0)
RhoEisoDurz = tanh( ParaPost.Sample.data(:,15) );
DurzGivenEiso.Tilt = RhoEisoDurz .* StdLogDurz ./ StdLogEiso;
DurzGivenEiso.Bias = AvgLogDurz - AvgLogEiso .* DurzGivenEiso.Tilt;

RadioType = {'Dark','Bright'};
for thisType = RadioType
    radioType = thisType{1};
    Radio.(radioType).LogEiso = log( Radio.(radioType).Eiso );
    Radio.(radioType).LogT90z = log( Radio.(radioType).T90z );
    Radio.(radioType).LogZPlusOne = log( Radio.(radioType).ZPlusOne );
    DurzGivenEiso.Avg.(radioType) = zeros(Radio.(radioType).count,length(AvgLogEiso));
    DurzGivenEiso.Std.(radioType) = zeros(Radio.(radioType).count,length(AvgLogEiso));
    for isample = 1:Radio.(radioType).count
        DurzGivenEiso.Avg.(radioType)(isample,:) = DurzGivenEiso.Bias + DurzGivenEiso.Tilt * Radio.(radioType).LogEiso(isample);
        DurzGivenEiso.Std.(radioType)(isample,:) = sqrt( 1.0 - RhoEisoDurz.^2 ) .* StdLogDurz; % this is needed in matrix form for vectorized normrnd generation
    end
end

for thisType = RadioType

    radioType = thisType{1};
    disp(['Simulating correlation for',radioType]);

    nsim = 100;
    CorEisoDurz.(radioType).Spearman.Rho  = zeros(ParaPost.Sample.count,nsim);
    CorEisoDurz.(radioType).Spearman.pval = zeros(ParaPost.Sample.count,nsim);
    for isim = 1:nsim
        disp(['sim # ',num2str(isim)]);
        LogDurz = normrnd( DurzGivenEiso.Avg.(radioType) , DurzGivenEiso.Std.(radioType) );
        for iparapost = 1:ParaPost.Sample.count
            [ CorEisoDurz.(radioType).Spearman.Rho(iparapost,isim) ...
            , CorEisoDurz.(radioType).Spearman.Pval(iparapost,isim) ...
            ] = corr( Radio.(radioType).LogEiso , LogDurz(:,iparapost) , 'Type' , 'Spearman' );
        end
    end

end

% plot the resulting
%figure; hold on;
