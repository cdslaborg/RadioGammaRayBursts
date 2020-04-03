%close all;
%clear all;
format compact; format long;
addpath(genpath('../../../../lib/matlab/')) % lib codes
addpath(genpath('../../')) % lib codes
filePath = mfilename('fullpath');
[scriptPath,fileName,fileExt] = fileparts(filePath); cd(scriptPath);

kfacType = 'OneThird';
if strcmp(kfacType,'OneThird')
    kfac = 0.66;
else
    error('Only kfacType=''OneThird'' is supported as input to readLloydRadioData()')
end
Nicole = readLloydRadioData('OneThird');
cd(scriptPath);

figure; hold on; box on;
h1 = histogram(Nicole.Dark.LogEiso.Val/log(10));
h2 = histogram(Nicole.Bright.LogEiso.Val/log(10));
binWidth = min(h1.BinWidth,h2.BinWidth);
h1.BinWidth = binWidth;
h2.BinWidth = binWidth;
hold off;

% nss stands for NicoleSynSam
Path.input = '../winx64/intel/release/static/serial/bin/out/kfacOneThird/';

nss.Dark = importdata([Path.input,'syntheticSampleB10Dark.csv']);
nss.Bright = importdata([Path.input,'syntheticSampleB10Bright.csv']);
%nss.SynSam = importdata([Path.input,'../../../SynSam.csv']);

% bring everything to log10 scale
%nss.Dark.data = nss.Dark.data / log(10);
%nss.Bright.data = nss.Bright.data / log(10);
RadioType = {'Dark','Bright'};
nss.Dark.ncol = length(nss.Dark.data(1,:));
nss.Bright.ncol = length(nss.Bright.data(1,:));
for icol = 1:nss.Dark.ncol

    figure; hold on; box on;
    hdark = histogram(nss.Dark.data(:,icol)/log(10));
    hbright = histogram(nss.Bright.data(:,icol)/log(10));
    binWidth = (hdark.BinWidth + hbright.BinWidth) / 2.0;
    hdark.BinWidth = binWidth;
    hbright.BinWidth = binWidth;
    xlabel(nss.Dark.textdata{icol});
    legend(RadioType,'location','northeast');
    hold off;
    
end

return

% Avg-Variance plot
nss.Dark.LogZone.Avg = nss.Dark.data(:,9);
nss.Dark.LogDurz.Avg = nss.Dark.data(:,3);
nss.Dark.LogEiso.Avg = nss.Dark.data(:,4);
nss.Dark.LogDurz.Std = nss.Dark.data(:,13);
nss.Dark.LogEiso.Std = nss.Dark.data(:,14);

nss.Bright.LogZone.Avg = nss.Bright.data(:,9);
nss.Bright.LogDurz.Avg = nss.Bright.data(:,3);
nss.Bright.LogEiso.Avg = nss.Bright.data(:,4);
nss.Bright.LogDurz.Std = nss.Bright.data(:,13);
nss.Bright.LogEiso.Std = nss.Bright.data(:,14);

Radio.Type = {'Dark','Bright'};
Radio.count = length(RadioType);
Var.Name = {'LogZone','LogDurz'}; %,'LogEiso'
Var.count = length(Var.Name);
for ivar = 1:Var.count-1
    for jvar = ivar+1:Var.count
        figure; hold on; box on;
        for irtype = 1:Radio.count
            rtype = Radio.Type{irtype};
            plot( exp(nss.(rtype).(Var.Name{ivar}).Avg/log(10)) , exp(nss.(rtype).(Var.Name{jvar}).Avg/log(10)) , '.' );
            xlabel(Var.Name{ivar}(4:end));
            ylabel(Var.Name{jvar}(4:end));
        end
        set(gca,'xscale','log');
        set(gca,'yscale','log');
        legend(Radio.Type)
        hold off;
    end
end

disp('DiffDurz Bright-Dark:')
mean(nss.Bright.LogDurz.Avg/log(10)) - mean(nss.Dark.LogDurz.Avg/log(10))

return


figure; hold on; box on;
scatter ( exp(nss.Dark.LogEiso.Avg) ...
        , exp(nss.Dark.LogDurz.Avg) ...
        , 5*ones(length(nss.Dark.LogEiso.Avg),1) ...
        , nss.Dark.LogEiso.Std ...
        ,'filled' );
set(gca,'xscale','log','yscale','log');
hold off;

figure; hold on; box on;
plot ( exp(nss.Dark.LogEiso.Avg) ...
     , nss.Dark.LogEiso.Std ...
     ,'.' );
set(gca,'xscale','log','yscale','log');
hold off;

figure; hold on; box on;
plot ( nss.Dark.LogEiso.Std ...
     , exp(nss.Dark.LogDurz.Avg) ...
     ,'.' );
set(gca,'xscale','log','yscale','log');
hold off;


nss.Mask.Dim = nss.Dark.LogEiso.Avg<log(4e52) & nss.Dark.LogEiso.Std<1.1;
nss.Mask.Bright = nss.Dark.LogEiso.Avg>log(8e52) & nss.Dark.LogEiso.Std>2.5;
exp(mean(nss.Dark.LogDurz.Avg(nss.Mask.Dim)))
exp(mean(nss.Dark.LogDurz.Avg(nss.Mask.Bright)))




figure;
scatter ( nss.Bright.data(:,9) ...
        , nss.Bright.data(:,7) ...
        , 5*ones(length(nss.Bright.data(:,9)),1) ...
        , nss.Bright.data(:,8) ...
        ,'filled' );
hold off;

