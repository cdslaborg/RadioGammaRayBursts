close all;
clear all;
format compact; format long;
path.projects.dir = getFullPath("../../../../../","-lean"); % getFullPath is called from libmatlab
addpath(genpath(fullfile(path.projects.dir,"libmatlab")),"-begin") % libmatlab codes
filePath = mfilename("fullpath");
[scriptPath,fileName,fileExt] = fileparts(filePath); cd(scriptPath);

path.git.dir = string(getFullPath(fullfile(path.projects.dir,"20190419_RadioGRBs","git")));
addpath(genpath(fullfile(path.git.dir)),"-begin") % libmatlab codes

% read Nicole data

kfacType = "OneThird";
if strcmp(kfacType,"OneThird")
    kfac = 0.66;
else
    error("Only kfacType=""OneThird"" is supported as input to readLloydRadioData()");
end
Nicole = readLloydRadioData("OneThird");
cd(scriptPath);

% setup path

path.git.synsam.dir = fullfile(path.git.dir,"SyntheticSample");
path.git.synsam.in.dir = fullfile(path.git.synsam.dir,"postproc");
path.git.synsam.in.batse = fullfile(path.git.synsam.in.dir,"batse_1366_lgrb_pbol_epk_sbol(0.001,20000).txt");
path.git.synsam.in.ghirlanda08 = fullfile(path.git.synsam.in.dir,"AmatiRelationGhirlanda2008.txt");
path.git.synsam.postproc.dir = fullfile(path.git.synsam.dir,"postproc");
path.git.synsam.postproc.src.dir = fullfile(path.git.synsam.postproc.dir,"src");
path.git.synsam.postproc.out.dir = fullfile(path.git.synsam.postproc.dir,"out");
path.git.synsam.out.dir = fullfile(path.git.synsam.dir,"winx64","intel","release","static","serial","bin","out",kfacType);
path.git.synsam.out.synsam = fullfile(path.git.synsam.out.dir,"syntheticSampleB10Dark.csv");
path.git.synsam.out.nssDarkB10 = fullfile(path.git.synsam.out.dir,"syntheticSampleB10Dark.csv");
path.git.synsam.out.nssBrightB10 = fullfile(path.git.synsam.out.dir,"syntheticSampleB10_B10detectionCriterion.csv");

figure; hold on; box on;
h1 = histogram(Nicole.Dark.LogEiso.Val/log(10));
h2 = histogram(Nicole.Bright.LogEiso.Val/log(10));
binWidth = min(h1.BinWidth,h2.BinWidth);
h1.BinWidth = binWidth;
h2.BinWidth = binWidth;
hold off;

% read Nichol-based synthetic sample. nss stands for NicoleSynSam

nss.Dark = importdata(path.git.synsam.out.nssDarkB10);
nss.Dark = importdata(path.git.synsam.out.nssBrightB10);
%nss.SynSam = importdata(path.git.synsam.out.synsam);

% bring everything to log10 scale
%nss.Dark.data = nss.Dark.data / log(10);
%nss.Bright.data = nss.Bright.data / log(10);

RadioType = ["Dark","Bright"];
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
    legend(RadioType,"location","northeast");
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

Radio.Type = {"Dark","Bright"};
Radio.count = length(RadioType);
Var.Name = {"LogZone","LogDurz"}; %,"LogEiso"
Var.count = length(Var.Name);
for ivar = 1:Var.count-1
    for jvar = ivar+1:Var.count
        figure; hold on; box on;
        for irtype = 1:Radio.count
            rtype = Radio.Type{irtype};
            plot( exp(nss.(rtype).(Var.Name{ivar}).Avg/log(10)) , exp(nss.(rtype).(Var.Name{jvar}).Avg/log(10)) , "." );
            xlabel(Var.Name{ivar}(4:end));
            ylabel(Var.Name{jvar}(4:end));
        end
        set(gca,"xscale","log");
        set(gca,"yscale","log");
        legend(Radio.Type)
        hold off;
    end
end

disp("DiffDurz Bright-Dark:")
mean(nss.Bright.LogDurz.Avg/log(10)) - mean(nss.Dark.LogDurz.Avg/log(10))

return


figure; hold on; box on;
scatter ( exp(nss.Dark.LogEiso.Avg) ...
        , exp(nss.Dark.LogDurz.Avg) ...
        , 5*ones(length(nss.Dark.LogEiso.Avg),1) ...
        , nss.Dark.LogEiso.Std ...
        ,"filled" );
set(gca,"xscale","log","yscale","log");
hold off;

figure; hold on; box on;
plot ( exp(nss.Dark.LogEiso.Avg) ...
     , nss.Dark.LogEiso.Std ...
     ,"." );
set(gca,"xscale","log","yscale","log");
hold off;

figure; hold on; box on;
plot ( nss.Dark.LogEiso.Std ...
     , exp(nss.Dark.LogDurz.Avg) ...
     ,"." );
set(gca,"xscale","log","yscale","log");
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
        ,"filled" );
hold off;

