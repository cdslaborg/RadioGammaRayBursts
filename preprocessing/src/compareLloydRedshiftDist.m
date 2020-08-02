% perform KS-test on the redshift distribution of the radio-loud and radio-quiet LGRBs.
%close all;
%clear all;
format compact; format long;
filePath = mfilename('fullpath');
[scriptPath,fileName,fileExt] = fileparts(filePath); cd(scriptPath);
addpath(genpath('../../../../lib/matlab/')) % lib codes

figExportRequested = 0;
fontSize = 13;
kfacType = 'OneThird';

% import Radio data of Lloyd 2019
Radio = readLloydRadioData(kfacType);

kstest.outPath = '../out/';
kstest.datPath = '../data/';

[kstest.h,kstest.p,kstest.ks2stat] = kstest2( Radio.Bright.Zone.Val, Radio.Dark.Zone.Val )

