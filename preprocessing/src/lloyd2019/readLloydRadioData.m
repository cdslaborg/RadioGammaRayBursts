function Radio = readLloydRadioData(kfacType)
% read L19 radio-bright and radio-quiet data from the corresponding in Excel files.
% input:
%   kfacType: a string with only two possible values: "OneThird", "none"
% output:
%   Radio: a structure with two major components (Dark,Bright) that store the properties of L19 sample.

    %close all;
    %clear all;
    format compact; format long;
    filePath = mfilename('fullpath');
    [scriptPath,fileName,fileExt] = fileparts(filePath); cd(scriptPath);
    addpath(genpath('../../../../lib/matlab/')) % lib codes

    if nargin<1
        warning('no kfacType was provided as input to readLloydRadioData(). Assuming kfacType=''OneThird''')
        kfac = 0.66;
    elseif strcmpi(kfacType,"OneThird")
        kfac = 0.66;
    elseif strcmpi(kfacType,"none")
        kfac = 1;
    else
        error('Only kfacType="OneThird" and kfacType="none" are supported as input to readLloydRadioData()')
    end


    Radio.Type.Name = {'Dark','Bright'};
    Radio.Type.count = length(Radio.Type.Name);


    fontSize = 13;
    Radio.MarkerSize.default = 25;
    Radio.MarkerSize.overlay = 13;


    Radio.Path.src = scriptPath;
    Radio.Path.input = '../data/';
    Radio.Dark = importdata([Radio.Path.input,'radioDark.xlsx']);
    Radio.Dark.Zone.Val = 1.0 + Radio.Dark.data(:,4);
    Radio.Dark.Eiso.Val = Radio.Dark.data(:,2) * 10.0^52;
    kcorrection = Radio.Dark.Zone.Val;
    if kfac~=1
        kcorrection = kcorrection .^ kfac;
    end
    Radio.Dark.Durz.Val = Radio.Dark.data(:,3) ./ kcorrection;
    Radio.Dark.count = length(Radio.Dark.Durz.Val);
    Radio.Dark.LogZone.Val = log( Radio.Dark.Zone.Val );
    Radio.Dark.LogEiso.Val = log( Radio.Dark.Eiso.Val );
    Radio.Dark.LogDurz.Val = log( Radio.Dark.Durz.Val );


    Radio.Bright = importdata([Radio.Path.input,'radioBright.xlsx']);
    Radio.Bright.Zone.Val = 1.0 + Radio.Bright.data(:,3);
    Radio.Bright.Eiso.Val = Radio.Bright.data(:,1) * 10.0^52;
    kcorrection = Radio.Bright.Zone.Val; if kfac~=1; kcorrection = kcorrection .^ kfac; end
    Radio.Bright.Durz.Val = Radio.Bright.data(:,2) ./ kcorrection;
    Radio.Bright.count = length(Radio.Bright.Durz.Val);
    Radio.Bright.LogZone.Val = log( Radio.Bright.Zone.Val );
    Radio.Bright.LogEiso.Val = log( Radio.Bright.Eiso.Val );
    Radio.Bright.LogDurz.Val = log( Radio.Bright.Durz.Val );


    %disp('Eiso-Durz');
    %figure; hold on;
    %plot(Radio.Dark.Eiso.Val,Radio.Dark.Durz.Val,'.','markerSize',Radio.MarkerSize.default);
    %plot(Radio.Bright.Eiso.Val,Radio.Bright.Durz.Val,'.','markerSize',Radio.MarkerSize.default);
    %set(gca,'xscale','log','fontsize',fontSize);
    %set(gca,'yscale','log','fontsize',fontSize);
    %xlabel('E_{iso}', 'Interpreter', 'tex', 'fontSize', fontSize);
    %ylabel('T_{90z}', 'Interpreter', 'tex', 'fontSize', fontSize);
    %hold off;


    %disp('z-Durz');
    %figure; hold on;
    %plot(Radio.Dark.Zone.Val,Radio.Dark.Durz.Val,'.','markerSize',Radio.MarkerSize.default);
    %plot(Radio.Bright.Zone.Val,Radio.Bright.Durz.Val,'.','markerSize',Radio.MarkerSize.default);
    %set(gca,'xscale','log','fontsize',fontSize);
    %set(gca,'yscale','log','fontsize',fontSize);
    %xlabel('(1+z)', 'Interpreter', 'tex', 'fontSize', fontSize);
    %ylabel('T_{90z}', 'Interpreter', 'tex', 'fontSize', fontSize);
    %hold off;
    %corr(log(Radio.Dark.Durz),log(Radio.Dark.Zone),'type','spearman')
    %corr(log(Radio.Bright.Durz),log(Radio.Bright.Zone),'type','spearman')


    % all data
    Radio.All.LogZone.Val = vertcat(Radio.Dark.LogZone.Val,Radio.Bright.LogZone.Val);
    Radio.All.LogEiso.Val = vertcat(Radio.Dark.LogEiso.Val,Radio.Bright.LogEiso.Val);
    Radio.All.LogDurz.Val = vertcat(Radio.Dark.LogDurz.Val,Radio.Bright.LogDurz.Val);
    Radio.All.Zone.Val = vertcat(Radio.Dark.Zone.Val,Radio.Bright.Zone.Val);
    Radio.All.Eiso.Val = vertcat(Radio.Dark.Eiso.Val,Radio.Bright.Eiso.Val);
    Radio.All.Durz.Val = vertcat(Radio.Dark.Durz.Val,Radio.Bright.Durz.Val);
    Radio.All.count = length(Radio.Bright.Durz.Val);


    RadioType = {'Dark','Bright','All'};
    VarType = {'LogZone','LogEiso','LogDurz'};
    CorType = {'ZoneEiso','ZoneDurz','EisoDurz'};
    VarTypeScalar = {'logZone','logEiso','logDurz'};
    for irtype = 1:length(RadioType)

        rtype = RadioType{irtype};

        %compute correlations

        for icor = 1:length(CorType)
            cor = CorType{icor};
            logxval = Radio.(rtype).(['Log',cor(1:4)]).Val;
            logyval = Radio.(rtype).(['Log',cor(5:8)]).Val;
            [ Radio.(rtype).Cor.Spearman.(cor).coef ...
            , Radio.(rtype).Cor.Spearman.(cor).pval ...
            ] = corr( logxval ... 
                    , logyval ...
                    , 'type', 'spearman' ...
                    );
        end

        % compute statistics

        for ivar = 1:length(VarType)
            var = VarType{ivar};
            Radio.(rtype).(var).avg = mean(Radio.(rtype).(var).Val);
            Radio.(rtype).(var).std = std(Radio.(rtype).(var).Val);
            Radio.(rtype).(var(4:7)).avg = exp( mean(Radio.(rtype).(var).Val) );
            Radio.(rtype).(['Log10',var(4:7)]).avg = mean(Radio.(rtype).(var).Val)/log(10);
            Radio.(rtype).(['Log10',var(4:7)]).std = std(Radio.(rtype).(var).Val)/log(10);
        end

    end

    %disp('All Eiso-Durz');
    %figure; hold on;
    %plot(Radio.All.Eiso,Radio.All.Durz,'.','markerSize',Radio.MarkerSize.default);
    %xlabel('E_{iso}', 'Interpreter', 'tex', 'fontSize', fontSize);
    %ylabel('T_{90z}', 'Interpreter', 'tex', 'fontSize', fontSize);
    %set(gca,'xscale','log','fontsize',fontSize);
    %set(gca,'yscale','log','fontsize',fontSize);
    %hold off;
    %corr(log(Radio.All.Eiso),log(Radio.All.Durz),'type','spearman')
    %
    %disp('All z-Durz');
    %figure; hold on;
    %plot(Radio.All.Zone,Radio.All.Durz,'.','markerSize',Radio.MarkerSize.default);
    %xlabel('(1+z)', 'Interpreter', 'tex', 'fontSize', fontSize);
    %ylabel('T_{90z}', 'Interpreter', 'tex', 'fontSize', fontSize);
    %set(gca,'xscale','log','fontsize',fontSize);
    %set(gca,'yscale','log','fontsize',fontSize);
    %hold off;
    %corr(log(Radio.All.Zone),log(Radio.All.Durz),'type','spearman')

end