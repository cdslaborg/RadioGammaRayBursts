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
outPath = '../out/';
if ~exist(outPath,'dir'), mkdir(outPath), end
figExportRequested = 0;
markerSize = 20;
fontSize = 13;
kfacType = 'OneThird';
bivarFigExportRequested = 1;

if strcmp(kfacType,'OneThird')
    kfacVal = 0.66;
else
    error(['The input kfacType (',kfacType,') is not supported.']);
end

% plot Chandra (2012) data in Figure 20
load([dataPath,'ChandraLog10EradLog10Eiso.mat']);
Chandra.LradEiso.Cor.Pearson.val = corr(ChandraLog10EradLog10Eiso(:,1),ChandraLog10EradLog10Eiso(:,2),'type','pearson');
Chandra.LradEiso.Cor.Kendall.val = corr(ChandraLog10EradLog10Eiso(:,1),ChandraLog10EradLog10Eiso(:,2),'type','kendall');
Chandra.LradEiso.Cor.Spearman.val = corr(ChandraLog10EradLog10Eiso(:,1),ChandraLog10EradLog10Eiso(:,2),'type','Spearman');
Chandra.LradEiso.Cor.Pearson.str = num2str( sprintf('%0.2f', Chandra.LradEiso.Cor.Pearson.val )   ) ;
Chandra.LradEiso.Cor.Kendall.str = num2str( sprintf('%0.2f', Chandra.LradEiso.Cor.Kendall.val ) ) ;
Chandra.LradEiso.Cor.Spearman.str = num2str( sprintf('%0.2f', Chandra.LradEiso.Cor.Spearman.val ) ) ;

if figExportRequested, close all, figure('visible','off','Color','none'), else, figure, end
    hold on; box on; %colormap('cool');
    plot(10.0.^ChandraLog10EradLog10Eiso(:,1),1e51*10.0.^ChandraLog10EradLog10Eiso(:,2),'.','markerSize',markerSize);

    set(gca,'xscale','log','fontsize',fontSize,'XDir','normal');
    set(gca,'yscale','log','fontsize',fontSize,'YDir','normal');

    %xlim(  );
    ylim( [1e48,1e55] );
    xlabel('Peak Radio Luminosity at 8.5 GHz: L_{rad} [ erg / s / Hz ]', 'Interpreter', 'tex', 'fontSize', fontSize);
    ylabel('Isotropic Gamma-Ray Emission: E_{iso} [ erg ]        ', 'Interpreter', 'tex', 'fontSize', fontSize);
    text( 0.05, 0.90, ['r = ',Chandra.LradEiso.Cor.Pearson.str], 'Units', 'normalized', 'HorizontalAlignment', 'left','Interpreter','tex', 'fontSize' , fontSize );
    text( 0.05, 0.80, ['\rho = ',Chandra.LradEiso.Cor.Spearman.str], 'Units', 'normalized', 'HorizontalAlignment', 'left','Interpreter','tex', 'fontSize' , fontSize );
    text( 0.05, 0.70, ['\tau = ',Chandra.LradEiso.Cor.Kendall.str], 'Units', 'normalized', 'HorizontalAlignment', 'left','Interpreter','tex', 'fontSize' , fontSize );
    set(gca, 'color', 'none', 'fontsize', fontSize);

if figExportRequested
    fileName = [outPath,'ChandraLog10EradLog10Eiso.png'];
    export_fig (fileName,'-m4 -transparent');
    hold off; close(gcf);
else
    hold off;
end


% Chandra Sbol-SoptR plot
ChandraTableOneKnownEisoEoptZ = importdata([dataPath,'ChandraTableOneKnownEisoEoptZ_clean.xlsx']);
Chandra.FradSgam.Cor.Pearson.val = corr(ChandraTableOneKnownEisoEoptZ.data(:,4),ChandraTableOneKnownEisoEoptZ.data(:,6),'type','pearson');
Chandra.FradSgam.Cor.Kendall.val = corr(ChandraTableOneKnownEisoEoptZ.data(:,4),ChandraTableOneKnownEisoEoptZ.data(:,6),'type','kendall');
Chandra.FradSgam.Cor.Spearman.val = corr(ChandraTableOneKnownEisoEoptZ.data(:,4),ChandraTableOneKnownEisoEoptZ.data(:,6),'type','Spearman');
Chandra.FradSgam.Cor.Pearson.str = num2str( sprintf('%0.2f', Chandra.FradSgam.Cor.Pearson.val )   ) ;
Chandra.FradSgam.Cor.Kendall.str = num2str( sprintf('%0.2f', Chandra.FradSgam.Cor.Kendall.val ) ) ;
Chandra.FradSgam.Cor.Spearman.str = num2str( sprintf('%0.2f', Chandra.FradSgam.Cor.Spearman.val ) ) ;

if figExportRequested, close all, figure('visible','off','Color','none'), else, figure, end
    hold on; box on; %colormap('cool');
    plot(ChandraTableOneKnownEisoEoptZ.data(:,4),ChandraTableOneKnownEisoEoptZ.data(:,6),'.','markerSize',20)

    set(gca,'xscale','log','fontsize',fontSize,'XDir','normal');
    set(gca,'yscale','log','fontsize',fontSize,'YDir','normal');

    %xlim(  );
    %ylim( [1e48,1e55] );
    xlabel('Gamma-Ray Fluence: S_{15-150} [ erg ]', 'Interpreter', 'tex', 'fontSize', fontSize);
    ylabel('Optical Flux Density at 11 hr: F^{11h}_{R} [ \mu Jy ]', 'Interpreter', 'tex', 'fontSize', fontSize);
    text( 0.05, 0.90, ['r = ',Chandra.FradSgam.Cor.Pearson.str], 'Units', 'normalized', 'HorizontalAlignment', 'left','Interpreter','tex', 'fontSize' , fontSize );
    text( 0.05, 0.80, ['\rho = ',Chandra.FradSgam.Cor.Spearman.str], 'Units', 'normalized', 'HorizontalAlignment', 'left','Interpreter','tex', 'fontSize' , fontSize );
    text( 0.05, 0.70, ['\tau = ',Chandra.FradSgam.Cor.Kendall.str], 'Units', 'normalized', 'HorizontalAlignment', 'left','Interpreter','tex', 'fontSize' , fontSize );
    set(gca, 'color', 'none', 'fontsize', fontSize);

if figExportRequested
    fileName = [outPath,'ChandraFoptSgam.png'];
    export_fig (fileName,'-m4 -transparent');
    hold off; close(gcf);
else
    hold off;
end


% match Chandra Radio data with Table1 data
ChandraTableOne = importdata([dataPath,'ChandraTableOne_clean.xlsx']);
ChandraRadioPeakFlux = importdata([dataPath,'ChandraRadioPeakFlux_clean.xlsx']);

CommonTrig = mapTriggers(ChandraTableOne.textdata(2:end,1),ChandraRadioPeakFlux.textdata(2:end,1),ChandraRadioPeakFlux.data(:,1));
DataX = ChandraRadioPeakFlux.data(CommonTrig.IndxRadio,2);
DataY = ChandraTableOne.data(CommonTrig.IndxTable1,9);

Chandra.FradSgam.Cor.Pearson.val = corr(log(DataX),log(DataY),'type','pearson', 'rows','complete');
Chandra.FradSgam.Cor.Kendall.val = corr(log(DataX),log(DataY),'type','kendall', 'rows','complete');
Chandra.FradSgam.Cor.Spearman.val = corr(log(DataX),log(DataY),'type','Spearman', 'rows','complete');
Chandra.FradSgam.Cor.Pearson.str = num2str( sprintf('%0.2f', Chandra.FradSgam.Cor.Pearson.val )   ) ;
Chandra.FradSgam.Cor.Kendall.str = num2str( sprintf('%0.2f', Chandra.FradSgam.Cor.Kendall.val ) ) ;
Chandra.FradSgam.Cor.Spearman.str = num2str( sprintf('%0.2f', Chandra.FradSgam.Cor.Spearman.val ) ) ;

if figExportRequested, close all, figure('visible','off','Color','none'), else, figure, end
    hold on; box on; %colormap('cool');
    plot(DataX,DataY,'.','markerSize',20)
    xlim([1e1,1e5]);
    ylim([1e-7,1e-3]);
    set(gca,'xscale','log','fontsize',fontSize);
    set(gca,'yscale','log','fontsize',fontSize);
    xlabel('Radio Flux Density at 8.46 GHz: F^{Radio} [ \muJy ]', 'Interpreter', 'tex', 'fontSize', fontSize);
    ylabel('Gamma-Ray Fluence: S_{15-150} [ erg ]', 'Interpreter', 'tex', 'fontSize', fontSize);
    text( 0.05, 0.90, ['r = '   ,Chandra.FradSgam.Cor.Pearson.str], 'Units', 'normalized', 'HorizontalAlignment', 'left','Interpreter','tex', 'fontSize' , fontSize );
    text( 0.05, 0.80, ['\rho = ',Chandra.FradSgam.Cor.Spearman.str], 'Units', 'normalized', 'HorizontalAlignment', 'left','Interpreter','tex', 'fontSize' , fontSize );
    text( 0.05, 0.70, ['\tau = ',Chandra.FradSgam.Cor.Kendall.str], 'Units', 'normalized', 'HorizontalAlignment', 'left','Interpreter','tex', 'fontSize' , fontSize );
    set(gca, 'color', 'none', 'fontsize', fontSize);

if figExportRequested
    fileName = [outPath,'ChandraFoptSgam.png'];
    export_fig (fileName,'-m4 -transparent');
    hold off; close(gcf);
else
    hold off;
end

% plot peak radio flux at 8.46 GHz
UniqueFreq = unique(ChandraRadioPeakFlux.data(:,1));   % unique values
UniqueFreq(:,2) = countmember(unique(ChandraRadioPeakFlux.data(:,1)),ChandraRadioPeakFlux.data(:,1)); % unique values counts in the sample


% create Radio-Loud and Radio-Quiet ChandraTableOne
ChandraTableOne.Mask.ZoneKnown = ~isnan(ChandraTableOne.data(:,8));
ChandraTableOne.Mask.T90zKnown = ~isnan(ChandraTableOne.data(:,7)) & ChandraTableOne.data(:,7)>1.0; % this observer-frame data, also exclude SGRBs
ChandraTableOne.Mask.EisoKnown = ~isnan(ChandraTableOne.data(:,10));
ChandraTableOne.Mask.IsRadioDark = strcmp(deblank(ChandraTableOne.textdata(2:end,8)),'N');
ChandraTableOne.Mask.IsRadioBright = strcmp(deblank(ChandraTableOne.textdata(2:end,8)),'Y');
ChandraTableOne.Mask.AllKnown.Dark = ChandraTableOne.Mask.ZoneKnown ...
                                    & ChandraTableOne.Mask.EisoKnown ...
                                    & ChandraTableOne.Mask.T90zKnown ...
                                    & ChandraTableOne.Mask.IsRadioDark;
ChandraTableOne.Mask.AllKnown.Bright= ChandraTableOne.Mask.ZoneKnown ...
                                    & ChandraTableOne.Mask.EisoKnown ...
                                    & ChandraTableOne.Mask.T90zKnown ...
                                    & ChandraTableOne.Mask.IsRadioBright;
ChandraTableOne.Dark.count = sum( ChandraTableOne.Mask.AllKnown.Dark );
ChandraTableOne.Bright.count = sum( ChandraTableOne.Mask.AllKnown.Bright );

Radio.Type.Name = {'Dark','Bright'};
Radio.Type.Label = {'Radio-Dark','Radio-Bright'};
Radio.Type.count = length(Radio.Type.Name);
for irtype = 1:Radio.Type.count
    radioType = Radio.Type.Name{irtype};
    Mask = ChandraTableOne.Mask.AllKnown.(radioType);
    ChandraTableOne.(radioType).LogT90   = log( ChandraTableOne.data(Mask,7) );
    ChandraTableOne.(radioType).LogZone  = log( ChandraTableOne.data(Mask,8) + 1.0 );
    ChandraTableOne.(radioType).LogEiso  = log( ChandraTableOne.data(Mask,10) );
    ChandraTableOne.(radioType).LogT90z  = ChandraTableOne.(radioType).LogT90 ...
                                         - ChandraTableOne.(radioType).LogZone * kfacVal;
end

Var.Name = {'LogZone','LogEiso','LogT90z'};
Var.Label = {'(z+1)','E_{iso}','T_{90z}'};
Var.count = length(Var.Name);

for ivar = 1:Var.count-1
    varNameI = Var.Name{ivar};
    for jvar = ivar+1:Var.count

        if bivarFigExportRequested, close all, figure('visible','off','Color','none'), else figure, end
        hold on; box on; colormap('cool');
            varNameJ = Var.Name{jvar};
            for irtype = 1:Radio.Type.count
                radioType = Radio.Type.Name{irtype};
                plot( exp( ChandraTableOne.(radioType).(varNameI) ) ...
                    , exp( ChandraTableOne.(radioType).(varNameJ) ) ...
                    , '.', 'markerSize', 20 );

            end
            set(gca,'xscale','log','fontsize',fontSize);
            set(gca,'yscale','log','fontsize',fontSize);
            xlabel(Var.Label{ivar}, 'Interpreter', 'Tex', 'fontSize', fontSize);
            ylabel(Var.Label{jvar}, 'Interpreter', 'Tex', 'fontSize', fontSize);
            legend  ( Radio.Type.Label ...
                    ..., 'location' , Legend.Loc{counter} ...
                    , 'location' , 'northwest' ...
                    , 'fontSize' , fontSize ...
                    , 'color' , 'none' ...
                    )
            legend boxoff;
        set(gca, 'color', 'none', 'fontsize', fontSize);
        if bivarFigExportRequested
            fileName = [outPath,'Chandra',Var.Name{ivar},Var.Name{jvar},'.png'];
            export_fig (fileName,'-m4 -transparent');
            hold off; close(gcf);
        else
            hold off;
        end

    end
end


% now simulate the effects of Eiso cutoff threshold on the variable correlations
Log10Eiso.lower = 48.;
Log10Eiso.upper = 53.;
Log10Eiso.delta = 0.02;
Log10Eiso.count = (Log10Eiso.upper-Log10Eiso.lower) / Log10Eiso.delta + 1;
corType = 'Spearman';
for i = Log10Eiso.count:-1:1
    Log10Eiso.Value(i) = Log10Eiso.lower + (i-1)*Log10Eiso.delta;
    for irtype = 1:Radio.Type.count
        radioType = Radio.Type.Name{irtype};
        EisoMask = ChandraTableOne.(radioType).LogEiso > log(10)*Log10Eiso.Value(i);
        for ivar = 1:Var.count-1
            varNameI = Var.Name{ivar};
            for jvar = ivar+1:Var.count
                varNameJ = Var.Name{jvar};
                corName = [varNameI,varNameJ];
                [ Cor.(corType).(corName).(radioType).Sval(i) , Cor.(corType).(corName).Pval(i) ] = ...
                corr( ChandraTableOne.(radioType).(varNameI)(EisoMask) ...
                    , ChandraTableOne.(radioType).(varNameJ)(EisoMask) ...
                    , 'type' , corType );
                
            end
        end
        Cor.(corType).(corName).(radioType).Count(i) = sum(EisoMask);
    end
end

% plot correlation vs threshold
for ivar = 1:Var.count-1
    varNameI = Var.Name{ivar};
    for jvar = ivar+1:Var.count
        varNameJ = Var.Name{jvar};
        corName = [varNameI,varNameJ];

        if bivarFigExportRequested, close all, figure('visible','off','Color','none'), else figure, end
        hold on; box on; colormap('cool');
            RangeY.min = +1000;
            RangeY.max = -1000;
            for irtype = 1:Radio.Type.count
                radioType = Radio.Type.Name{irtype};
                plot( Log10Eiso.Value ...
                    , Cor.(corType).(corName).(radioType).Sval ...
                    , '-.', 'markerSize', 20 ...
                    , 'lineWidth', 2.5 );
                RangeY.min = min( RangeY.min , min(Cor.(corType).(corName).(radioType).Sval) );
                RangeY.max = max( RangeY.max , max(Cor.(corType).(corName).(radioType).Sval) );
            end
            minMaxDiff = RangeY.max - RangeY.min;
            ylim( [ max(-1,RangeY.min-0.5*minMaxDiff) , min(1,RangeY.max+0.5*minMaxDiff) ] );
            set(gca,'xscale','log','fontsize',fontSize);
            %set(gca,'yscale','log','fontsize',fontSize);
            xlabel('Isotropic Energy Cutoff on Sample: E_{iso}^{cut} [ ergs ]', 'Interpreter', 'Tex', 'fontSize', fontSize);
            ylabel( ['Spearman''s \rho:  ',Var.Label{ivar},' - ',Var.Label{jvar}] , 'Interpreter', 'Tex', 'fontSize', fontSize);
            legend  ( Radio.Type.Label ...
                    ..., 'location' , Legend.Loc{counter} ...
                    , 'location' , 'southwest' ...
                    , 'fontSize' , fontSize ...
                    , 'color' , 'none' ...
                    )
            legend boxoff;
            set(gca, 'color', 'none', 'fontsize', fontSize);
        if bivarFigExportRequested
            fileName = [outPath,'ChandraEisoCutoff',Var.Name{ivar},Var.Name{jvar},'.png'];
            export_fig (fileName,'-m4 -transparent');
            hold off; close(gcf);
        else
            hold off;
        end
    end
end



