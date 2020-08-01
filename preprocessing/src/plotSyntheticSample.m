%close all;
%clear all;
format compact; format long;
filePath = mfilename('fullpath');
[scriptPath,fileName,fileExt] = fileparts(filePath); cd(scriptPath);
%addpath(genpath('../../../../../../lib/matlab/')) % lib codes
addpath(genpath('../../../../../libmatlab/')) % added by josh

% read Swift time table
datPath = '../data/';
Swift = importdata([datPath,'swiftTimeTable.xlsx']);
% plot histogram of log10(T90/Tr45)
figure; box on; hold on;
histogram( log10( Swift.data.Sheet2(:,1) ./ Swift.data.Sheet2(:,3) ) );
hold off;
butlerCeffmin = 10^(-0.59);
avgSwiftPartialCoding = 10^(-0.23);
avgSwiftT90OverTr45 = exp(mean(log( Swift.data.Sheet2(:,1) ./ Swift.data.Sheet2(:,3) ))) ;
medSwiftT90OverTr45 = exp(median(log( Swift.data.Sheet2(:,1) ./ Swift.data.Sheet2(:,3) ))) ;
disp( [ 'Butler''s Ceffmin = ', num2str(butlerCeffmin,15) ]);
disp( [ 'Swift mean partial coding = ', num2str(avgSwiftPartialCoding,15) ]);
disp( [ 'mean(T90/Tr45) = ',    num2str(avgSwiftT90OverTr45,15) ]);
disp( [ 'median(T90/Tr45) = ',  num2str(medSwiftT90OverTr45,15) ]);


outPath = '../out/';
kfacType = 'OneThird';
ChandraDataPlotRequested = 0;
if ChandraDataPlotRequested
    % read Chandra Data
    readChandraData;
%    close all;
    % restructure data into Lloyd Radio struct to avoid code redundancy
    for irtype = 1:Radio.Type.count
        radioType = Radio.Type.Name{irtype};
        Radio.(radioType).Zone = exp( ChandraTableOne.(radioType).LogZone );
        Radio.(radioType).Eiso = exp( ChandraTableOne.(radioType).LogEiso );
        Radio.(radioType).Durz = exp( ChandraTableOne.(radioType).LogDurz );
        Radio.MarkerSize.default = 25;
        Radio.MarkerSize.overlay = 12;
    end
    datasetName = 'Chandra';
else
    % import Radio data of Lloyd 2019
    Radio = readLloydRadioData(kfacType);
    datasetName = 'Lloyd';
%    close all;
end

ZoneRequested = 1;
bivarPlotsRequested = 1;
histFigExportRequested = 0;
bivarFigExportRequested = 0;
fontSize = 13;
myYellow = [1,1,0];
myGrey = [0 0 0]; % dtermines the grey color strength, higher means whiter
%myGrey = [0.6 0.6 0.6]; % dtermines the grey color strength, higher means whiter
%SynSam.Path.root = '../../../../20181213_BatseLgrbRedshift/git/SyntheticSample/';
SynSam.Path.root = '../../../../../BatseLgrbRedshiftCatalog/git/___SyntheticSample___/'; %added by Josh
SynSam.Path.output = [SynSam.Path.root,'winx64/intel/release/static/serial/bin/out/kfac',kfacType,'/'];
SynSam.Path.input = [SynSam.Path.root,'in/'];

% import BATSE data
Dummy = importdata([SynSam.Path.input,'batse_1366_lgrb_pbol_epk_sbol(0.001,20000).txt']);
Batse.LogData.Obs = [ Dummy.data(:,2) ... % logPbol
                    , Dummy.data(:,4) ... % logEpk
                    , Dummy.data(:,3) ... % logSbol
                    , Dummy.data(:,8) ... % logDur
                    ];
Batse.Data.Obs = exp(Batse.LogData.Obs);
Batse.ngrb = length(Batse.Data.Obs(:,1));
Batse.Trigger = Dummy.data(:,1);

% import Amati relation data
Ghirlanda08 = importdata([SynSam.Path.input,'AmatiRelationGhirlanda2008.txt']);

% read synthetic data
ZModel.ID = {'H06','L08','B10'};    %,'B04','F00','Y04'};
%ZModel.ID = {'H06'};    %,'B04','F00','Y04'};
%ZModel.ID = {'L08'};    %,'B04','F00','Y04'};
%ZModel.ID = {'B10'};    %,'B04','F00','Y04'};
ZModel.count = length(ZModel.ID);
ZModel.Ref = { 'This Work (with H06 Rate)' ...
             , 'This Work (with L08 Rate)' ...
             , 'This Work (with B10 Rate)' ...
             ..., 'Band (2004)' ...
             ..., 'Fenimore (2000)' ...
             ..., 'Yonetoku (2004)' ...
             };
ZPath = cell(ZModel.count,1);
for imodel = 1:ZModel.count
    ZPath{imodel} = [SynSam.Path.root,'../zestimation/winx64/intel/release/static/serial/kfac',kfacType,'/',ZModel.ID{imodel},'/bin/out/'];
end

nvar = 4;
VarPair = zeros(20,2);
VarPair(1:6,1:2)    = [ 1, 2 ... logLiso-logEpkz
                      ; 1, 3 ... logLiso-logEiso
                      ; 1, 4 ... logLiso-logDurz
                      ; 3, 2 ... logEiso-logEpkz
                      ; 3, 4 ... logEiso-logDurz
                      ; 4, 2 ... logDurz-logEpkz
                      ];
VarPair(7:12,1:2) = VarPair(1:6,1:2) + 4;
VarPair(13:20,1:2)  = [ 9, 1 ... redshift-logLiso
                      ; 9, 2 ... redshift-logEpkz
                      ; 9, 3 ... redshift-logEiso
                      ; 9, 4 ... redshift-logDurz
                      ; 9, 5 ... redshift-logPbol
                      ; 9, 6 ... redshift-logEpk
                      ; 9, 7 ... redshift-logSbol
                      ; 9, 8 ... redshift-logDur
                      ];

VarName = {'Liso','Epkz','Eiso','Durz','Pbol','Epk','Sbol','T90','Redshift','RedshiftPlusOne'};%,'LumDis'};
AxisLabel = { 'Bolometric Peak Luminosity: L_{iso} [ ergs/s ]' ...
            , 'Intrinsic Spectral Peak Energy: E_{pz} [ keV ]' ...
            , 'Isotropic Radiated Energy: E_{iso} [ ergs ]' ...
            , 'Intrinsic Duration: T_{90z} [ s ]' ...
            , 'Bolometric Peak Flux: P_{bol} [ ergs/s/cm^2 ]' ...
            , 'Observed Spectral Peak Energy: E_{p} [ keV ]' ...
            , 'Bolometric Fluence: S_{bol} [ ergs/cm^2 ]' ...
            , 'Observed Duration: T_{90} [ s ]' ...
            , 'Redshift: z' ...
            , 'z + 1' ...
            };
Log10VarLim =   [ 46, 56 ... log10Liso
                ; 0 , 5  ... log10Epkz
                ; 46, 56 ... log10Eiso
                ; -1, 4  ... log10Durz
                ; -13, -2 ... log10Pbol
                ; -1 , 4 ... log10Epk
                ; -15, 0 ... log10Sbol
                ; -1, 4  ... log10Dur
                ; -1, 1.5  ... redshift
                ; 0, 2  ... redshift+1
                ];
VarLim =        [ 5.e46, 1.e55  ... log10Liso
                ; 1.e0 , 3.e4  ... log10Epkz
                ; 1.e46, 3.e56 ... log10Eiso
                ; 5.e-2, 2.e3  ... log10Durz
                ; 1.e-13, 1.e-2 ... log10Pbol
                ; 2.e-1 , 1.e4 ... log10Epk
                ; 1.e-13, 1.e0 ... log10Sbol
                ; 1.e-1, 5.e3  ... log10Dur
                ; 1.e-1, 3.5e1  ... redshift
                ; 1.e0, 3.6e1  ... redshift+1
                ];

outPathSynSam = [outPath,'SynSam_kfac',kfacType,'/'];
if ~exist(['kfac',kfacType],'dir'); mkdir(outPathSynSam); end

nVarPair = length(VarPair);
for iVarPair = [5,15,16] %1:nVarPair

    disp(['processing variable pairs # ', num2str(iVarPair)]);

    for imodel = 1:ZModel.count

        if ~isfield(ZModel,ZModel.ID{imodel})

            ZModel.(ZModel.ID{imodel}).Synthetic = importdata([SynSam.Path.output,'syntheticSample',ZModel.ID{imodel},'.csv']);
            ZModel.(ZModel.ID{imodel}).Synthetic.data(:,1:8) = exp( ZModel.(ZModel.ID{imodel}).Synthetic.data(:,1:8) );
% 
%             % read redshift grid data
%             MPC2CM = 3.09e24; % 1 Mega Parsec = MPC2CM centimeters.
%             LOGMPC2CMSQ4PI = log(4.0*pi) + 2.0*log(MPC2CM);  % log(MegaParsec2centimeters).
%             ZModel.(ZModel.ID{imodel}).zstat = importdata([ZPath{imodel},'batse_zstat.txt']);
%             ZModel.(ZModel.ID{imodel}).zstat.count = length(ZModel.(ZModel.ID{imodel}).zstat.data(:,1));
%             % terms to map observer-frame data to rest-frame, equivalent to logRestFrameProp-logObserverFrameProp
%             ZModel.(ZModel.ID{imodel}).zstat.IntObsDiff = zeros(ZModel.(ZModel.ID{imodel}).zstat.count,4);
%             ZModel.(ZModel.ID{imodel}).zstat.IntObsDiff(:,1) = LOGMPC2CMSQ4PI + 2*getLogLumDisWicMPC(ZModel.(ZModel.ID{imodel}).zstat.data(:,2)+1.0);
%             ZModel.(ZModel.ID{imodel}).zstat.IntObsDiff(:,2) = log(ZModel.(ZModel.ID{imodel}).zstat.data(:,2)+1);
%             ZModel.(ZModel.ID{imodel}).zstat.IntObsDiff(:,3) = ZModel.(ZModel.ID{imodel}).zstat.IntObsDiff(:,1) - ZModel.(ZModel.ID{imodel}).zstat.IntObsDiff(:,2);
%             ZModel.(ZModel.ID{imodel}).zstat.IntObsDiff(:,4) = -ZModel.(ZModel.ID{imodel}).zstat.IntObsDiff(:,2)*0.66666666666;
% 
%             % compute Batse GRB Intrinsic Prob
%             Batse.LogData.Int = Batse.LogData.Obs + ZModel.(ZModel.ID{imodel}).zstat.IntObsDiff;
%             Batse.Data.Int = exp( Batse.LogData.Int );

        end

        synBegin = 1;
        synEnd = length(ZModel.(ZModel.ID{imodel}).Synthetic.data(:,1));%13660;
        ZModel.(ZModel.ID{imodel}).Synthetic.count = synEnd-synBegin+1;

        if bivarPlotsRequested

            if bivarFigExportRequested
                figure('visible','off','Color','none');
            else
                figure;
            end
            hold on; box on;

            colormap('cool');
            % get the matrix containing that colormap, then flip the matrix to invert the colormap.
            %cmap = colormap; cmap = flipud(cmap); colormap(cmap);

            %Mask = zeros( length(ZModel.(ZModel.ID{imodel}).Synthetic.data(:,1)) , 1 ); Mask(synBegin:synEnd) = 1;
            Mask = ZModel.(ZModel.ID{imodel}).Synthetic.data(synBegin:synEnd,14)>-1; % Josh made this change 10 to 14, plots using SWift detection threshold
            %Mask = ZModel.(ZModel.ID{imodel}).Synthetic.data(synBegin:synEnd,end)>0.5;
            %Mask = log10 (     ZModel.(ZModel.ID{imodel}).Synthetic.data(:,11) ) ... Nbol
            %     - log10 ( 0.5*ZModel.(ZModel.ID{imodel}).Synthetic.data(:,8) ) > log10(3);

            VarPairX = VarPair(iVarPair,1);
            DataX = ( ZModel.(ZModel.ID{imodel}).Synthetic.data( Mask , VarPairX ) );
            if VarPairX==9 && ZoneRequested % corresponding to redshift variable
                DataX = DataX + 1.0;
                VarPairX = VarPairX + 1;
            end

            VarPairY = VarPair(iVarPair,2);
            DataY = ( ZModel.(ZModel.ID{imodel}).Synthetic.data( Mask , VarPairY ) );
            if VarPairY==9 && ZoneRequested % corresponding to redshift variable
                DataY = DataY + 1.0;
                VarPairY = VarPairY + 1;
            end

            
            if iVarPair ~=16
                scatter ( DataX ...
                        , DataY ...
                        ..., 0.75*ones(ZModel.(ZModel.ID{imodel}).Synthetic.count,1) ...
                        , 0.75*ones(sum(Mask),1) ...
                        ..., ZModel.(ZModel.ID{imodel}).Synthetic.data(Mask,10) ...
                        ..., log10 ( ZModel.(ZModel.ID{imodel}).Synthetic.data(Mask,11) ) ... Nbol
                        ...- log10 ( 0.5*ZModel.(ZModel.ID{imodel}).Synthetic.data(Mask,8) ) ... duration
                        ..., ZModel.(ZModel.ID{imodel}).Synthetic.data(Mask,13) ...
                        , ZModel.(ZModel.ID{imodel}).Synthetic.data(Mask,14) ... %josh changed 13 to 14
                        , '.' ..., 'filled' ...
                        )
            else
                scatter ( DataX ...
                        , DataY ...
                        ..., 0.75*ones(ZModel.(ZModel.ID{imodel}).Synthetic.count,1) ...
                        , 3.25*ones(sum(Mask),1) ...
                        ..., ZModel.(ZModel.ID{imodel}).Synthetic.data(Mask,10) ...
                        ..., log10 ( ZModel.(ZModel.ID{imodel}).Synthetic.data(Mask,11) ) ... Nbol
                        ...- log10 ( 0.5*ZModel.(ZModel.ID{imodel}).Synthetic.data(Mask,8) ) ... duration
                        ..., ZModel.(ZModel.ID{imodel}).Synthetic.data(Mask,13) ...
                        , ZModel.(ZModel.ID{imodel}).Synthetic.data(Mask,14) ... %josh changed 13 to 14
                        , '.' ..., 'filled' ...
                        )
            end    
            
            
            CBar = colorbar;
            %CBar.Label.String = 'Probability of Detection by BATSE LADs';
            %CBar.Label.String = 'Probability of Detection by Swift Based on B07';
            CBar.Label.String = 'Probability of Detection by Swift Based on B10';
            CBar.Label.Interpreter = 'tex';
            CBar.Label.FontSize = fontSize;
            %plot( ZModel.(ZModel.ID{imodel}).Synthetic.data(synBegin:synEnd,3) ...
            %    , ZModel.(ZModel.ID{imodel}).Synthetic.data(synBegin:synEnd,2) ...
            %    , '.', 'MarkerSize', 0.5 ...
            %    )
            %if iVarPair<7
            %    % plot BATSE data
            %    scatter ( Batse.Data.Int(:,VarPair(iVarPair,1)) ...
            %            , Batse.Data.Int(:,VarPair(iVarPair,2)) ...
            %            , 30*ones(Batse.ngrb,1) ...
            %            , 'black' ...
            %            ..., ColorMapCloud(1:instanceCounter,1:3) ...
            %            , '.' ..., 'filled' ...
            %            );
            %    set(gca,'xscale','log','fontsize',fontSize);
            %    set(gca,'yscale','log','fontsize',fontSize,'YDir','normal');

            % add Amati relation data to Eiso-Epkz plot
            if iVarPair==4
                errorbar( Ghirlanda08.data(:,5) ...
                        , Ghirlanda08.data(:,3) ...
                        , Ghirlanda08.data(:,4) ...
                        , Ghirlanda08.data(:,4) ...
                        , Ghirlanda08.data(:,6) ...
                        , Ghirlanda08.data(:,6) ...
                        , '.', 'MarkerSize', 11 ...
                        , 'LineWidth', 1 ...
                        , 'color', [0.2 0.2 0.2] ... 'black' ...
                        , 'CapSize', 0 ...
                        );
            %elseif iVarPair>6
            %    plot( Batse.Data.Obs(:,VarPair(iVarPair-6,1)) ...
            %        , Batse.Data.Obs(:,VarPair(iVarPair-6,2)) ...
            %        , '.', 'MarkerSize', 6 ...
            %        , 'color', 'black' ...
            %        );
            end

            % add Lloyd data to Eiso-Durz plot
            if iVarPair==5
               % 'orange' [0.9100 0.4100 0.1700]
                plot( Radio.Bright.Eiso.Val, Radio.Bright.Durz.Val, '.', 'MarkerSize', Radio.MarkerSize.overlay, 'color', myYellow);
                plot(   Radio.Dark.Eiso.Val,   Radio.Dark.Durz.Val, '.', 'MarkerSize', Radio.MarkerSize.overlay, 'color', myGrey);
            end

            % add Lloyd data to redshift-Eiso plot
            if iVarPair==15
                plot( Radio.Bright.Zone.Val, Radio.Bright.Eiso.Val, '.', 'MarkerSize', Radio.MarkerSize.overlay, 'color', myYellow);
                plot(   Radio.Dark.Zone.Val,   Radio.Dark.Eiso.Val, '.', 'MarkerSize', Radio.MarkerSize.overlay, 'color', myGrey);
            end

            % add Lloyd data to redshift-Durz plot
            if iVarPair==16
                plot( Radio.Bright.Zone.Val, Radio.Bright.Durz.Val, '.', 'MarkerSize', Radio.MarkerSize.overlay, 'color', myYellow);
                plot(   Radio.Dark.Zone.Val,   Radio.Dark.Durz.Val, '.', 'MarkerSize', Radio.MarkerSize.overlay, 'color', myGrey);
            end

            set(gca,'xscale','log','fontsize',fontSize);
            set(gca,'yscale','log','fontsize',fontSize,'YDir','normal');

            xlim( VarLim(VarPairX,:) );
            ylim( VarLim(VarPairY,:) );
            xlabel(AxisLabel{VarPairX}, 'Interpreter', 'tex', 'fontSize', fontSize);
            ylabel(AxisLabel{VarPairY}, 'Interpreter', 'tex', 'fontSize', fontSize);
            legend  ( {['Simulated LGRB: ',ZModel.ID{imodel},' Rate']} ...
                    , 'location' , 'southwest' ...
                    , 'fontSize' , fontSize ...
                    , 'color' , 'none' ...
                    , 'box', 'off' ...
                    );

            set(gca, 'color', 'none', 'fontsize', fontSize);
            if bivarFigExportRequested
                fileName = [outPathSynSam,VarName{VarPairX},VarName{VarPairY},ZModel.ID{imodel},'_',datasetName,'.png'];
                export_fig (fileName,'-m4 -transparent');
                hold off; close(gcf);
            else
                hold off;
            end

        end % bivarPlotsRequested

    end % imodel

end


% find the moving average of Durz vs. Eiso for the detectable sample
if ~exist('Boot','var'), bootstrap, end   % get the standard deviations of the sample averages
Mask = ZModel.(ZModel.ID{imodel}).Synthetic.data(:,14)>0.5; % detectable sample %%% josh changed 13 to 14
Detectable.count = sum(Mask);
Detectable.sampleSize = 5000;
Detectable.LogEiso = log(ZModel.(ZModel.ID{imodel}).Synthetic.data(Mask,3));
Detectable.LogDurz = log(ZModel.(ZModel.ID{imodel}).Synthetic.data(Mask,4));
[ Detectable.Sorted.LogEiso , Indx ] = sort(Detectable.LogEiso);
Detectable.Sorted.LogDurz = Detectable.LogDurz(Indx);
Detectable.MovMean.LogEiso = movmean(Detectable.Sorted.LogEiso,Detectable.sampleSize);
Detectable.MovMean.LogDurz = movmean(Detectable.Sorted.LogDurz,Detectable.sampleSize);
figure; hold on; box on;
plot( exp(Detectable.MovMean.LogEiso), exp(Detectable.MovMean.LogDurz) ...
    , 'linewidth', 2 ...
    )
%plot( exp(Detectable.Sorted.LogEiso), exp(Detectable.MovMean.LogDurz), '.')
xscale = 'log';
yscale = 'log';
if ~exist('Nicole','var'), Nicole = readLloydRadioData; end
RadioType = {'Dark','Bright'};
RadioTypeColor = {'black','red'};
%VarType = {'Zone','Eiso','Durz'};
%for ivar = 1:length(VarType)
%    var = VarType{ivar};
    for irtype = 1:length(RadioType)
        rtype = RadioType{irtype};
        %yneg = exp( Nicole.(rtype).LogDurz.avg - Nicole.(rtype).LogDurz.std );
        %ypos = exp( Nicole.(rtype).LogDurz.avg + Nicole.(rtype).LogDurz.std );
        %xneg = exp( Nicole.(rtype).LogEiso.avg - Nicole.(rtype).LogEiso.std );
        %xpos = exp( Nicole.(rtype).LogEiso.avg + Nicole.(rtype).LogEiso.std );
        yneg = exp( Nicole.(rtype).LogDurz.avg ) - exp( Nicole.(rtype).LogDurz.avg - std(Boot.Var.LogDurz.(rtype).Avg) );
        ypos = exp( Nicole.(rtype).LogDurz.avg ) - exp( Nicole.(rtype).LogDurz.avg + std(Boot.Var.LogDurz.(rtype).Avg) ); ypos = -ypos
        xneg = exp( Nicole.(rtype).LogEiso.avg ) - exp( Nicole.(rtype).LogEiso.avg - std(Boot.Var.LogEiso.(rtype).Avg) );
        xpos = exp( Nicole.(rtype).LogEiso.avg ) - exp( Nicole.(rtype).LogEiso.avg + std(Boot.Var.LogEiso.(rtype).Avg) ); xpos = -xpos
        errorbar( exp( Nicole.(rtype).LogEiso.avg ) ...
                , exp( Nicole.(rtype).LogDurz.avg ) ...
                , yneg, ypos, xneg, xpos ...
                , '.', 'markersize', 20 ...
                , 'color', RadioTypeColor{irtype} );
    end
%end
set(gca,'xscale',xscale);
set(gca,'yscale',yscale);
xlim([1e51 1e55]);
ylim([5 50]);
hold off;





return





% get the hitogram data for z+1
ZDataAll.X = cell(ZModel.count,1);
ZDataAll.Y = cell(ZModel.count,1);
ZDataDetected.X = cell(ZModel.count,1);
ZDataDetected.Y = cell(ZModel.count,1);
ZRange = [-1 2]*log(10);
%ZRange = [0.1 100];
zBinWidth = 0.1;
for imodel = 1:ZModel.count

    % get histogram data for all data first:
    Mask = ZModel.(ZModel.ID{imodel}).Synthetic.data(:,end)>-0.5;
   %zdatadum = log( ZModel.(ZModel.ID{imodel}).Synthetic.data( Mask , end-1 ) + 1 );
    zdatadum = log( ZModel.(ZModel.ID{imodel}).Synthetic.data( Mask , end-1 ) );
    h = histogram( zdatadum );%, 'Normalization' , 'probability' );
    %h.BinLimits = [-1 2]*log(10);
    h.BinLimits = ZRange;
    h.BinWidth = zBinWidth;
    ZDataAll.X{imodel} = h.BinEdges(1:end-1);
    ZDataAll.Y{imodel} = h.Values * 1366 / length(ZModel.(ZModel.ID{imodel}).Synthetic.data(:,1)); %h.BinCounts;
    close(gcf);

    % get histogram data for BATSE detected LGRBs:
    Mask = ZModel.(ZModel.ID{imodel}).Synthetic.data(:,end)>0.5;
   %zdatadum = log( ZModel.(ZModel.ID{imodel}).Synthetic.data( Mask , end-1 ) + 1 );
    zdatadum = log( ZModel.(ZModel.ID{imodel}).Synthetic.data( Mask , end-1 ) );
    h = histogram( zdatadum );%, 'Normalization' , 'probability' );
    h.BinLimits = ZRange;
    h.BinWidth = zBinWidth;
    ZDataDetected.X{imodel} = h.BinEdges(1:end-1);
    ZDataDetected.Y{imodel} = h.Values * 1366 / length(ZModel.(ZModel.ID{imodel}).Synthetic.data(:,1)); %h.BinCounts;
    close(gcf);

end

% now make the plot of redshift histogram

if histFigExportRequested
    figure('visible','off','Color','none');
else
    figure;
end
hold on; box on;

    for imodel = 1:ZModel.count
        %disp(['plotting histogram for variable # ', num2str(iVarName)]);
        % reset MATLAB default color order
        ax = gca; ax.ColorOrderIndex = imodel;
        plot( exp(ZDataAll.X{imodel}) ...
            , ZDataAll.Y{imodel} ...
            , '-' , 'linewidth' , 2.5 );
        ax = gca; ax.ColorOrderIndex = imodel;
        plot( exp(ZDataDetected.X{imodel}) ...
            , ZDataDetected.Y{imodel} ...
            , ':', 'linewidth' , 2.5 );
        %hold on;
        %[f,xi] = ksdensity(zdatadum);
        %plot(exp(xi),f,'linewidth',3);
        %VarLim
    end

    %xlim( [1 50] );
    xlim( [0.1 20] );
    %ylim( VarLim(VarPair(iVarPair,2),:) );
    %xlabel('Redshift: z + 1', 'Interpreter', 'tex', 'fontSize', fontSize);
    xlabel('Redshift: z', 'Interpreter', 'tex', 'fontSize', fontSize);
    %ylabel('Normalized LGRB Rate: dN / d(log_{10}(z+1))  ', 'Interpreter', 'tex', 'fontSize', fontSize);
    ylabel('Normalized LGRB Rate: dN / d(log_{10}(z))  ', 'Interpreter', 'tex', 'fontSize', fontSize);
    legend  ( { [ZModel.ID{1},' Cosmic Rate'] ...
              , [ZModel.ID{1},' Observed Rate'] ...
              %, [ZModel.ID{2},' Cosmic Rate'] ...
              %, [ZModel.ID{2},' Observed Rate'] ...
              %, [ZModel.ID{3},' Cosmic Rate'] ...
              %, [ZModel.ID{3},' Observed Rate'] ...
              } ...
            , 'location' , 'northwest' ...
            , 'fontSize' , fontSize ...
            , 'color' , 'none' ...
            , 'box', 'off' ...
            );
    set(gca,'xscale','log','fontsize',fontSize);
    set(gca,'yscale','linear','fontsize',fontSize,'YDir','normal');
    set(gca,'color','none');

zdistOutPath = 'zdist/';
if ~exist(zdistOutPath,'dir'); mkdir(zdistOutPath); end
if histFigExportRequested
    fileName = [zdistOutPath,'zdist.png'];
    export_fig (fileName,'-m4 -transparent');
    hold off; close(gcf);
else
    hold off;
end

% get statistics of different terms
MPC2CM = 3.09e24;   % 1 Mega Parsec = MPC2CM centimeters.
LOGMPC2CMSQ4PI = log(4.0*pi) + 2.0*log(MPC2CM);     % log(MegaParsec2centimeters).
Intrinsic = cell(ZModel.count,1);
Observed = cell(ZModel.count,1);
StaDev = cell(ZModel.count,1);
Resolution = zeros(ZModel.count,1);
Accuracy = cell(ZModel.count,1);
for imodel = 1:ZModel.count

    % cosmic rate events and statistics
    Intrinsic{imodel}.LogZone = ( ZModel.(ZModel.ID{imodel}).Synthetic.data( : , end-1 ) + 1.0 );
    Intrinsic{imodel}.LogLisoPbolDiff = LOGMPC2CMSQ4PI + 2 * getLogLumDisWicMPC( Intrinsic{imodel}.LogZone );
    Intrinsic{imodel}.LogZone = log( Intrinsic{imodel}.LogZone );
    Intrinsic{imodel}.LogEisoSbolDiff = Intrinsic{imodel}.LogLisoPbolDiff - Intrinsic{imodel}.LogZone;
    for ivar = 1:8
        fieldName = ['Log',VarName{ivar}];
        Intrinsic{imodel}.(fieldName) = log( ZModel.(ZModel.ID{imodel}).Synthetic.data( : , ivar ) );
    end

    % observed rate events and statistics
    Mask = ZModel.(ZModel.ID{imodel}).Synthetic.data(:,end)>0.5;
    Observed{imodel}.LogZone = Intrinsic{imodel}.LogZone(Mask);
    Observed{imodel}.LogLisoPbolDiff = Intrinsic{imodel}.LogLisoPbolDiff(Mask);
    Observed{imodel}.LogEisoSbolDiff = Intrinsic{imodel}.LogEisoSbolDiff(Mask);
    for ivar = 1:8
        fieldName = ['Log',VarName{ivar}];
        Observed{imodel}.(fieldName) = Intrinsic{imodel}.(fieldName)(Mask);
    end

    StaDev{imodel}.Int.log10Zone     = log10(exp(1)) * std(Intrinsic{imodel}.LogZone);
    StaDev{imodel}.Obs.log10Zone     = log10(exp(1)) * std(Observed{imodel}.LogZone);
    StaDev{imodel}.Int.log10LisoPbolDiff = log10(exp(1)) * std(Intrinsic{imodel}.LogLisoPbolDiff);
    StaDev{imodel}.Obs.log10LisoPbolDiff = log10(exp(1)) * std(Observed{imodel}.LogLisoPbolDiff);
    StaDev{imodel}.Int.log10EisoSbolDiff = log10(exp(1)) * std(Intrinsic{imodel}.LogEisoSbolDiff);
    StaDev{imodel}.Obs.log10EisoSbolDiff = log10(exp(1)) * std(Observed{imodel}.LogEisoSbolDiff);
    for ivar = 1:8
        fieldName = ['Log',VarName{ivar}];
        fieldNameStaDev = ['log10',VarName{ivar}];
        StaDev{imodel}.Int.(fieldNameStaDev) = log10(exp(1)) * std(Intrinsic{imodel}.(fieldName));
        StaDev{imodel}.Obs.(fieldNameStaDev) = log10(exp(1)) * std(Observed{imodel}.(fieldName));
    end

    % divide the standard deviation of the overall redshift distribution to
    % the mean of the predicted 90% intervals of individual Batse redshifts
    Resolution(imodel) = std(ZModel.(ZModel.ID{imodel}).Synthetic.data( Mask , end-1 )) / mean( ZModel.(ZModel.ID{imodel}).zstat.data(:,8) - ZModel.(ZModel.ID{imodel}).zstat.data(:,4) );
    Accuracy{imodel}.Obs.log10Liso = StaDev{imodel}.Obs.log10LisoPbolDiff / StaDev{imodel}.Obs.log10Liso;
    Accuracy{imodel}.Obs.log10Eiso = StaDev{imodel}.Obs.log10EisoSbolDiff / StaDev{imodel}.Obs.log10Eiso;
    Accuracy{imodel}.Obs.log10Epkz = StaDev{imodel}.Obs.log10Zone     / StaDev{imodel}.Obs.log10Epkz;
    Accuracy{imodel}.Obs.log10Durz = StaDev{imodel}.Obs.log10Zone     / StaDev{imodel}.Obs.log10Durz;
    Accuracy{imodel}.Int.log10Liso = StaDev{imodel}.Int.log10LisoPbolDiff / StaDev{imodel}.Int.log10Liso;
    Accuracy{imodel}.Int.log10Eiso = StaDev{imodel}.Int.log10EisoSbolDiff / StaDev{imodel}.Int.log10Eiso;
    Accuracy{imodel}.Int.log10Epkz = StaDev{imodel}.Int.log10Zone     / StaDev{imodel}.Int.log10Epkz;
    Accuracy{imodel}.Int.log10Durz = StaDev{imodel}.Int.log10Zone     / StaDev{imodel}.Int.log10Durz;

end

