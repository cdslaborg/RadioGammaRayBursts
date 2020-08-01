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
kstest.filename = 'Corksteststrap.mat';

return

% perform ksteststrapping
kstest.count = 100000;
%kstest.Cor.Type.Name = {'Pearson','Spearman','Kendall'};
kstest.Cor.Type.Name = {'Spearman'};
kstest.Cor.Type.count = length(kstest.Cor.Type.Name);

% construct correlation names
kstest.Var.Name = {'LogZone','LogEiso','LogDurz'};
kstest.Var.count = length(kstest.Var.Name);

if ksteststrapRequested

    kstest.Cor.Name = cell((kstest.Var.count^2-kstest.Var.count)/2,1);
    counter = 0;
    for ivar = 1:kstest.Var.count
        for jvar = ivar+1:kstest.Var.count
            counter = counter + 1;
            kstest.Cor.Name{counter} = [kstest.Var.Name{ivar},kstest.Var.Name{jvar}];
        end
    end
    kstest.Cor.count = counter;
    
    for icorType = 1:kstest.Cor.Type.count
    
        corType = kstest.Cor.Type.Name{icorType};
        disp([ 'kstest.Cor.Type: ' , corType ]);
    
        % Generate sample
        tic
        for irtype = 1:Radio.Type.count
            radioType = Radio.Type.Name{irtype};
            disp([ 'RadioType: ' , radioType ]);
            Indx = randi([1,Radio.(radioType).count],Radio.(radioType).count,kstest.count);
            % compute ksteststrap correlations
            counter = 0;
            for ivar = 1:kstest.Var.count
                ivarName = kstest.Var.Name{ivar};
                kstest.Var.(ivarName).(radioType).Avg = zeros(kstest.count,1);
                kstest.Var.(ivarName).(radioType).Std = zeros(kstest.count,1);
                for ikstest = 1:kstest.count
                    kstest.Var.(ivarName).(radioType).Avg(ikstest) = mean( Radio.(radioType).(ivarName).Val(Indx(:,ikstest)) );
                    kstest.Var.(ivarName).(radioType).Std(ikstest) =  std( Radio.(radioType).(ivarName).Val(Indx(:,ikstest)) );
                end
                for jvar = ivar+1:kstest.Var.count
                    jvarName = kstest.Var.Name{jvar};
                    counter = counter + 1;
                    corName = [ivarName(4:end),jvarName(4:end)];
                    kstest.Cor.(corType).(corName).(radioType).Sval = zeros(kstest.count,1);
                    kstest.Cor.(corType).(corName).(radioType).Pval = zeros(kstest.count,1);
                    for ikstest = 1:kstest.count
                        [ kstest.Cor.(corType).(corName).(radioType).Sval(ikstest) ...
                        , kstest.Cor.(corType).(corName).(radioType).Pval(ikstest) ] ...
                        = corr  ( Radio.(radioType).(ivarName).Val(Indx(:,ikstest)) ...
                                , Radio.(radioType).(jvarName).Val(Indx(:,ikstest)) ...
                                , 'type' , corType );
                    end
                end
            end
        end
        toc
    
    end

    save([kstest.outPath,kstest.filename],'kstest');
    
else

    %if ~exist('kstest','var')
        disp('loading ksteststrap Mat file, containing Structure ''kstest'' ...');
        load([kstest.outPath,kstest.filename])
    %else
    %    warning('Variable ''kstest'' already exists in MATLAB environment. Skipping kstest load from Hard Drive...');
    %end
end


% plot histograms and plots
kstest.Var.Label = {'(z+1)','E_{iso} [ ergs ]','T_{90z} [ s ]'};

kstest.Cor.Spearman.ZoneEiso.Hist.textLoc = [ 0.03, 0.77 ];
kstest.Cor.Spearman.ZoneDurz.Hist.textLoc = [ 0.97, 0.77 ];
kstest.Cor.Spearman.EisoDurz.Hist.textLoc = [ 0.03, 0.77 ];
kstest.Cor.Spearman.ZoneEiso.Hist.textAlign = 'left';
kstest.Cor.Spearman.ZoneDurz.Hist.textAlign = 'right';
kstest.Cor.Spearman.EisoDurz.Hist.textAlign = 'left';
kstest.Cor.Spearman.ZoneEiso.Hist.legendLoc = 'northwest';
kstest.Cor.Spearman.ZoneDurz.Hist.legendLoc = 'northeast';
kstest.Cor.Spearman.EisoDurz.Hist.legendLoc = 'northwest';

kstest.Var.LogZone.Hist.Avg.xlabel = 'Sample Mean of Log10 ( ';
kstest.Var.LogZone.Hist.Avg.textLoc = [ 0.97, 0.77 ];
kstest.Var.LogEiso.Hist.Avg.textLoc = [ 0.03, 0.77 ];
kstest.Var.LogDurz.Hist.Avg.textLoc = [ 0.03, 0.77 ];
kstest.Var.LogZone.Hist.Avg.textAlign = 'right';
kstest.Var.LogEiso.Hist.Avg.textAlign = 'left';
kstest.Var.LogDurz.Hist.Avg.textAlign = 'left';
kstest.Var.LogZone.Hist.Avg.legendLoc = 'northeast';
kstest.Var.LogEiso.Hist.Avg.legendLoc = 'northwest';
kstest.Var.LogDurz.Hist.Avg.legendLoc = 'northwest';
kstest.Var.LogZone.Hist.Avg.Limits = [0.40 0.70] * log(10);
kstest.Var.LogEiso.Hist.Avg.Limits = [51.0 53.6] * log(10);
kstest.Var.LogDurz.Hist.Avg.Limits = [0.30 1.70] * log(10);

kstest.Var.LogZone.Hist.Std.xlabel = 'Sample Standard Deviation of Log10 ( ';
kstest.Var.LogZone.Hist.Std.textLoc = [ 0.97, 0.77 ];
kstest.Var.LogEiso.Hist.Std.textLoc = [ 0.03, 0.77 ];
kstest.Var.LogDurz.Hist.Std.textLoc = [ 0.03, 0.77 ];
kstest.Var.LogZone.Hist.Std.textAlign = 'right';
kstest.Var.LogEiso.Hist.Std.textAlign = 'left';
kstest.Var.LogDurz.Hist.Std.textAlign = 'left';
kstest.Var.LogZone.Hist.Std.legendLoc = 'northeast';
kstest.Var.LogEiso.Hist.Std.legendLoc = 'northwest';
kstest.Var.LogDurz.Hist.Std.legendLoc = 'northwest';
kstest.Var.LogZone.Hist.Std.Limits = [0.13 0.35] * log(10);
kstest.Var.LogEiso.Hist.Std.Limits = [0.18 0.74] * log(10);
kstest.Var.LogDurz.Hist.Std.Limits = [0.20 0.65] * log(10);

for icorType = 1:1%1:kstest.Cor.Type.count
    corType = kstest.Cor.Type.Name{icorType};
    for ivar = 1:kstest.Var.count

        ivarName = kstest.Var.Name{ivar};

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        StatName = {'Avg','Std'};

        % plot Avg, Std histograms
        for iStat = 1:length(StatName)

            iStatName = StatName{iStat};

            if figExportRequested, close all, figure('visible','off','Color','none'), else, figure, end
                hold on; box on; %colormap('cool');
                LegendName = cell(Radio.Type.count,1);
                for irtype = 1:Radio.Type.count
                    radioType = Radio.Type.Name{irtype};
                    LegendName{irtype} = ['Radio-',radioType];
                    h = histogram( kstest.Var.(ivarName).(radioType).(iStatName) / log(10) ,'normalization','probability');
                    %h.BinWidth = 0.02;
                    h.BinLimits = kstest.Var.(ivarName).Hist.(iStatName).Limits / log(10);
                    h.NumBins = 100;
                    h.EdgeAlpha = 0; % remove bin edges
                    kstest.Var.(ivarName).(radioType).Hist.Height.Values = h.Values;
                    kstest.Var.(ivarName).(radioType).Hist.Height.CumSum = cumsum(h.Values);
                    kstest.Var.(ivarName).(radioType).avg = mean( kstest.Var.(ivarName).(radioType).(iStatName) / log(10) );
                    kstest.Var.(ivarName).(radioType).std =  std( kstest.Var.(ivarName).(radioType).(iStatName) / log(10) );
                    %xlabel([corType,' Correlation Strength: ',corName], 'fontSize', fontSize);
                end
                xlabel([kstest.Var.LogZone.Hist.(iStatName).xlabel, kstest.Var.Label{ivar},' )'], 'Interpreter', 'Tex', 'fontSize', fontSize);
                ylabel('Normalized ksteststrapping Count', 'Interpreter', 'Tex', 'fontSize', fontSize);
                %xlim([-0.82 0.82]);
                legend  ( LegendName ...
                        , 'location' , kstest.Var.(ivarName).Hist.(iStatName).legendLoc ...
                        , 'fontSize' , fontSize ...
                        , 'color' , 'none' ...
                        )
                legend boxoff;
    
                % Compute the odds of Radio-Bright correlations being higher than Radio-Dark correlations
                prob = sum  ( kstest.Var.(ivarName).Bright.Hist.Height.Values ...
                            .*kstest.Var.(ivarName).Dark.Hist.Height.CumSum );
                text( kstest.Var.(ivarName).Hist.(iStatName).textLoc(1) ...
                    , kstest.Var.(ivarName).Hist.(iStatName).textLoc(2) ...
                    , ['\pi ( {Bright} \geq {Dark} ) = ',sprintf('%.2f',round(prob,2))] ...
                    , 'Units', 'normalized', 'HorizontalAlignment', kstest.Var.(ivarName).Hist.(iStatName).textAlign ...
                    , 'Interpreter', 'tex', 'fontSize' , fontSize );
                text( kstest.Var.(ivarName).Hist.(iStatName).textLoc(1) ...
                    , kstest.Var.(ivarName).Hist.(iStatName).textLoc(2) - 0.1 ...
                    , [ '{Bright} = ' ...
                    , sprintf('%.2f',round(kstest.Var.(ivarName).Bright.avg,2)) ...
                    , ' \pm ' ...
                    , sprintf('%.2f',round(kstest.Var.(ivarName).Bright.std,2)) ...
                    ] ...
                    , 'Units', 'normalized', 'HorizontalAlignment', kstest.Var.(ivarName).Hist.(iStatName).textAlign ...
                    , 'Interpreter', 'tex', 'fontSize' , fontSize );
                text( kstest.Var.(ivarName).Hist.(iStatName).textLoc(1) ...
                    , kstest.Var.(ivarName).Hist.(iStatName).textLoc(2) - 0.2 ...
                    , [ '{Dark} = ' ...
                    , sprintf('%.2f',round(kstest.Var.(ivarName).Dark.avg,2)) ...
                    , ' \pm ' ...
                    , sprintf('%.2f',round(kstest.Var.(ivarName).Dark.std,2)) ...
                    ] ...
                    , 'Units', 'normalized', 'HorizontalAlignment', kstest.Var.(ivarName).Hist.(iStatName).textAlign ...
                    , 'Interpreter', 'tex', 'fontSize' , fontSize );
    
            set(gca, 'color', 'none', 'fontsize', fontSize);
            if figExportRequested
                fileName = [kstest.outPath,'kstestHist_',iStatName,'Log10',ivarName(4:end),'.png'];
                disp(['Generating Figure ',fileName]);
                export_fig(fileName,'-m4 -transparent');
                hold off; close(gcf);
            else
                hold off;
            end

        end % istat

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for jvar = ivar+1:kstest.Var.count

            jvarName = kstest.Var.Name{jvar};
            corName = [ivarName(4:end),jvarName(4:end)];

            % Correlation strength histograms
            if corFigExportRequested, close all, figure('visible','off','Color','none'), else, figure, end
                hold on; box on; %colormap('cool');
                LegendName = cell(Radio.Type.count,1);
                for irtype = 1:Radio.Type.count
                    radioType = Radio.Type.Name{irtype};
                    LegendName{irtype} = ['Radio-',radioType];
                    h = histogram(kstest.Cor.(corType).(corName).(radioType).Sval,'normalization','probability');
                    h.BinWidth = 0.02;
                    h.BinLimits = [-1.,1.];
                    h.NumBins = 100;
                    h.EdgeAlpha = 0; % remove bin edges
                    kstest.Cor.(corType).(corName).(radioType).Hist.Height.Values = h.Values;
                    kstest.Cor.(corType).(corName).(radioType).Hist.Height.CumSum = cumsum(h.Values);
                    kstest.Cor.(corType).(corName).(radioType).avg = mean( kstest.Cor.(corType).(corName).(radioType).Sval);
                    kstest.Cor.(corType).(corName).(radioType).std =  std( kstest.Cor.(corType).(corName).(radioType).Sval);
                    %xlabel([corType,' Correlation Strength: ',corName], 'fontSize', fontSize);
                end
                xlabel([corType,'''s Correlation Coefficient: ',kstest.Var.Label{ivar},' - ',kstest.Var.Label{jvar}], 'Interpreter', 'Tex', 'fontSize', fontSize);
                ylabel('Normalized ksteststrapping Count', 'Interpreter', 'Tex', 'fontSize', fontSize);
                xlim([-0.82 0.82]);
                legend  ( LegendName ...
                        , 'location' , kstest.Cor.(corType).(corName).Hist.legendLoc ...
                        , 'fontSize' , fontSize ...
                        , 'color' , 'none' ...
                        )
                legend boxoff;

                % Compute the odds of Radio-Bright correlations being higher than Radio-Dark correlations
                prob = sum  ( kstest.Cor.(corType).(corName).Bright.Hist.Height.Values ...
                            .*kstest.Cor.(corType).(corName).Dark.Hist.Height.CumSum );
                text( kstest.Cor.(corType).(corName).Hist.textLoc(1) ...
                    , kstest.Cor.(corType).(corName).Hist.textLoc(2) ...
                    , ['\pi ( \rho_{Bright} \geq \rho_{Dark} ) = ',sprintf('%.2f',round(prob,2))] ...
                    , 'Units', 'normalized', 'HorizontalAlignment', kstest.Cor.(corType).(corName).Hist.textAlign ...
                    , 'Interpreter', 'tex', 'fontSize' , fontSize );
                text( kstest.Cor.(corType).(corName).Hist.textLoc(1) ...
                    , kstest.Cor.(corType).(corName).Hist.textLoc(2) - 0.1 ...
                    , [ '\rho_{Bright} = ' ...
                      , sprintf('%.2f',round(kstest.Cor.(corType).(corName).Bright.avg,2)) ...
                      , ' \pm ' ...
                      , sprintf('%.2f',round(kstest.Cor.(corType).(corName).Bright.std,2)) ...
                      ] ...
                    , 'Units', 'normalized', 'HorizontalAlignment', kstest.Cor.(corType).(corName).Hist.textAlign ...
                    , 'Interpreter', 'tex', 'fontSize' , fontSize );
                text( kstest.Cor.(corType).(corName).Hist.textLoc(1) ...
                    , kstest.Cor.(corType).(corName).Hist.textLoc(2) - 0.2 ...
                    , [ '\rho_{Dark} = ' ...
                      , sprintf('%.2f',round(kstest.Cor.(corType).(corName).Dark.avg,2)) ...
                      , ' \pm ' ...
                      , sprintf('%.2f',round(kstest.Cor.(corType).(corName).Dark.std,2)) ...
                      ] ...
                    , 'Units', 'normalized', 'HorizontalAlignment', kstest.Cor.(corType).(corName).Hist.textAlign ...
                    , 'Interpreter', 'tex', 'fontSize' , fontSize );

            set(gca, 'color', 'none', 'fontsize', fontSize);
            if corFigExportRequested
                fileName = [kstest.outPath,'kstestHistCor_',ivarName(4:end),jvarName(4:end),'.png'];
                disp(['Generating Figure ',fileName]);
                export_fig(fileName,'-m4 -transparent');
                hold off; close(gcf);
            else
                hold off;
            end

        end % jvar

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end % ivar
end % icorType

