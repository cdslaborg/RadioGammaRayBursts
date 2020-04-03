% get the confidence intervals on the radio-loud and radio-quiet GRB
% correlations: logEiso-logDurz, logZone-logEiso, logZone-logDurz
%close all;
%clear all;
format compact; format long;
filePath = mfilename('fullpath');
[scriptPath,fileName,fileExt] = fileparts(filePath); cd(scriptPath);
addpath(genpath('../../../../lib/matlab/')) % lib codes

bootstrapRequested = 0;
varFigExportRequested = 0;
corFigExportRequested = 0;
fontSize = 13;
kfacType = 'OneThird';

% import Radio data of Lloyd 2019
Radio = readLloydRadioData(kfacType);

Boot.outPath = '../out/';
Boot.datPath = '../data/';
Boot.filename = 'CorBootstrap.mat';

% perform bootstrapping
Boot.count = 100000;
%Boot.Cor.Type.Name = {'Pearson','Spearman','Kendall'};
Boot.Cor.Type.Name = {'Spearman'};
Boot.Cor.Type.count = length(Boot.Cor.Type.Name);

% construct correlation names
Boot.Var.Name = {'LogZone','LogEiso','LogDurz'};
Boot.Var.count = length(Boot.Var.Name);

if bootstrapRequested

    Boot.Cor.Name = cell((Boot.Var.count^2-Boot.Var.count)/2,1);
    counter = 0;
    for ivar = 1:Boot.Var.count
        for jvar = ivar+1:Boot.Var.count
            counter = counter + 1;
            Boot.Cor.Name{counter} = [Boot.Var.Name{ivar},Boot.Var.Name{jvar}];
        end
    end
    Boot.Cor.count = counter;
    
    for icorType = 1:Boot.Cor.Type.count
    
        corType = Boot.Cor.Type.Name{icorType};
        disp([ 'Boot.Cor.Type: ' , corType ]);
    
        % Generate sample
        tic
        for irtype = 1:Radio.Type.count
            radioType = Radio.Type.Name{irtype};
            disp([ 'RadioType: ' , radioType ]);
            Indx = randi([1,Radio.(radioType).count],Radio.(radioType).count,Boot.count);
            % compute bootstrap correlations
            counter = 0;
            for ivar = 1:Boot.Var.count
                ivarName = Boot.Var.Name{ivar};
                Boot.Var.(ivarName).(radioType).Avg = zeros(Boot.count,1);
                Boot.Var.(ivarName).(radioType).Std = zeros(Boot.count,1);
                for iboot = 1:Boot.count
                    Boot.Var.(ivarName).(radioType).Avg(iboot) = mean( Radio.(radioType).(ivarName).Val(Indx(:,iboot)) );
                    Boot.Var.(ivarName).(radioType).Std(iboot) =  std( Radio.(radioType).(ivarName).Val(Indx(:,iboot)) );
                end
                for jvar = ivar+1:Boot.Var.count
                    jvarName = Boot.Var.Name{jvar};
                    counter = counter + 1;
                    corName = [ivarName(4:end),jvarName(4:end)];
                    Boot.Cor.(corType).(corName).(radioType).Sval = zeros(Boot.count,1);
                    Boot.Cor.(corType).(corName).(radioType).Pval = zeros(Boot.count,1);
                    for iboot = 1:Boot.count
                        [ Boot.Cor.(corType).(corName).(radioType).Sval(iboot) ...
                        , Boot.Cor.(corType).(corName).(radioType).Pval(iboot) ] ...
                        = corr  ( Radio.(radioType).(ivarName).Val(Indx(:,iboot)) ...
                                , Radio.(radioType).(jvarName).Val(Indx(:,iboot)) ...
                                , 'type' , corType );
                    end
                end
            end
        end
        toc
    
    end

    save([Boot.outPath,Boot.filename],'Boot');
    
else

    %if ~exist('Boot','var')
        disp('loading bootstrap Mat file, containing Structure ''Boot'' ...');
        load([Boot.outPath,Boot.filename])
    %else
    %    warning('Variable ''Boot'' already exists in MATLAB environment. Skipping Boot load from Hard Drive...');
    %end
end


% plot histograms and plots
Boot.Var.Label = {'(z+1)','E_{iso} [ ergs ]','T_{90z} [ s ]'};

Boot.Cor.Spearman.ZoneEiso.Hist.textLoc = [ 0.03, 0.77 ];
Boot.Cor.Spearman.ZoneDurz.Hist.textLoc = [ 0.97, 0.77 ];
Boot.Cor.Spearman.EisoDurz.Hist.textLoc = [ 0.03, 0.77 ];
Boot.Cor.Spearman.ZoneEiso.Hist.textAlign = 'left';
Boot.Cor.Spearman.ZoneDurz.Hist.textAlign = 'right';
Boot.Cor.Spearman.EisoDurz.Hist.textAlign = 'left';
Boot.Cor.Spearman.ZoneEiso.Hist.legendLoc = 'northwest';
Boot.Cor.Spearman.ZoneDurz.Hist.legendLoc = 'northeast';
Boot.Cor.Spearman.EisoDurz.Hist.legendLoc = 'northwest';

Boot.Var.LogZone.Hist.Avg.xlabel = 'Sample Mean of Log10 ( ';
Boot.Var.LogZone.Hist.Avg.textLoc = [ 0.97, 0.77 ];
Boot.Var.LogEiso.Hist.Avg.textLoc = [ 0.03, 0.77 ];
Boot.Var.LogDurz.Hist.Avg.textLoc = [ 0.03, 0.77 ];
Boot.Var.LogZone.Hist.Avg.textAlign = 'right';
Boot.Var.LogEiso.Hist.Avg.textAlign = 'left';
Boot.Var.LogDurz.Hist.Avg.textAlign = 'left';
Boot.Var.LogZone.Hist.Avg.legendLoc = 'northeast';
Boot.Var.LogEiso.Hist.Avg.legendLoc = 'northwest';
Boot.Var.LogDurz.Hist.Avg.legendLoc = 'northwest';
Boot.Var.LogZone.Hist.Avg.Limits = [0.40 0.70] * log(10);
Boot.Var.LogEiso.Hist.Avg.Limits = [51.0 53.6] * log(10);
Boot.Var.LogDurz.Hist.Avg.Limits = [0.30 1.70] * log(10);

Boot.Var.LogZone.Hist.Std.xlabel = 'Sample Standard Deviation of Log10 ( ';
Boot.Var.LogZone.Hist.Std.textLoc = [ 0.97, 0.77 ];
Boot.Var.LogEiso.Hist.Std.textLoc = [ 0.03, 0.77 ];
Boot.Var.LogDurz.Hist.Std.textLoc = [ 0.03, 0.77 ];
Boot.Var.LogZone.Hist.Std.textAlign = 'right';
Boot.Var.LogEiso.Hist.Std.textAlign = 'left';
Boot.Var.LogDurz.Hist.Std.textAlign = 'left';
Boot.Var.LogZone.Hist.Std.legendLoc = 'northeast';
Boot.Var.LogEiso.Hist.Std.legendLoc = 'northwest';
Boot.Var.LogDurz.Hist.Std.legendLoc = 'northwest';
Boot.Var.LogZone.Hist.Std.Limits = [0.13 0.35] * log(10);
Boot.Var.LogEiso.Hist.Std.Limits = [0.18 0.74] * log(10);
Boot.Var.LogDurz.Hist.Std.Limits = [0.20 0.65] * log(10);

for icorType = 1:1%1:Boot.Cor.Type.count
    corType = Boot.Cor.Type.Name{icorType};
    for ivar = 1:Boot.Var.count

        ivarName = Boot.Var.Name{ivar};

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        StatName = {'Avg','Std'};

        % plot Avg, Std histograms
        for iStat = 1:length(StatName)

            iStatName = StatName{iStat};

            if varFigExportRequested, close all, figure('visible','off','Color','none'), else, figure, end
                hold on; box on; %colormap('cool');
                LegendName = cell(Radio.Type.count,1);
                for irtype = 1:Radio.Type.count
                    radioType = Radio.Type.Name{irtype};
                    LegendName{irtype} = ['Radio-',radioType];
                    h = histogram( Boot.Var.(ivarName).(radioType).(iStatName) / log(10) ,'normalization','probability');
                    %h.BinWidth = 0.02;
                    h.BinLimits = Boot.Var.(ivarName).Hist.(iStatName).Limits / log(10);
                    h.NumBins = 100;
                    h.EdgeAlpha = 0; % remove bin edges
                    Boot.Var.(ivarName).(radioType).Hist.Height.Values = h.Values;
                    Boot.Var.(ivarName).(radioType).Hist.Height.CumSum = cumsum(h.Values);
                    Boot.Var.(ivarName).(radioType).avg = mean( Boot.Var.(ivarName).(radioType).(iStatName) / log(10) );
                    Boot.Var.(ivarName).(radioType).std =  std( Boot.Var.(ivarName).(radioType).(iStatName) / log(10) );
                    %xlabel([corType,' Correlation Strength: ',corName], 'fontSize', fontSize);
                end
                xlabel([Boot.Var.LogZone.Hist.(iStatName).xlabel, Boot.Var.Label{ivar},' )'], 'Interpreter', 'Tex', 'fontSize', fontSize);
                ylabel('Normalized Bootstrapping Count', 'Interpreter', 'Tex', 'fontSize', fontSize);
                %xlim([-0.82 0.82]);
                legend  ( LegendName ...
                        , 'location' , Boot.Var.(ivarName).Hist.(iStatName).legendLoc ...
                        , 'fontSize' , fontSize ...
                        , 'color' , 'none' ...
                        )
                legend boxoff;
    
                % Compute the odds of Radio-Bright correlations being higher than Radio-Dark correlations
                prob = sum  ( Boot.Var.(ivarName).Bright.Hist.Height.Values ...
                            .*Boot.Var.(ivarName).Dark.Hist.Height.CumSum );
                text( Boot.Var.(ivarName).Hist.(iStatName).textLoc(1) ...
                    , Boot.Var.(ivarName).Hist.(iStatName).textLoc(2) ...
                    , ['\pi ( {Bright} \geq {Dark} ) = ',sprintf('%.2f',round(prob,2))] ...
                    , 'Units', 'normalized', 'HorizontalAlignment', Boot.Var.(ivarName).Hist.(iStatName).textAlign ...
                    , 'Interpreter', 'tex', 'fontSize' , fontSize );
                text( Boot.Var.(ivarName).Hist.(iStatName).textLoc(1) ...
                    , Boot.Var.(ivarName).Hist.(iStatName).textLoc(2) - 0.1 ...
                    , [ '{Bright} = ' ...
                    , sprintf('%.2f',round(Boot.Var.(ivarName).Bright.avg,2)) ...
                    , ' \pm ' ...
                    , sprintf('%.2f',round(Boot.Var.(ivarName).Bright.std,2)) ...
                    ] ...
                    , 'Units', 'normalized', 'HorizontalAlignment', Boot.Var.(ivarName).Hist.(iStatName).textAlign ...
                    , 'Interpreter', 'tex', 'fontSize' , fontSize );
                text( Boot.Var.(ivarName).Hist.(iStatName).textLoc(1) ...
                    , Boot.Var.(ivarName).Hist.(iStatName).textLoc(2) - 0.2 ...
                    , [ '{Dark} = ' ...
                    , sprintf('%.2f',round(Boot.Var.(ivarName).Dark.avg,2)) ...
                    , ' \pm ' ...
                    , sprintf('%.2f',round(Boot.Var.(ivarName).Dark.std,2)) ...
                    ] ...
                    , 'Units', 'normalized', 'HorizontalAlignment', Boot.Var.(ivarName).Hist.(iStatName).textAlign ...
                    , 'Interpreter', 'tex', 'fontSize' , fontSize );
    
            set(gca, 'color', 'none', 'fontsize', fontSize);
            if varFigExportRequested
                fileName = [Boot.outPath,'BootHist_',iStatName,'Log10',ivarName(4:end),'.png'];
                disp(['Generating Figure ',fileName]);
                export_fig(fileName,'-m4 -transparent');
                hold off; close(gcf);
            else
                hold off;
            end

        end % istat

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for jvar = ivar+1:Boot.Var.count

            jvarName = Boot.Var.Name{jvar};
            corName = [ivarName(4:end),jvarName(4:end)];

            % Correlation strength histograms
            if corFigExportRequested, close all, figure('visible','off','Color','none'), else, figure, end
                hold on; box on; %colormap('cool');
                LegendName = cell(Radio.Type.count,1);
                for irtype = 1:Radio.Type.count
                    radioType = Radio.Type.Name{irtype};
                    LegendName{irtype} = ['Radio-',radioType];
                    h = histogram(Boot.Cor.(corType).(corName).(radioType).Sval,'normalization','probability');
                    h.BinWidth = 0.02;
                    h.BinLimits = [-1.,1.];
                    h.NumBins = 100;
                    h.EdgeAlpha = 0; % remove bin edges
                    Boot.Cor.(corType).(corName).(radioType).Hist.Height.Values = h.Values;
                    Boot.Cor.(corType).(corName).(radioType).Hist.Height.CumSum = cumsum(h.Values);
                    Boot.Cor.(corType).(corName).(radioType).avg = mean( Boot.Cor.(corType).(corName).(radioType).Sval);
                    Boot.Cor.(corType).(corName).(radioType).std =  std( Boot.Cor.(corType).(corName).(radioType).Sval);
                    %xlabel([corType,' Correlation Strength: ',corName], 'fontSize', fontSize);
                end
                xlabel([corType,'''s Correlation Coefficient: ',Boot.Var.Label{ivar},' - ',Boot.Var.Label{jvar}], 'Interpreter', 'Tex', 'fontSize', fontSize);
                ylabel('Normalized Bootstrapping Count', 'Interpreter', 'Tex', 'fontSize', fontSize);
                xlim([-0.82 0.82]);
                legend  ( LegendName ...
                        , 'location' , Boot.Cor.(corType).(corName).Hist.legendLoc ...
                        , 'fontSize' , fontSize ...
                        , 'color' , 'none' ...
                        )
                legend boxoff;

                % Compute the odds of Radio-Bright correlations being higher than Radio-Dark correlations
                prob = sum  ( Boot.Cor.(corType).(corName).Bright.Hist.Height.Values ...
                            .*Boot.Cor.(corType).(corName).Dark.Hist.Height.CumSum );
                text( Boot.Cor.(corType).(corName).Hist.textLoc(1) ...
                    , Boot.Cor.(corType).(corName).Hist.textLoc(2) ...
                    , ['\pi ( \rho_{Bright} \geq \rho_{Dark} ) = ',sprintf('%.2f',round(prob,2))] ...
                    , 'Units', 'normalized', 'HorizontalAlignment', Boot.Cor.(corType).(corName).Hist.textAlign ...
                    , 'Interpreter', 'tex', 'fontSize' , fontSize );
                text( Boot.Cor.(corType).(corName).Hist.textLoc(1) ...
                    , Boot.Cor.(corType).(corName).Hist.textLoc(2) - 0.1 ...
                    , [ '\rho_{Bright} = ' ...
                      , sprintf('%.2f',round(Boot.Cor.(corType).(corName).Bright.avg,2)) ...
                      , ' \pm ' ...
                      , sprintf('%.2f',round(Boot.Cor.(corType).(corName).Bright.std,2)) ...
                      ] ...
                    , 'Units', 'normalized', 'HorizontalAlignment', Boot.Cor.(corType).(corName).Hist.textAlign ...
                    , 'Interpreter', 'tex', 'fontSize' , fontSize );
                text( Boot.Cor.(corType).(corName).Hist.textLoc(1) ...
                    , Boot.Cor.(corType).(corName).Hist.textLoc(2) - 0.2 ...
                    , [ '\rho_{Dark} = ' ...
                      , sprintf('%.2f',round(Boot.Cor.(corType).(corName).Dark.avg,2)) ...
                      , ' \pm ' ...
                      , sprintf('%.2f',round(Boot.Cor.(corType).(corName).Dark.std,2)) ...
                      ] ...
                    , 'Units', 'normalized', 'HorizontalAlignment', Boot.Cor.(corType).(corName).Hist.textAlign ...
                    , 'Interpreter', 'tex', 'fontSize' , fontSize );

            set(gca, 'color', 'none', 'fontsize', fontSize);
            if corFigExportRequested
                fileName = [Boot.outPath,'BootHistCor_',ivarName(4:end),jvarName(4:end),'.png'];
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

