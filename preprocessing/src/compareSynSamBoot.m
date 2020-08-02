format compact; format long;
filePath = mfilename('fullpath');
[scriptPath,fileName,fileExt] = fileparts(filePath); cd(scriptPath);
addpath(genpath('../../../../lib/matlab/')) % lib codes

bootstrap
comCorSim
close all
need radio type use radioType = Radio.Type.Name{irtype} where irtype is
an integer of either 1 or 2
Boot.filename = 'CorNicoleBootstrap.mat';


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
            Indx = randi([1,Radio.(radioType).count],Radio.(radioType).count,Nicole.(radioType).count);
            % compute bootstrap correlations
            counter = 0;
            for ivar = 1:Boot.Var.count
                ivarName = Boot.Var.Name{ivar};
                NicoleBoot.Var.(ivarName).(radioType).Avg = zeros(Nicole.(radioType).count,1);
                NicoleBoot.Var.(ivarName).(radioType).Std = zeros(Nicole.(radioType).count,1);
                for iboot = 1:Nicole.(radioType).count
                    NicoleBoot.Var.(ivarName).(radioType).Avg(iboot) = mean( Radio.(radioType).(ivarName).Val(Indx(:,iboot)) );
                    NicoleBoot.Var.(ivarName).(radioType).Std(iboot) =  std( Radio.(radioType).(ivarName).Val(Indx(:,iboot)) );
                end
                for jvar = ivar+1:Boot.Var.count
                    jvarName = Boot.Var.Name{jvar};
                    counter = counter + 1;
                    corName = [ivarName(4:end),jvarName(4:end)];
                    NicoleBoot.Cor.(corType).(corName).(radioType).Sval = zeros(Nicole.(radioType).count,1);
                    NicoleBoot.Cor.(corType).(corName).(radioType).Pval = zeros(Nicole.(radioType).count,1);
                    for iboot = 1:Nicole.(radioType).count
                        [ NicoleBoot.Cor.(corType).(corName).(radioType).Sval(iboot) ...
                        , NicoleBoot.Cor.(corType).(corName).(radioType).Pval(iboot) ] ...
                        = corr  ( Radio.(radioType).(ivarName).Val(Indx(:,iboot)) ...
                                , Radio.(radioType).(jvarName).Val(Indx(:,iboot)) ...
                                , 'type' , corType );
                    end
                end
            end
        end
        toc
    
    end

CorType = {'ZoneEiso','ZoneDurz','EisoDurz'};

QuietET90cor=0.45798;
QuietEZcor=0.33913;
Quietdeviation=0.00145;


QuietET90=Cor.Spearman.EisoDurz;
QuietZT90=Cor.Spearman.ZoneDurz;
QuietZE=Cor.Spearman.ZoneEiso;


QuietET90res=QuietET90(QuietET90>QuietET90cor-Quietdeviation & QuietET90<QuietET90cor+Quietdeviation);
QuietZEres=QuietZE(QuietZE>QuietEZcor-Quietdeviation & QuietZE<QuietEZcor+Quietdeviation);
QuietZT90res=QuietZT90(QuietET90>QuietET90cor-Quietdeviation & QuietET90<QuietET90cor+Quietdeviation & QuietZE>QuietEZcor-Quietdeviation & QuietZE<QuietEZcor+Quietdeviation);





LoudET90cor=0.23352;
LoudEZcor=0.076638;
Louddeviation=0.0115;

LoudET90=Cor.Spearman.EisoDurz;
LoudZT90=Cor.Spearman.ZoneDurz;
LoudZE=Cor.Spearman.ZoneEiso;


LoudET90res=LoudET90(LoudET90>LoudET90cor-Louddeviation & LoudET90<LoudET90cor+Louddeviation);
LoudZEres=LoudZE(LoudZE>LoudEZcor-Louddeviation & LoudZE<LoudEZcor+Louddeviation);
LoudZT90res=LoudZT90(LoudET90>LoudET90cor-Louddeviation & LoudET90<LoudET90cor+Louddeviation & ... 
LoudZE>LoudEZcor-Louddeviation & LoudZE<LoudEZcor+Louddeviation);

Loud={LoudZEres,LoudZT90res,LoudET90res};
Quiet={QuietZEres,QuietZT90res,QuietET90res};
    
    counter=0;
  normType='pdf';
  colorName = "white";
            stats={'Avg','Std'};
            colorPairs={{'Cyan','#7E2F8E'},{'yellow','#D95319'}};
            for iVarName1 = 1:Boot.Var.count      
                varName1 = Boot.Var.Name{iVarName1};
                CorName = CorType{iVarName1};
                    counter = counter+1;
                    varName2 = Boot.Var.Name{ivar};
                    
                    figure("color",colorName); hold on; box on;
                        for irtype = 1:Radio.Type.count
                            radioType = Radio.Type.Name{irtype};
                             if isequal(radioType,'Dark')
                                 histogram(cell2mat(Quiet(iVarName1)),50 ,'normalization',normType...
                                 ,'FaceColor',colorPairs{irtype}{1},'EdgeAlpha',0);
                             else
                                 histogram(cell2mat(Loud(iVarName1)),50 ,'normalization',normType...
                                 ,'FaceColor',colorPairs{irtype}{1},'EdgeAlpha',0);
                             end

                            histogram(Boot.Cor.Spearman.(CorName).(radioType).Sval,'normalization',normType...
                            ,'FaceColor',colorPairs{irtype}{2},'EdgeAlpha',0);

                            xlabel(['Spearman Correlation Coefficient: ',CorName], 'Interpreter', 'Tex', 'fontSize', fontSize);
                            ylabel('Normalized Count', 'Interpreter', 'Tex', 'fontSize', fontSize);
                            legend  ( {'Synthetic Radio-Dark','Boot Radio-Dark','Synthetic Radio-Bright','Boot Radio-Bright'} ...
                                , 'location' , Boot.Var.(ivarName).Hist.(iStatName).legendLoc ...
                                , 'fontSize' , fontSize ...
                                , 'color' , colorName ...
                                )
                        end  
                        set(gca,"color",colorName)
                        fileName = ['../out/SynSamCompareBoot/',CorName,'.png'];
                        export_fig(fileName,'-m4 -transparent');
                    hold off;
                
            end
            


    %save([Boot.outPath,Boot.filename],'Boot');
    