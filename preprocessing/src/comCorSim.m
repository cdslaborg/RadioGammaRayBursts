% read synthetic data
plotSyntheticSample
Nicole = readLloydRadioData;
close all

simRunRequested = 0;
figExportRequested = 1;

corType = 'Spearman';
%corType = 'Kendall';

% generate random samples
SynSam.Mask = log(ZModel.B10.Synthetic.data(:,3)) > 1.197344248356904e+02;
SynSam.LogEiso = log( ZModel.B10.Synthetic.data( SynSam.Mask , 3 ) );
SynSam.LogDurz = log( ZModel.B10.Synthetic.data( SynSam.Mask , 4 ) );
SynSam.LogZone = log( ZModel.B10.Synthetic.data( SynSam.Mask , 9 ) + 1.0 );
SynSam.DetProb = ZModel.B10.Synthetic.data( SynSam.Mask , 14 ); %josh changed 10 to 14
SynSam.size = length(SynSam.LogEiso);

%figure; hold on; box on; plot( SynSam.LogEiso , SynSam.LogDurz , '.' ); set(gca,'color','none'); hold off;
%figure; hold on; box on; plot( SynSam.LogZone , SynSam.LogDurz , '.' ); set(gca,'color','none'); hold off;
%figure; hold on; box on; plot( SynSam.LogZone , SynSam.LogEiso , '.' ); set(gca,'color','none'); hold off;

% repeatable random numbers
rng(0,'twister');

corFileName = [outPathSynSam,'comCorSimFixedCutoff.mat'];
if simRunRequested
    Cor.Sample.size = 50;
    Cor.nsim = 20000000;
    Cor.(corType).EisoDurz = zeros(Cor.nsim,1);
    Cor.(corType).ZoneDurz = zeros(Cor.nsim,1);
    Cor.(corType).ZoneEiso = zeros(Cor.nsim,1);
    Cor.Sample.LogEiso.Avg = zeros(Cor.nsim,1);
    Cor.Sample.LogDurz.Avg = zeros(Cor.nsim,1);
    Cor.Sample.LogZone.Avg = zeros(Cor.nsim,1);
    Cor.Sample.LogEiso.std = zeros(Cor.nsim,1);
    Cor.Sample.LogDurz.std = zeros(Cor.nsim,1);
    Cor.Sample.LogZone.std = zeros(Cor.nsim,1);
    %Pearson.Cor = zeros(Cor.nsim,1);
    for isim = 1:Cor.nsim
        if mod(isim,1000)==0, disp(['isim = ', num2str(isim)]), end
        isample = 1;
        while isample<=Cor.Sample.size
            %Indx = randi([1 SynSam.size],Cor.Sample.size,1);
            Indx(isample) = randi([1 SynSam.size]);
            if SynSam.DetProb(Indx(isample)) > rand()
                isample = isample + 1;
            end
        end
        Dum.LogEiso = SynSam.LogEiso(Indx);
        Dum.LogDurz = SynSam.LogDurz(Indx);
        Dum.LogZone = SynSam.LogZone(Indx);
        Cor.(corType).EisoDurz(isim) = corr( Dum.LogEiso , Dum.LogDurz , 'type' , corType );
        Cor.(corType).ZoneDurz(isim) = corr( Dum.LogZone , Dum.LogDurz , 'type' , corType );
        Cor.(corType).ZoneEiso(isim) = corr( Dum.LogZone , Dum.LogEiso , 'type' , corType );
        % compute sample statistics
        Cor.Sample.LogEiso.Avg(isim) = mean( Dum.LogEiso );
        Cor.Sample.LogDurz.Avg(isim) = mean( Dum.LogDurz );
        Cor.Sample.LogZone.Avg(isim) = mean( Dum.LogZone );
        Cor.Sample.LogEiso.Std(isim) = std( Dum.LogEiso );
        Cor.Sample.LogDurz.Std(isim) = std( Dum.LogDurz );
        Cor.Sample.LogZone.Std(isim) = std( Dum.LogZone );
    end
else
    if exist('Cor','var')
        warning('Variable Cor already exists in MATLAB. Skipping the reload...');
    else
        disp('Loading 1.1Gb Cor object from file...');
        tic
        load(corFileName); % loads struct Cor
        toc
    end
end


Cor.Name = {'EisoDurz','ZoneEiso','ZoneDurz'};
Cor.Label = {'E_{iso} - T_{90z}','(z+1) - E_{iso}','(z+1) - T_{90z}'};
Cor.Triple = [ 1,2,3 ; 2,3,1 ; 1,3,2 ]';
close all;
skip = Cor.nsim / 1000000;
for itriple = 1:length(Cor.Triple(1,:))
    xid = Cor.Triple(1,itriple);
    yid = Cor.Triple(2,itriple);
    zid = Cor.Triple(3,itriple);
    plotName = ['comCorSimFixedCutoff_',Cor.Name{xid},'_',Cor.Name{yid},'_',Cor.Name{zid}];
    disp(['Plotting ', plotName]);
    tic
    if figExportRequested, figure('visible','off','Color','none'), else, figure, end, hold on, box on;
        scatter ( Cor.(corType).(Cor.Name{xid})(1:skip:end) ...
                , Cor.(corType).(Cor.Name{yid})(1:skip:end) ...
                , 5*ones(length(Cor.(corType).(Cor.Name{yid})(1:skip:end)),1) ...
                , Cor.(corType).(Cor.Name{zid})(1:skip:end) ...
                ,'filled' );
        plot( Nicole.Dark.Cor.Spearman.(Cor.Name{xid}).coef ...
            , Nicole.Dark.Cor.Spearman.(Cor.Name{yid}).coef ...
            , '.', 'markersize', 20, 'color', 'black' );
        plot( Nicole.Bright.Cor.Spearman.(Cor.Name{xid}).coef ...
            , Nicole.Bright.Cor.Spearman.(Cor.Name{yid}).coef ...
            , '.', 'markersize', 20, 'color', 'red' );
        plots = flip(findall(gcf,'Type','Line'));
        leg = legend(plots,{'Quiet','Loud'}, 'Location', 'northwest'); %added by josh
        set(leg,'Box','off')
        xlabel(['Spearman''s \rho:  ',Cor.Label{xid}], 'Interpreter', 'tex', 'fontSize', fontSize);
        ylabel(['Spearman''s \rho:  ',Cor.Label{yid}], 'Interpreter', 'tex', 'fontSize', fontSize);
        cbar = colorbar;
        cbar.Label.String = ['Spearman''s \rho:  ',Cor.Label{zid}];
        cbar.Label.FontSize = fontSize;
        set(gca,'color','none','XMinorTick','on','YMinorTick','on', 'fontSize', fontSize);
    if figExportRequested
        fileName = [outPathSynSam,plotName,'.png'];
        export_fig (fileName,'-m4 -transparent');
        hold off; close(gcf);
    else
        hold off;
    end
    toc
end


tic
if figExportRequested, figure('visible','off','Color','none'), else, figure, end, hold on, box on;
    plotName = 'comCorSimFixedCutoff_MeanEiso_CorEisoDurz_Durz';
    disp(['Plotting ', plotName]);
    scatter ( (Cor.Sample.LogEiso.Avg(1:skip:end))/log(10) ...
            , Cor.(corType).EisoDurz(1:skip:end) ...
            , 5*ones(length(Cor.(corType).EisoDurz(1:skip:end)),1) ...
            , (Cor.Sample.LogDurz.Avg(1:skip:end))/log(10) ...
            , 'filled',  'displayname', '' );
	plot( Nicole.Dark.LogEiso.avg/log(10) ...
        , Nicole.Dark.Cor.Spearman.EisoDurz.coef ...
        , '.', 'markersize', 20, 'color', 'black', 'displayname', 'Quiet');
    plot( Nicole.Bright.LogEiso.avg/log(10) ...
        , Nicole.Bright.Cor.Spearman.EisoDurz.coef ...
        , '.', 'markersize', 20, 'color', 'red', 'displayname', 'Loud' );
    plots = flip(findall(gcf,'Type','Line'));
    leg = legend(plots,{'Quiet','Loud'}, 'Location', 'northwest'); %added by josh
    set(leg,'Box','off')
    xlabel(['Mean LogEiso:  E_{iso}'], 'Interpreter', 'tex', 'fontSize', fontSize);
    ylabel(['Spearman''s \rho:  E_{iso} - T_{90z}'], 'Interpreter', 'tex', 'fontSize', fontSize);
    cbar = colorbar;
    cbar.Label.String = ['Mean LogDurz:  T_{90z}'];
    cbar.Label.FontSize = fontSize;
    set(gca,'color','none','XMinorTick','on','YMinorTick','on', 'fontSize', fontSize);
if figExportRequested
    fileName = [outPathSynSam,plotName,'.png'];
    export_fig (fileName,'-m4 -transparent');
    hold off; close(gcf);
else
    hold off;
end
toc


tic
if figExportRequested, figure('visible','off','Color','none'), else, figure, end, hold on, box on;
    plotName = 'comCorSimFixedCutoff_MeanEiso_Durz_CorEisoDurz';
    disp(['Plotting ', plotName]);
    scatter ( (Cor.Sample.LogEiso.Avg(1:skip:end))/log(10) ...
            , (Cor.Sample.LogDurz.Avg(1:skip:end))/log(10) ...
            , 5*ones(length(Cor.(corType).EisoDurz(1:skip:end)),1) ...
            , Cor.(corType).EisoDurz(1:skip:end) ...
            , 'filled' );
	plot( Nicole.Dark.LogEiso.avg/log(10) ...
        , Nicole.Dark.LogDurz.avg/log(10) ...
        , '.', 'markersize', 20, 'color', 'black', 'displayname', 'Quiet' );
    plot( Nicole.Bright.LogEiso.avg/log(10) ...
        , Nicole.Bright.LogDurz.avg/log(10) ...
        , '.', 'markersize', 20, 'color', 'red', 'displayname', 'Loud' );
    leg = legend; %added by josh
    set(leg,'Box','off')
    xlabel(['Mean LogEiso:  E_{iso}'], 'Interpreter', 'tex', 'fontSize', fontSize);
    ylabel(['Mean LogDurz:  T_{90z}'], 'Interpreter', 'tex', 'fontSize', fontSize);
    cbar = colorbar;
    cbar.Label.String = ['Spearman''s \rho:  E_{iso} - T_{90z}'];
    cbar.Label.FontSize = fontSize;
    set(gca,'color','none','XMinorTick','on','YMinorTick','on', 'fontSize', fontSize);
if figExportRequested
    fileName = [outPathSynSam,plotName,'.png'];
    export_fig (fileName,'-m4 -transparent');
    hold off; close(gcf);
else
    hold off;
end
toc


%figure; hold on; box on; plot(SynSam.LogZone,SynSam.LogEiso, '.'); set(gca,'color','none'); hold off;
if exist(corFileName,'file')
    warning(['Variable Cor has been already exported to external file: ', corFileName, '. Skipping the export...']);
else
    save(corFileName,'Cor');
end


if ~exist('Nicole','var'), Nicole = readLloydRadioData; end

RadioType = {'Dark','Bright'};
VarType = {'Zone','Eiso','Durz'};
CorType = {'ZoneEiso','ZoneDurz','EisoDurz'};
for ivar = 1:length(VarType)
    var = VarType{ivar};
    for irtype = 1:length(RadioType)
        rtype = RadioType{irtype};
        disp([ 'Mean ', var, ' ', rtype, ': ', num2str(Nicole.(rtype).(['Log',var]).avg), ' ', num2str(num2str(Nicole.(rtype).(['Log10',var]).avg)), ' ', num2str(Nicole.(rtype).(var).avg) ]);
    end
end
for icor = 1:length(CorType)
    cor = CorType{icor};
    for irtype = 1:length(RadioType)
        rtype = RadioType{irtype};
        disp([ 'Corr ', cor, ' ', rtype, ': ', num2str( Nicole.(rtype).Cor.Spearman.(cor).coef ) ]);
    end
end

