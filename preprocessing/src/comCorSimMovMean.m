%comCorSim % get data
fontSize = 13;

movingWindow = 10000;
[Sorted.LogEiso.Avg,Sorted.LogEiso.Indx] = sort(Cor.Sample.LogEiso.Avg);
Sorted.LogDurz.Avg = Cor.Sample.LogDurz.Avg(Sorted.LogEiso.Indx);
Sorted.CorEisoDurz = Cor.(corType).EisoDurz(Sorted.LogEiso.Indx);
MovMean.LogDurz = movmean(Sorted.LogDurz.Avg,movingWindow);
MovMean.CorEisoDurz = movmean(Sorted.CorEisoDurz,movingWindow);

figure;
plot( Sorted.LogEiso.Avg/log(10) , MovMean.LogDurz/log(10) );
set(gca,'xscale','log');
set(gca,'yscale','log');

figure;
plot( exp(Sorted.LogEiso.Avg) , exp(MovMean.LogDurz) );

figure;
plot( Sorted.LogEiso.Avg/log(10) , MovMean.CorEisoDurz );


figure;
xscale = 'log';
yscale = 'log';
colormap('default');
skip = 10000;
x = exp(Sorted.LogEiso.Avg(1:skip:end));
y = exp(MovMean.LogDurz(1:skip:end));
c = MovMean.CorEisoDurz(1:skip:end);
xx=[x x];           %// create a 2D matrix based on "X" column
yy=[y y];           %// same for Y
zz=zeros(size(xx)); %// everything in the Z=0 plane
cc =[c c] ;         %// matrix for "CData"
s = surf( xx ...
        , yy ...
        , zz ...
        , cc ...
        , 'EdgeColor','interp','FaceColor','none','Marker','.', 'linewidth', 2 );
xlabel(['Mean LogEiso:  E_{iso}'], 'Interpreter', 'tex', 'fontSize', fontSize);
ylabel(['Mean LogDurz:  T_{90z}'], 'Interpreter', 'tex', 'fontSize', fontSize);
xrange = [5e52 2e53];
yrange = [13 32];
xlim(xrange);
ylim(yrange);
%if strcmp(xscale,'log'), xlim(log10(xrange)); else, xlim(xrange); end
%if strcmp(yscale,'log'), ylim(log10(yrange)); else, ylim(yrange); end
set(gca,'color','none','XMinorTick','on','YMinorTick','on', 'fontSize', fontSize);
set(gca,'xscale',xscale);
set(gca,'yscale',yscale);
cbar = colorbar;
cbar.Label.String = ['Spearman''s \rho:  E_{iso} - T_{90z}'];
cbar.Label.FontSize = fontSize;
cbar.Label.Interpreter = 'tex';
shading flat                    %// so each line segment has a plain color
view(2) %// view(0,90)          %// set view in X-Y plane


if ~exist('Nicole','var'), Nicole = readLloydRadioData; end

disp([ 'Mean Eiso Dark: ', num2str( exp(mean(Nicole.Dark.LogEiso)) )]);
disp([ 'Mean Eiso Bright: ', num2str( exp(mean(Nicole.Bright.LogEiso)) )]);
disp([ 'Mean Durz Dark: ', num2str( exp(mean(Nicole.Dark.LogDurz)) )]);
disp([ 'Mean Durz Bright: ', num2str( exp(mean(Nicole.Bright.LogDurz)) )]);
disp([ 'Corr Eiso-Durz Dark: ', num2str( corr( Nicole.Dark.LogEiso , Nicole.Dark.LogDurz , 'type', 'sp' ) )]);
disp([ 'Corr Eiso-Durz Bright: ', num2str( corr( Nicole.Bright.LogEiso , Nicole.Bright.LogDurz , 'type', 'sp' ) )]);
disp([ 'Corr Zone-Eiso Dark: ', num2str( corr( Nicole.Dark.LogZone , Nicole.Dark.LogEiso , 'type', 'sp' ) )]);
disp([ 'Corr Zone-Eiso Bright: ', num2str( corr( Nicole.Bright.LogZone , Nicole.Bright.LogEiso , 'type', 'sp' ) )]);
disp([ 'Corr Zone-Durz Dark: ', num2str( corr( Nicole.Dark.LogZone , Nicole.Dark.LogDurz , 'type', 'sp' ) )]);
disp([ 'Corr Zone-Durz Bright: ', num2str( corr( Nicole.Bright.LogZone , Nicole.Bright.LogDurz , 'type', 'sp' ) )]);

exp( mean( Nicole.Dark.LogDurz.Val(Nicole.Dark.LogEiso.Val>min(Nicole.Bright.LogEiso.Val)) ) )
exp( mean( Nicole.Dark.LogEiso.Val(Nicole.Dark.LogEiso.Val>min(Nicole.Bright.LogEiso.Val)) ) )

% Nicole.Mean.LogEiso.dark = mean(Nicole.Dark.LogEiso);
% Nicole.Mean.Log10Eiso.dark = Nicole.Mean.LogEiso / log(10);
% Nicole.Mean.LogDurz.dark = mean(Nicole.Dark.LogDurz);
% Nicole.Mean.Log10Durz.dark = Nicole.Mean.LogDurz / log(10);
% 
% Nicole.Mean.LogEiso.bright = mean(Nicole.Dark.LogEiso);
% Nicole.Mean.Log10Eiso.bright = Nicole.Mean.LogEiso / log(10);
% Nicole.Mean.LogDurz.bright = mean(Nicole.Dark.LogDurz);
% Nicole.Mean.Log10Durz.bright = Nicole.Mean.LogDurz / log(10);

