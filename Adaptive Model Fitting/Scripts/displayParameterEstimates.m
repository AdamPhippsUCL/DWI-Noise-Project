
%% Load parameter estimates

outputfolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Code\DWI-Noise-Project\Adaptive Model Fitting\Outputs\Parameter Estimates");

patnum = '20240628_Patient10';

geom = 'Original';

sigma0 = load([outputfolder '/' patnum '/' geom '/sigma0.mat']).sigma0;
T2 = load([outputfolder '/' patnum '/' geom '/T2.mat']).T2;
b0fromhighb = load([outputfolder '/' patnum '/' geom '/b0fromhighb.mat']).b0fromhighb;
META = load([outputfolder '/' patnum '/' geom '/META.mat']).META;

%% Display figures

slindx = 8;

xmin=50;
xmax=xmin+55;

ymin=50;
ymax=ymin+55;

f=figure;
f.Position = [300   200   409   300];
tiledlayout(1,1, "TileSpacing","compact");
nexttile;
imshow(b0fromhighb(ymin:ymax,xmin:xmax, slindx), []);
c=colorbar;
c.Label.String = 'b=0 signal';
ax = gca();
ax.FontSize = 12;
% title('b=0')
saveas(f, [outputfolder '/' patnum '/' geom '/b0fromhighb.fig'])

f=figure;
f.Position = [300   200   411   300];
tiledlayout(1,1, "TileSpacing","compact");
nexttile;
imshow(sigma0(ymin:ymax,xmin:xmax, slindx), [0 0.1]);
c=colorbar;
c.Label.String='\sigma_0';
ax = gca();
ax.FontSize = 12;
% title('Estimated \sigma_0' )
saveas(f, [outputfolder '/' patnum '/' geom '/sigma0.fig'])

f=figure;
f.Position = [300   200   401   300];
tiledlayout(1,1, "TileSpacing","compact");
nexttile;
imshow(T2(ymin:ymax,xmin:xmax, slindx), [0 200]);
c=colorbar;
c.Label.String = 'T2 (ms)';
ax = gca();
ax.FontSize = 12;
% title('Estimated T2')
saveas(f, [outputfolder '/' patnum '/' geom '/T2.fig'])


% Add to Meta data
META.xs = [xmin, xmax];
META.ys = [ymin, ymax];
META.slice = slindx;
save([outputfolder '/' patnum '/' geom '/META.mat'], "META");
