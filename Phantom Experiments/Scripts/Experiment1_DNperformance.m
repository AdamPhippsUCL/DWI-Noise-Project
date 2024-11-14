% Test to see how sigma0 estimate changes after denoising

Results = load("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Code\DWI-Noise-Project\Phantom Experiments\Outputs\MAT Ratio (8)\FittingResults.mat").FittingResults;
ResultsDN = load("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Code\DWI-Noise-Project\Phantom Experiments\Outputs\MAT DN Ratio (8)\FittingResults.mat").FittingResults;

sigma0changes = ([ResultsDN(:).sigma0fit] - [Results(:).sigma0fit])./[Results(:).sigma0fit] ;


ImageNames = {Results(:).ImageName};
ROINames = {Results(:).ROIName};
T2s = [Results(:).T2];
uniqueT2s = unique(T2s);

NROI = length(unique(T2s));
Nimg = length(unique(ImageNames));


% Load image details
load("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Code\DWI-Noise-Project\Phantom Experiments\Imaging Data\ImageDetails.mat");





%% Plot changes

colordict = dictionary([1,2,3,4,5,6,7,8], [	"#0072BD", 	"#D95319",	"#EDB120", 	"#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F", "#e981cd"]);
concentrations = dictionary(["ROI1","ROI2","ROI3","ROI4","ROI5","ROI6","ROI7","ROI8"], ["50%", "40%", "30%", "20%", "10%", "5%", "2.5%", "1mM"] );

Sigma0Changes = zeros(NROI, Nimg);

figure;

for ROIindx = 1:NROI

    ROIbools = (T2s == uniqueT2s(ROIindx));
    
    Sigma0Changes(ROIindx, :) = sigma0changes(ROIbools);
    scatter(ROIindx*ones(Nimg, 1), sigma0changes(ROIbools), '*', MarkerEdgeColor = colordict(ROIindx), DisplayName = concentrations(['ROI' num2str(ROIindx)])); 
    hold on

end

boxplot(transpose(Sigma0Changes))
ylim([-1 0])


% Get TE and b values
imgnames = {Results(ROIbools).ImageName};
TEs = zeros(size(imgnames));
bvals = zeros(size(imgnames));
for imgindx = 1:length(imgnames)
    imgname = imgnames{imgindx};
    bool = ([ImageDetails.ProtocolName{:}] == imgname);
    TEs(imgindx) = ImageDetails.TE{bool};
    bvals(imgindx) = ImageDetails.bvalue{bool};
end

TEs(1) = 451;
[~,ITE]=sort(TEs);
[~,Ib]=sort(bvals);


f1=figure;
for ROIindx = 1:NROI
    plot(bvals(Ib), Sigma0Changes(ROIindx,Ib), '--*', color= colordict(ROIindx), DisplayName = concentrations(['ROI' num2str(ROIindx)]));
    hold on
end

boxplot(Sigma0Changes(:,Ib), Positions=bvals(Ib), Widths=80)
ylim([-0.8 -0.2])
xlim([1250 3750])
xticks([1500 2000 2500 3000 3500])
xticklabels([1500 2000 2500 3000 3500])
legend
xlabel('b value (s/mm^2)')
ylabel('Fractional Change in \sigma_0')


saveas(f1, [char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Code\DWI-Noise-Project\Phantom Experiments\Outputs") '/sigma0 DN change (bval).fig']);


f2=figure;
for ROIindx = 1:NROI
    plot(TEs(ITE), Sigma0Changes(ROIindx,ITE), '--*', color= colordict(ROIindx), DisplayName = concentrations(['ROI' num2str(ROIindx)]));
    hold on
end

boxplot(Sigma0Changes(:,ITE), Positions=TEs(ITE), Widths=16)
ylim([-0.8 -0.2])
xlim([50 500])
xticks([100 200 300 400 500])
xticklabels([100 200 300 400 500])
xlabel('TE (ms)')
ylabel('Fractional Change in \sigma_0')
% legend

saveas(f2, [char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Code\DWI-Noise-Project\Phantom Experiments\Outputs") '/sigma0 DN change (TE).fig']);
