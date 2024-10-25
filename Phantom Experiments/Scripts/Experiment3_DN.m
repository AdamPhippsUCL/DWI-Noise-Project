% Experiment applying DN to images and see how histograms change



%% Initial Settings


% Define necessary folders
ImageDatafolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\DWI-Noise-Project\Phantom Experiments\Imaging Data");
imagefolder = [ImageDatafolder '/MAT'];
imagefolderDN = [ImageDatafolder '/Denoised'];
ROIfolder = [ImageDatafolder '/ROIs'];

% Load image and ROI details
load([ImageDatafolder '/ImageDetails.mat'])
load([ImageDatafolder '/ROIDetails.mat'])

% Define image names to use
ImageNames = {
    'b1500_Ex1',...
    'b2000_Ex1',...
    'b2100_Ex1',...
    'b2400_Ex1',...
    'b2500_Ex1',...
    'b2750_Ex1',...
    'b3000_Ex1',...
    'b3500_Ex1',...
    'b4000_Ex1'...
    };

Nimg = length(ImageNames);

% Define ROI names to use
ROINames = {
    % "ROI1",...
    % "ROI2",...
    % "ROI3",...
    % "ROI4",...
    "ROI5",...
    % "ROI6",...
    % "ROI7",...
    % "ROI8"
    };

NROI = length(ROINames);

% Initial structure for results
FittingResults = struct();
FittingResultsDN = struct();

%% Run Experiment

% Define noise model type
noisetype = 'Ratio';

for ROIindx = 1:NROI

    % Get ROI details
    ROIName = ROINames{ROIindx};
    whereROI = [ROIDetails.ROIName{:}] == ROIName;
    T2 = ROIDetails.T2{whereROI};
    D = ROIDetails.D{whereROI};



    for imgindx = 1:Nimg

        % Get image details
        ImageName = ImageNames{imgindx};
        whereimg = [ImageDetails.ProtocolName{:}] == ImageName;
        b = ImageDetails.bvalue{whereimg};
        TE = ImageDetails.TE{whereimg};
        TR = ImageDetails.TR{whereimg};
        NSA = ImageDetails.NSA{whereimg};
        Rav = ImageDetails.Rav{whereimg};
        

        % Load ROI
        ROIstruct = load([ROIfolder '/' char(ImageName) '/' char(ROIName) '.mat']);
        for name = fieldnames(ROIstruct)
            ROI = getfield(ROIstruct, name{1});
            ROI = (ROI==1);
        end


        % Load image and denoised image
        img = double(load([imagefolder '/' char(ImageName) '/ImageArray.mat']).ImageArray);
        imgDN = double(load([imagefolderDN '/' char(ImageName) '/ImageArray.mat']).ImageArray);

        % Normalised images
        normimg = img(:,:,:,2)./img(:,:,:,1);
        normimgDN = imgDN(:,:,:,2)./imgDN(:,:,:,1);

        % Remove nans and infinities
        normimg(isnan(normimg)) = 0;
        normimg(isinf(normimg)) = 0;
        normimgDN(isnan(normimgDN)) = 0;
        normimgDN(isinf(normimgDN)) = 0;

        % Extract ROI values
        ROIvalues = normimg(ROI);
        ROIvaluesDN = normimgDN(ROI);



        % Define bin edges and centers
        binmin = 0;
        binmax = 2;

        if max(ROIvalues)<1
            binmax = 1;
        end

        % Temporary histogram to get optimal bin spacing
        ftemp=figure;
        h = histogram(ROIvalues, 100);
        binedges = h.BinEdges;
        binspacing = binedges(2)-binedges(1);

        % Max 1200, min 200
        nbin =  min( [round((binmax-binmin)/binspacing), 1200] );
        if nbin<200
            nbin=200;
        end
        close(ftemp)

        
        binedges = linspace(binmin, binmax, nbin+1);
        bincentres = (binedges(:,1:end-1) + binedges(:,2:end)) / 2 ;
        binspacing = bincentres(2)-bincentres(1);


        % == Fitting to original image

        f=figure;
        h = histogram(ROIvalues, binedges);
        % binedges = h.BinEdges;
        % binspacing = binedges(2)-binedges(1);
        % bincentres = (1/2)*(binedges(1:end-1)+binedges(2:end));
        counts = h.Values;
        xlim([0 2])
        ylim([0 1.05*max(counts)]);

        % fitting

        % First guess
        beta0guess = [0.02, T2, exp(-D*b)];

        % Bounds
        lb = [0.001, T2-50, exp(-(D+0.5e-4)*b)];
        ub = [0.4, T2+50, exp(-(D-0.5e-4)*b)];

        % Apply fitting
        [coeffs, resnorm] = fitDistToHist( ...
            counts, ...
            bincentres, ...
            disttype=noisetype,...
            TE=TE, ...
            N0 = NSA, ...
            Nb = NSA*Rav,...
            beta0guess=beta0guess,...
            lb = lb, ...
            ub = ub);
        
        % Fitted parameters
        sigma0fit = coeffs(1);
        T2fit= coeffs(2);
        fdfit = coeffs(3);

        Dfit = (-1/b)*log(fdfit);


        % == Plot fitted distribution
        b0signal = exp(-TE/T2fit);
        bsignal = fdfit*b0signal;

        switch noisetype
            case 'Ratio'
                [dist, signals] = RatioDistRician(b0signal, bsignal, sigma0fit, N0 = NSA, Nb = NSA*Rav, zs = linspace(0,2,200));
            case 'Rice'
                [dist, signals] = RiceDist(b0signal, bsignal, sigma0fit, zs = linspace(0,2,200));
        end
        
        hold on

        xlabel('Normalized signal')
        ylabel('Counts')
        switch noisetype
            case 'Ratio'
                plot(signals, dist*binspacing*sum(counts), LineWidth = 2, DisplayName = [noisetype ' distribution'], color = "#D95319");
            case 'Rice'
                plot(signals, dist*binspacing*sum(counts), LineWidth = 2, DisplayName = [noisetype ' distribution'], color = "#7E2F8E");
        end

        ax = gca();
        ax.FontSize = 12;
        legend;
        pause(1)



        % == Fitting to denoised image

        fDN=figure;
        h = histogram(ROIvaluesDN, binedges);
        % binedges = h.BinEdges;
        % binspacing = binedges(2)-binedges(1);
        % bincentres = (1/2)*(binedges(1:end-1)+binedges(2:end));
        countsDN = h.Values;
        xlim([0 2])
        ylim([0 1.05*max(countsDN)]);

        % fitting

        % First guess
        beta0guess = [0.02, T2, exp(-D*b)];

        % Bounds
        lb = [0.001, T2-50, exp(-(D+0.5e-4)*b)];
        ub = [0.4, T2+50, exp(-(D-0.5e-4)*b)];

        % Apply fitting
        [coeffs, resnormDN] = fitDistToHist( ...
            countsDN, ...
            bincentres, ...
            disttype=noisetype,...
            TE=TE, ...
            N0 = NSA, ...
            Nb = NSA*Rav,...
            beta0guess=beta0guess,...
            lb = lb, ...
            ub = ub);
        
        % Fitted parameters
        sigma0fitDN = coeffs(1);
        T2fitDN= coeffs(2);
        fdfitDN = coeffs(3);

        DfitDN = (-1/b)*log(fdfitDN);


        % == Plot fitted distribution
        b0signal = exp(-TE/T2fitDN);
        bsignal = fdfitDN*b0signal;

        switch noisetype
            case 'Ratio'
                [dist, signals] = RatioDistRician(b0signal, bsignal, sigma0fitDN, N0 = NSA, Nb = NSA*Rav, zs = linspace(0,2,200));
            case 'Rice'
                [dist, signals] = RiceDist(b0signal, bsignal, sigma0fitDN, zs = linspace(0,2,200));
        end
        
        hold on

        xlabel('Normalized signal')
        ylabel('Counts')
        switch noisetype
            case 'Ratio'
                plot(signals, dist*binspacing*sum(counts), LineWidth = 2, DisplayName = [noisetype ' distribution'], color = "#D95319");
            case 'Rice'
                plot(signals, dist*binspacing*sum(counts), LineWidth = 2, DisplayName = [noisetype ' distribution'], color = "#7E2F8E");
        end

        ax = gca();
        ax.FontSize = 12;
        legend;
        pause(1)

        close(f);
        close(fDN)

        % Fill in results
        FittingResults((ROIindx-1)*Nimg + imgindx).ImageName = ImageName;
        FittingResults((ROIindx-1)*Nimg + imgindx).ROIName = ROIName;        
        FittingResults((ROIindx-1)*Nimg + imgindx).sigma0fit = sigma0fit;
        FittingResults((ROIindx-1)*Nimg + imgindx).T2 = T2;
        FittingResults((ROIindx-1)*Nimg + imgindx).T2fit = T2fit;
        FittingResults((ROIindx-1)*Nimg + imgindx).D = D;        
        FittingResults((ROIindx-1)*Nimg + imgindx).Dfit = Dfit;
        FittingResults((ROIindx-1)*Nimg + imgindx).resnorm = resnorm;


         % Fill in results
        FittingResultsDN((ROIindx-1)*Nimg + imgindx).ImageName = ImageName;
        FittingResultsDN((ROIindx-1)*Nimg + imgindx).ROIName = ROIName;        
        FittingResultsDN((ROIindx-1)*Nimg + imgindx).sigma0fit = sigma0fitDN;
        FittingResultsDN((ROIindx-1)*Nimg + imgindx).T2 = T2;
        FittingResultsDN((ROIindx-1)*Nimg + imgindx).T2fit = T2fitDN;
        FittingResultsDN((ROIindx-1)*Nimg + imgindx).D = D;        
        FittingResultsDN((ROIindx-1)*Nimg + imgindx).Dfit = DfitDN;
        FittingResultsDN((ROIindx-1)*Nimg + imgindx).resnorm = resnormDN;       
    end
end



% Sort images by echo time
TEs = ImageDetails.TE;
theseTEs = [TEs{1:9}];
[~,I] = sort(theseTEs);

figure;
scatter(theseTEs(I), [FittingResultsDN(I).sigma0fit])
hold on
scatter(theseTEs(I), [FittingResults(I).sigma0fit])
ylim([0 0.05])

% figure
% for imgindx = 1:Nimg
%     imgname = ImageNames{imgindx};
%     scatter(linspace(1,NROI,NROI), [FittingResults(string({FittingResults.ImageName})==imgname).sigma0fit], '*', MarkerEdgeColor='b')
%     hold on
%     scatter(linspace(1,NROI,NROI), [FittingResultsDN(string({FittingResultsDN.ImageName})==imgname).sigma0fit], '*', MarkerEdgeColor='r')
% end
% 
% figure
% for imgindx = 1:Nimg
%     imgname = ImageNames{imgindx};
%     scatter(linspace(1,NROI,NROI), [FittingResultsDN(string({FittingResultsDN.ImageName})==imgname).sigma0fit]./[FittingResults(string({FittingResults.ImageName})==imgname).sigma0fit])
%     hold on
% end
% 
% figure
% for imgindx = 1:Nimg
%     scatter(linspace(1,NROI,NROI), [FittingResults(:).T2fit], '*', MarkerEdgeColor='b')
%     hold on
%     scatter(linspace(1,NROI,NROI), [FittingResultsDN(:).T2fit], '*', MarkerEdgeColor='r')
% end
% 
% figure
% for imgindx = 1:Nimg
%     scatter(linspace(1,NROI,NROI), [FittingResults(:).Dfit], '*', MarkerEdgeColor='b')
%     hold on
%     scatter(linspace(1,NROI,NROI), [FittingResultsDN(:).Dfit], '*', MarkerEdgeColor='r')
% end
% 

% %% Plot results
% 
% figfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\Noise-Statistics-Project\New\Phantom Experiments\Figures";
% 
% colordict = dictionary([1,2,3,4,5,6,7,8], [	"#0072BD", 	"#D95319",	"#EDB120", 	"#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F", "#e981cd"]);
% concentrations = ["50% PVP", "40% PVP", "30% PVP", "20% PVP", "10% PVP", "5% PVP", "2.5% PVP", "1mM NiCl_2"];
% 
% 
% % == Figure 1: sigma0 estimation
% 
% fig1 = figure;
% 
% sigma0s = zeros(NROI, Nimg);
% 
% for ROIindx = 1:NROI
% 
%     % Get ROI details
%     ROIName = ROINames{ROIindx};
%     whereROI = [ROIDetails.ROIName{:}] == ROIName;
%     T2 = ROIDetails.T2{whereROI};
%     D = ROIDetails.D{whereROI};
% 
% 
%     % ROI bools
%     ROIbools = ([FittingResults.ROIName] == ROIName);
% 
%     % Extract sigma0 measurements
%     sigma0s(ROIindx, :) = [FittingResults(ROIbools).sigma0fit];
%     scatter(ROIindx*ones(Nimg, 1), sigma0s(ROIindx, :), '*', MarkerEdgeColor = colordict(ROIindx), DisplayName = concentrations(ROIindx)) 
%     hold on
% 
% end
% 
% boxplot(transpose(sigma0s))
% xticks(linspace(1,NROI,NROI))
% xticklabels(ROINames)
% ylabel('Estimated \sigma_0')
% % xlabel('Concentration')
% ylim([0, 0.25]);%opts.sigma0max])
% grid on
% title(noisetype)
% legend;
% ax = gca();
% ax.FontSize = 12;
% % saveas(fig1, [char(figfolder) '/sigma0 (' noisetype ').fig'])
% 
% 
% % === Figure 2: T2 estimation
% 
% fig2 = figure;
% 
% T2s = zeros(NROI, Nimg);
% 
% for ROIindx = 1:NROI
% 
%     % Get ROI details
%     ROIName = ROINames{ROIindx};
%     whereROI = [ROIDetails.ROIName{:}] == ROIName;
%     T2 = ROIDetails.T2{whereROI};
%     D = ROIDetails.D{whereROI};
% 
%     % ROI bools
%     ROIbools = ([FittingResults.ROIName] == ROIName);
% 
%     trueT2 = T2;
%     T2s(ROIindx, : ) = ([FittingResults(ROIbools).T2fit]);
%     scatter(ROIindx*ones(Nimg, 1), T2s(ROIindx, :), '*', MarkerEdgeColor = colordict(ROIindx), DisplayName = concentrations(ROIindx)) 
%     hold on
%     scatter([ROIindx+0.4], [trueT2], MarkerEdgeColor = 'black', MarkerFaceColor='black', Marker = 'o', HandleVisibility = 'off')
% 
% end
% 
% boxplot(transpose(T2s))
% xticks(linspace(1,NROI,NROI))
% xticklabels(ROINames)
% ylabel('Estimated T2')
% % xlabel('Concentration')
% ylim([0,1000])
% grid on
% title(noisetype)
% % legend;
% ax = gca();
% ax.FontSize = 12;
% % saveas(fig2, [char(figfolder) '/T2 (' noisetype ').fig'])
% 
% 
% % === Figure 3: D estimation
% 
% fig3 = figure;
% 
% Ds = zeros(NROI, Nimg);
% 
% for ROIindx = 1:NROI
% 
%     % Get ROI details
%     ROIName = ROINames{ROIindx};
%     whereROI = [ROIDetails.ROIName{:}] == ROIName;
%     T2 = ROIDetails.T2{whereROI};
%     D = ROIDetails.D{whereROI};
% 
%     % ROI bools
%     ROIbools = ([FittingResults.ROIName] == ROIName);
% 
%     trueD = D;
%     Ds(ROIindx, : ) = ([FittingResults(ROIbools).Dfit]);
%     scatter(ROIindx*ones(Nimg, 1), Ds(ROIindx, :), '*', MarkerEdgeColor = colordict(ROIindx), DisplayName = concentrations(ROIindx)) 
%     hold on
%     scatter([ROIindx+0.4], [trueD], MarkerEdgeColor = 'black', MarkerFaceColor='black', Marker = 'o', HandleVisibility = 'off')
% 
% end
% 
% boxplot(transpose(Ds))
% xticks(linspace(1,NROI,NROI))
% xticklabels(ROINames)
% ylabel('Estimated D')
% % xlabel('Concentration')
% ylim([0,1.2e-3])
% grid on
% title(noisetype)
% % legend;
% ax = gca();
% ax.FontSize = 12;
% % saveas(fig3, [char(figfolder) '/D (' noisetype ').fig'])