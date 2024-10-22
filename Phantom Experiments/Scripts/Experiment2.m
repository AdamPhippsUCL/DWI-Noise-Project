% Script to perform analysis in experiment 1

%% Initial Settings


% Define necessary folders
ImageDatafolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\Noise-Statistics-Project\New\Phantom Experiments\Imaging Data");
imagefolder = [ImageDatafolder '/MAT'];
ROIfolder = [ImageDatafolder '/ROIs'];

% Load image and ROI details
load([ImageDatafolder '/ImageDetails.mat'])
load([ImageDatafolder '/ROIDetails.mat'])

% Define image names to use
ImageNames = {
    "b2500_Avg1_Ex3",...
    "b2500_Avg2_Ex3",...
    "b2500_Avg4_Ex3",...
    "b2500_Avg6_Ex3",...
    "b2500_Avg8_Ex3",...  
    };

Nimg = length(ImageNames);

% Define ROI names to use
ROINames = {
    "ROI1",...
    % "ROI2",...
    % "ROI3",...
    % "ROI4",...
    % "ROI5",...
    % "ROI6",...
    % "ROI7",...
    % "ROI8"
    };

NROI = length(ROINames);

% Initial structure for results
FittingResults = struct();


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


        % Load image
        img = double(load([imagefolder '/' char(ImageName) '/ImageArray.mat']).ImageArray);


        % Normalised image
        normimg = img(:,:,:,2)./img(:,:,:,1);

        % Remove nans and infinities
        normimg(isnan(normimg)) = 0;
        normimg(isinf(normimg)) = 0;

        % % Check proper loading
        % figure;
        % imshow(normimg(:,:,6), [])
        % figure;
        % imshow( double( ROI(:,:,6)), [])


        % Extract ROI values
        ROIvalues = normimg(ROI);


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

        % Max 100, min 200
        nbin =  min( [round((binmax-binmin)/binspacing), 1000] );
        if nbin<200
            nbin=200;
        end
        close(ftemp)


        
        binedges = linspace(binmin, binmax, nbin+1);
        bincentres = (binedges(:,1:end-1) + binedges(:,2:end)) / 2 ;
        binspacing = bincentres(2)-bincentres(1);

        f=figure;
        h = histogram(ROIvalues, binedges);
        % binedges = h.BinEdges;
        % binspacing = binedges(2)-binedges(1);
        % bincentres = (1/2)*(binedges(1:end-1)+binedges(2:end));
        counts = h.Values;
        xlim([0 2])
        ylim([0 1.05*max(counts)]);


        % == Fitting

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
        close(f);

        % Fill in results
        FittingResults((ROIindx-1)*Nimg + imgindx).ImageName = ImageName;
        FittingResults((ROIindx-1)*Nimg + imgindx).Rav = Rav;
        FittingResults((ROIindx-1)*Nimg + imgindx).ROIName = ROIName;        
        FittingResults((ROIindx-1)*Nimg + imgindx).sigma0fit = sigma0fit;
        FittingResults((ROIindx-1)*Nimg + imgindx).T2 = T2;
        FittingResults((ROIindx-1)*Nimg + imgindx).T2fit = T2fit;
        FittingResults((ROIindx-1)*Nimg + imgindx).D = D;        
        FittingResults((ROIindx-1)*Nimg + imgindx).Dfit = Dfit;
        FittingResults((ROIindx-1)*Nimg + imgindx).resnorm = resnorm;

    end
end


%% Plot results

figfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\Noise-Statistics-Project\New\Phantom Experiments\Figures";

colordict = dictionary([1,2,3,4,5,6,7,8], [	"#0072BD", 	"#D95319",	"#EDB120", 	"#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F", "#e981cd"]);
concentrations = ["50% PVP", "40% PVP", "30% PVP", "20% PVP", "10% PVP", "5% PVP", "2.5% PVP", "1mM NiCl_2"];

fig1 = figure;

sigma0s = zeros(NROI, Nimg);

% Scatter plot for each ROI
for ROIindx = 1:NROI

    % Get ROI details
    ROIName = ROINames{ROIindx};
    whereROI = [ROIDetails.ROIName{:}] == ROIName;
    T2 = ROIDetails.T2{whereROI};
    D = ROIDetails.D{whereROI};

    % ROI bools
    ROIbools = ([FittingResults.ROIName] == ROIName);

    % Extract sigma0 measurements
    sigma0s(ROIindx, :) = [FittingResults(ROIbools).sigma0fit];
    Ravs = [FittingResults(ROIbools).Rav];

    plot(Ravs, sigma0s, '*-', DisplayName = ['ROI ' num2str(ROIindx)], Color = colordict(ROIindx))
    hold on
    xlabel('N_{b>0}')
    ylabel('\sigma_0')
    ylim([0,0.1])
    legend(NumColumns = 2)
    ax = gca();
    ax.FontSize = 12;

end
% 
% figure;
% 
% % Scatter plot for each ROI
% for ROIindx = opts.ROInums
% 
%     trueT2 = ROIinfo(ROIindx).T2;
%     bools = [FittingResults.ROIindx] == ROIindx;
%     T2s = [FittingResults(bools).T2fit];
%     Ravs = [ImageSettings(:).Nav_ratio];
% 
%     plot(Ravs, T2s, '*-', DisplayName = ['ROI ' num2str(ROIindx)], Color = colordict(ROIindx))
%     hold on
%     plot(Ravs, trueT2*ones(size(Ravs)), '--', Color = colordict(ROIindx))
%     xlabel('N_{b>0}')
%     ylabel('T2 (ms)')
%     ylim([0,1000])
%     xlim([1, 25])
%     ax = gca();
%     ax.FontSize = 12;       
%     % legend
% end
% 
% 
% figure;
% 
% % Scatter plot for each ROI
% for ROIindx = opts.ROInums
% 
%     trueADC = ROIinfo(ROIindx).ADC;
%     bools = [FittingResults.ROIindx] == ROIindx;
%     ADCs = [FittingResults(bools).ADCfit];
%     Ravs = [ImageSettings(:).Nav_ratio];
% 
%     plot(Ravs, ADCs, '*-', DisplayName = ['ROI ' num2str(ROIindx)], Color = colordict(ROIindx))
%     hold on
%     plot(Ravs, trueADC*ones(size(Ravs)), '--', Color = colordict(ROIindx))
%     % plot(Ravs, )
%     xlabel('N_{b>0}')
%     ylabel('D (mm^2/s)')
%     ylim([0,1.2e-3])
%     xlim([1, 25])
%     ax = gca();
%     ax.FontSize = 12;       
%     % legend
% end