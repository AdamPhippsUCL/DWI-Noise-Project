% Script to apply MP-PCA denoising + bias removal

% Imaging Data folder
ImagingDataFolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Code\DWI-Noise-Project\Phantom Experiments\Imaging Data");

% Folder path for original and denoised MAT files
MATfolder = fullfile(ImagingDataFolder, 'MAT');
Denoisedfolder = fullfile(ImagingDataFolder, 'MAT DN (window 5)');

% Read contents
contents = dir(MATfolder);

% folder names
foldernames = {contents(3:end).name};


% Denoising window
window = [5 5];

% For bias removal
sigma = 1;
SNRs = linspace(0,10,201);
SNRhats = zeros(size(SNRs));
for indx = 1:length(SNRs)
    SNRhats(indx) = eqn3(SNRs(indx));
end


for findx = 1:length(foldernames)

    imgname = foldernames{findx}

    folder = [MATfolder '/' foldernames{findx}];

    % Load image
    img = load([folder '/ImageArray.mat']).ImageArray;
    img = double(img);

    % Load b matrix
    bMatrix = load([folder '/bMatrix.mat']).bMatrix;

    % Load image details
    load(fullfile(ImagingDataFolder, 'ImageDetails.mat'));

    imgbools = [ImageDetails.ProtocolName{:}]==foldernames{findx};

    if all(imgbools==0)
        continue
    end

    NSA = ImageDetails.NSA{imgbools};
    Rav = ImageDetails.Rav{imgbools};
    Nb = NSA*Rav; 
    
    % Load alpha beta calibration curves for b>0 image
    try
        calfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Code\DWI-Noise-Project\Functions\Calibration Curves";
        load(fullfile(calfolder, ['NSA ' num2str(Nb)], 'MeanCalibration.mat'))
        load(fullfile(calfolder, ['NSA ' num2str(Nb)], 'SigmaCalibration.mat'))
    catch
        disp(['Calibration curves not generated for NSA = ' num2str(Nb)])
        continue
    end

    dictSNRs = keys(MeanDict);
    alphas = values(MeanDict);
    betas = values(SigmaDict);
    AlphaBetaRatios = betas./alphas;
    AlphaBetaLine = AlphaBetaRatios.*dictSNRs;

    % Image dimensions
    Nrow = size(img, 1);
    Ncol = size(img, 2);
    Nslice = size(img, 3);

    % Initialise empty array for denoised image
    imgDN = zeros(size(img));

    % Initial empty array for signal bisa (alpha)
    Rescale = zeros(size(img(:,:,:,1)));

    % Denoise slicewise
    for slindx = 1:Nslice

        slindx
        
        % Select slice
        slice = squeeze(img(:,:,slindx,:) );

        % Apply denoising
        [sliceDN, S2, P] = denoise(slice, window);

        imgDN(:,:,slindx,:) = sliceDN;
    

        %% Adding bias removal (experimental)

        patchwidth = 3;

        for indx = round((patchwidth-1)/2)+1:Nrow-round((patchwidth-1)/2)
            for jndx = round((patchwidth-1)/2)+1:Ncol-round((patchwidth-1)/2)

                imgpatch = slice(indx-round((patchwidth-1)/2):indx+round((patchwidth-1)/2), jndx-round((patchwidth-1)/2):jndx+round((patchwidth-1)/2), 2);
                imgDNpatch = sliceDN(indx-round((patchwidth-1)/2):indx+round((patchwidth-1)/2), jndx-round((patchwidth-1)/2):jndx+round((patchwidth-1)/2), 2);
                S2patch = S2( indx-round((patchwidth-1)/2):indx+round((patchwidth-1)/2), jndx-round((patchwidth-1)/2):jndx+round((patchwidth-1)/2) );

                patchmuhat = mean(imgpatch(:));
                patchDNmuhat = mean(imgDNpatch(:));
                patchsigmahat = sqrt(mean(S2patch(:)));
                patchSNRhat = patchmuhat/patchsigmahat;

                [~,I] = min( abs( patchSNRhat-SNRhats));
                patchSNRprime = SNRs(I);

                patchsigmaprime = inveqn1(patchmuhat, patchSNRprime);
                patchmuprime = patchSNRprime*patchsigmaprime;

                [~,I] = min( abs( patchSNRprime-AlphaBetaLine.*dictSNRs));
                SNR = dictSNRs(I);

                patchalpha = alphas(I);

                rescale = (1/patchalpha)*(patchmuprime/patchDNmuhat); 

                % Catch errors
                if and(rescale <= 1, rescale ~=0)
                    Rescale(indx, jndx, slindx) = rescale;
                else
                    Rescale(indx, jndx, slindx) = 1;
                end

            end

        end




    end


    % Save denoised image
    DNfolder = [Denoisedfolder '/' foldernames{findx}];
    mkdir(DNfolder)
    ImageArray = imgDN;
    save([DNfolder '/ImageArray.mat'], 'ImageArray')
    save([DNfolder '/Rescale.mat'], 'Rescale')
    save([DNfolder '/bMatrix.mat'], 'bMatrix')

% 
% sl = 4;
% figure;
% imshow(img(:,:,sl, 1), [])
% colorbar;
% figure;
% imshow(img(:,:,sl, 2), [])
% colorbar;
% figure;
% imshow(imgDN(:,:,sl,1), [])
% colorbar;
% figure;
% imshow(imgDN(:,:,sl,2), [])
% colorbar;


end



%% Define functions

function out = Lhalf(in)
    out = exp(in/2)*( (1-in)*besseli(0,-in/2) - in*besseli(1,-in/2)  );
end

function out = muhat(nu, sigma)
    out = sigma*sqrt(pi/2)*Lhalf( - (nu^2)/(2*(sigma^2)) );
end


function out = sigmahat(nu, sigma)
    outsqrd = 2*(sigma^2) + nu^2 - (muhat(nu,sigma))^2;
    out = sqrt( outsqrd  );
end


function out = eqn3(SNR)
    nu = SNR;
    sigma=1;
    out = muhat(nu, sigma)/sigmahat(nu, sigma);
end


function out = inveqn1(mu, SNR)
    out = mu/( sqrt(pi/2)*Lhalf( -0.5*(SNR^2) )  );
end
