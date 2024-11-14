% Script to test new bias removal method


% == Original images 

% Folder path for MAT files
MATfolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\DWI-Noise-Project\Phantom Experiments\Imaging Data\MAT");

% Read contents
contents = dir(MATfolder);

% folder names
foldernames = {contents(3:end).name};


% First image
findx = 11;

% Original Image
folder = [MATfolder '/' foldernames{findx}];
img = load([folder '/ImageArray.mat']).ImageArray;
img = double(img);


%% Apply denoising

window = [3 3];

imgDN = zeros(size(img));
S2 = zeros(size(img, 1:3));

% Denoise slicewise
Nslice = size(img, 3);
for slindx = 1:Nslice
    
    slice = squeeze(img(:,:,slindx,:) );

    % Apply denoising
    [sliceDN, sliceS2, P] = denoise(slice, window);

    imgDN(:,:,slindx,:) = sliceDN;
    S2(:,:,slindx) = sliceS2;

end

figure;
imshow(imgDN(:,:,5,2)./imgDN(:,:,5,1), [])


%% Setting up functions

% Load calibration curves for this NSA
calfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\DWI-Noise-Project\Functions\Calibration Curves";
NSA = 3;

load(fullfile(calfolder, ['NSA ' num2str(NSA)], 'MeanCalibration.mat'))
load(fullfile(calfolder, ['NSA ' num2str(NSA)], 'SigmaCalibration.mat'))

SNRs = keys(MeanDict);
alphas = values(MeanDict);
betas = values(SigmaDict);


xprimes = linspace(0.1,5,101);
Line3 = zeros(size(xprimes));
for indx = 1:length(xprimes)
    xprime = xprimes(indx);
    Line3(indx) = eqn3(xprime);
end


% Test patch;
row = 50;
col = 75;
width = 3;
slice = 5;

imgpatch = img(row:row+width, col:col+width, slice, 2);
imgDNpatch = imgDN(row:row+width, col:col+width, slice, 2);
S2patch = S2(row:row+width, col:col+width, slice);

figure;
tiledlayout(1,3);

nexttile;
imshow(imgpatch, []);
colorbar;

nexttile;
imshow(imgDNpatch, []);
colorbar;

nexttile;
imshow(S2patch, []);
colorbar;


muhat = mean(imgDNpatch(:));
sigmahat = sqrt(mean(S2patch(:)));

yhat = sigmahat/muhat;

% Find xprime
[~,I] = min( abs( yhat-Line3));
xprime = 1/SNRs(I);




%% Defining functions


function sigmaprime = eqn1(muhat, xprime)

    sigmaprime = muhat/( sqrt(pi/2)*Lhalf( -1/(2*xprime^2) ));


end



function yhat = eqn3(xprime)

    yhatsqrd = (2/pi)*( 1/( Lhalf( -1/(2*xprime^2) )^2) )*(2 + 1/(xprime^2)) - 1;

    yhat = sqrt(abs(yhatsqrd));

end


function out = Lhalf(in)
        
    out = exp(in/2)*( (1-in)*besselj(0,-in/2) - in*besselj(1,-in/2) );

end



% % Load calibration curves for this NSA
% calfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\DWI-Noise-Project\Functions\Calibration Curves";
% NSA = 3;
% 
% load(fullfile(calfolder, ['NSA ' num2str(NSA)], 'MeanCalibration.mat'))
% load(fullfile(calfolder, ['NSA ' num2str(NSA)], 'SigmaCalibration.mat'))
% 
% SNRs = keys(MeanDict);
% alphas = values(MeanDict);
% betas = values(SigmaDict);
% 
% figure;
% plot(SNRs, alphas)
% hold on
% plot(SNRs, betas)
% 
% BetaAlphaRatios = betas./alphas;
% 
% 
% Line = BetaAlphaRatios.*SNRs;
% % 
% % figure;
% % plot(SNRs, Line, '*')
% 
% % Solve for sigma and mu
% 
% [~,I] = min(abs(yprime-Line));
% SNR = SNRs(I);
% 
% 
% % Patch scalings
% patchalpha = MeanDict(SNR)
% patchbeta = SigmaDict(SNR)







