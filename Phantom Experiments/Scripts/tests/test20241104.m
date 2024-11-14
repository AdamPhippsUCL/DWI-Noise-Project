%% Trying to reproduce curves seen in Yakovleva et al.

sigma = 1;

SNRs = linspace(0,10,201);


muhats = zeros(size(SNRs));
sigmahats = zeros(size(SNRs));
SNRhats = zeros(size(SNRs));
for indx = 1:length(SNRs)
    muhats(indx) = muhat(SNRs(indx), sigma);
    sigmahats(indx) = sigmahat(SNRs(indx), sigma);
    SNRhats(indx) = eqn3(SNRs(indx));
end


f1=figure;
plot(SNRs, muhats);
hold on
plot(SNRs, sigmahats);

f2=figure;
plot(SNRs, SNRhats);
ax2=gca();

%%%%%%%%%%%%%%%%%%%%%% BIAS REMOVAL TEST

%% Load image

% == Original image

% Folder path for MAT files
MATfolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\DWI-Noise-Project\Phantom Experiments\Imaging Data\MAT");

% Read contents
contents = dir(MATfolder);

% folder names
foldernames = {contents(3:end).name};

% Image index
findx = 16;

% Original Image
folder = [MATfolder '/' foldernames{findx}];
img = load([folder '/ImageArray.mat']).ImageArray;
img = double(img);

%% Apply denoising

window = [5 5];

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

f3=figure;
imshow(imgDN(:,:,5,2)./sqrt(S2(:,:,5)), [0 10])


%% Testing on patch

% Test patch;
row = 50;
col = 75;
width = 3;
slice = 5;

imgpatch = img(row:row+width, col:col+width, slice, 2);
imgDNpatch = imgDN(row:row+width, col:col+width, slice, 2);
S2patch = S2(row:row+width, col:col+width, slice);

% Find measured mu, sigma, SNR
patchmuhat = mean(imgDNpatch(:));
patchsigmahat = mean(sqrt(S2patch(:)));
patchSNRhat = patchmuhat/patchsigmahat;


% Find SNR (prime)

yline(ax2,patchSNRhat);
[~,I] = min( abs( patchSNRhat-SNRhats));
SNRprime = SNRs(I);
xline(ax2,SNRprime);

sigmaprime = inveqn1(patchmuhat, SNRprime);
muprime = SNRprime*sigmaprime;



% Load alpha and beta 

% Load calibration curves for this NSA
calfolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\DWI-Noise-Project\Functions\Calibration Curves";
NSA = 3;

load(fullfile(calfolder, ['NSA ' num2str(NSA)], 'MeanCalibration.mat'))
load(fullfile(calfolder, ['NSA ' num2str(NSA)], 'SigmaCalibration.mat'))

theseSNRs = keys(MeanDict);
alphas = values(MeanDict);
betas = values(SigmaDict);

AlphaBetaLine = alphas./betas;

% Plot alpha\beta line
figure;
ax3 = gca();
plot(theseSNRs, AlphaBetaLine.*theseSNRs, '-*');


yline(ax3, SNRprime);
[~,I] = min( abs( SNRprime-AlphaBetaLine.*theseSNRs));
SNR = theseSNRs(I);
xline(ax3,SNR);


% Final results for alpha and beta
alpha = alphas(I) % Rescale image by this amount!
beta = betas(I);





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