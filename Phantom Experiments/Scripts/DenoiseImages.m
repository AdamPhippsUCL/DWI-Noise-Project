% Script to apply MP-PCA denoising 

% Folder path for MAT files
MATfolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\DWI-Noise-Project\Phantom Experiments\Imaging Data\MAT");

Denoisedfolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\DWI-Noise-Project\Phantom Experiments\Imaging Data\MAT DN");

% Read contents
contents = dir(MATfolder);

% folder names
foldernames = {contents(3:end).name};


% Denoise each image
window = [3 3];

for findx = 1:length(foldernames)

    foldernames{findx}

    folder = [MATfolder '/' foldernames{findx}];

    % Load image
    img = load([folder '/ImageArray.mat']).ImageArray;
    img = double(img);

    % Load b matrix
    bMatrix = load([folder '/bMatrix.mat']).bMatrix;


    imgDN = zeros(size(img));

    % Denoise slicewise
    Nslice = size(img, 3);
    for slindx = 1:Nslice
        
        slice = squeeze(img(:,:,slindx,:) );

        % Apply denoising
        [sliceDN, S2, P] = denoise(slice, window);

        imgDN(:,:,slindx,:) = sliceDN;

    end


    % Save denoised image
    DNfolder = [Denoisedfolder '/' foldernames{findx}];
    mkdir(DNfolder)
    ImageArray = imgDN;
    save([DNfolder '/ImageArray.mat'], 'ImageArray')
    save([DNfolder '/bMatrix.mat'], 'bMatrix')


sl = 4;
figure;
imshow(img(:,:,sl, 1), [])
colorbar;
figure;
imshow(img(:,:,sl, 2), [])
colorbar;
figure;
imshow(imgDN(:,:,sl,1), [])
colorbar;
figure;
imshow(imgDN(:,:,sl,2), [])
colorbar;


end
