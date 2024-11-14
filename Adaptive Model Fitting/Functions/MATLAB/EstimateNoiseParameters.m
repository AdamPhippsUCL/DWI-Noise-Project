% MATLAB Function to estimate Ratio noise model parameters


function [sigma0, T2, resnorm] = EstimateNoiseParameters(Echo1fname, Echo2fname, opts)

arguments
    Echo1fname % DICOM fname of first echo image (TE1)
    Echo2fname % DICOM fname of first echo image (TE1)



    % OPTIONS
    
    % registration
    opts.register
    opts.maskhw
    opts.b0fromhighb


    % image info
    opts.multiframe = true

    % Python script for normalising and saving diffusion images
    opts.pyscript = "C:\Users\adam\OneDrive - University College London\UCL PhD\Image-Processing\DICOM\Python\normaliseDiffusionImage.py";

    % Fitting
    opts.patchsize = 5
    opts.disttype = 'DTIRatio'
    opts.nbin = 20;
    opts.mask = [];
    opts.smoothsigma = 1 %for Gaussian smoothing

    % Folder to save mat images in
    opts.save = false
    opts.OutputFolder = ""

end

%% STEPS

% 1. Read DICOM information for each echo image (Need their TE values)
% 2. Save/Load DICOM images as MAT files (CAREFUL WITH IMAGE ORIENTATION!)
% 2. Find ratio of echo images
% 3. Define mask/ROI for estimation
% 4. For each voxel, estimate sigma0 and T2 over surrounding patch (starting values from quick
% median calculation)
% 5. Return maps of sigma0 and T2


opts.ImageOutputFolder = [char(opts.OutputFolder)  '/Images'];

%% 1. Read DICOM information for each echo image (Need their TE values)

dinfo1 = dfparse(Echo1fname);
dinfo2 = dfparse(Echo2fname);

try
    TE1 = dinfo1(1).EffectiveEchoTime;
catch
    TE1 = dinfo1(1).EchoTime;
end

try
    TE2 = dinfo2(1).EffectiveEchoTime;
catch
    TE2 = dinfo2(1).EchoTime;
end

function [b, a] = swap(a, b)
end

% Make sure img2 has longer echo time
if TE2<TE1
    swap(dinfo1, dinfo2)
    swap(Echo1fname, Echo2fname)
    swap(TE1, TE2)
end



SeriesDescription1 = dinfo1(1).ProtocolName;
Ndirec1 = max([dinfo1.DiffGradOrientIdentifier]);
bval1 = dinfo1(1).DiffusionBValue;


SeriesDescription2 = dinfo2(1).ProtocolName;
Ndirec2 = max([dinfo2.DiffGradOrientIdentifier]);
bval2 = dinfo2(1).DiffusionBValue;
try
    TE2 = dinfo2(1).EffectiveEchoTime;
catch
    TE2 = dinfo2(1).EchoTime;
end


if Ndirec1 ~= Ndirec2
    error("Images have different number of DTI directions!");
else
    Ndirec = Ndirec1;
end



%% 2. Save/Load DICOM images as MAT files

RS1 = dinfo1(1).RescaleSlope;
RI1 = dinfo1(1).RescaleIntercept;
IMG1 = RS1*double(squeeze(dicomread(Echo1fname))) + RI1;

RS2 = dinfo1(1).RescaleSlope;
RI2 = dinfo1(1).RescaleIntercept;
IMG2 = RS2*double(squeeze(dicomread(Echo2fname))) + RI2;

% Separate images with diffusion direction
IMG1 = IMG1(:,:,[dinfo1.DiffusionDirectionality]==1);
IMG2 = IMG2(:,:,[dinfo2.DiffusionDirectionality]==1);

%% Register 'average' images to b0fromhighb


if opts.register
    
    % == Find average images
    for direcIndx = 1:Ndirec
        img1 = IMG1(:,:,[dinfo1(:).DiffGradOrientIdentifier] == direcIndx);
        img2 = IMG2(:,:,[dinfo2(:).DiffGradOrientIdentifier] == direcIndx);
    
        if direcIndx == 1
            avgIMG1 = img1;
            avgIMG2 = img2;
        else
            avgIMG1=avgIMG1+img1;
            avgIMG2=avgIMG2+img2;
        end
    end
    
    % normalise
    avgIMG1=mat2gray(avgIMG1, [0 double(prctile(avgIMG1(:),98))]);
    avgIMG2=mat2gray(avgIMG2, [0 double(prctile(avgIMG2(:),98))]);
    b0fromhighb=mat2gray(opts.b0fromhighb, [0 double(prctile(opts.b0fromhighb(:),98))]);
    
    
    % == Find regitstration transforms between slices
    
    [ny, nx, nz] = size(b0fromhighb) ;
    
    % % mask
    % maskhw=opts.maskhw;
    % maskcentre = [ceil( (ny+1)/2 )  ceil((nx+1)/2) ] ;
    % maskc = { max(1,maskcentre(1)-maskhw) : min(ny,maskcentre(1)+maskhw) , ...
    %     max(1,maskcentre(2)-maskhw) : min(nx,maskcentre(2)+maskhw) } ;
    % reg_slice = ceil(nz+1)/2;
    
    
    % % fixed and moving
    % fixed = b0fromhighb(maskc{1}, maskc{2}, reg_slice);
    % moving1 = avgIMG1(maskc{1}, maskc{2}, reg_slice);
    % moving2 = avgIMG2(maskc{1}, maskc{2}, reg_slice);
    
    
    for sliceIndx=1:nz
    
        fixed = b0fromhighb(:,:, sliceIndx);
        moving1 = avgIMG1(:,:, sliceIndx);
        moving2 = avgIMG2(:,:, sliceIndx);
    
        [optimizer,metric] = imregconfig('monomodal');
        tform1 = imregtform(moving1,fixed,'translation',optimizer,metric);
        tform2 = imregtform(moving2,fixed,'translation',optimizer,metric);
    
        for direcIndx = 1:Ndirec
            IMG1(:,:,Ndirec*(sliceIndx-1)+direcIndx) = imwarp(IMG1(:,:,Ndirec*(sliceIndx-1)+direcIndx),tform1,"OutputView",imref2d(size( b0fromhighb(:,:,1) ))) ;
            IMG2(:,:,Ndirec*(sliceIndx-1)+direcIndx) = imwarp(IMG2(:,:,Ndirec*(sliceIndx-1)+direcIndx),tform2,"OutputView",imref2d(size( b0fromhighb(:,:,1) ))) ;
        end
    end

end


%% Create stacked normalised image

for direcIndx1 = 1:Ndirec
        
    img1 = IMG1(:,:,[dinfo1(:).DiffGradOrientIdentifier] == direcIndx1);
    
    for direcIndx2 = 1:Ndirec

        img2 = IMG2(:,:,[dinfo2(:).DiffGradOrientIdentifier] == direcIndx2);

        % Remove zeros
        img1(img1==0) = eps;
        img2(img2==0) = eps;
    
    
        % Calculate ratio
        ratioimg = double(img2)./double(img1);
    
    
        % Configure array for ratio image stack
        if and(direcIndx1 == 1, direcIndx2 == 1)
            RatioImageStack = zeros([size(ratioimg) Ndirec^2]);
        end
    
        RatioImageStack(:,:,:, (direcIndx1-1)*Ndirec+direcIndx2 ) = ratioimg;

    end
end

clear IMG1 IMG2 


%% PARAMETER ESTIMATION

% Define blank arrays for T2 and sigma0 estimates and error
T2 = squeeze( zeros(size(RatioImageStack(:,:,:,1))) );
sigma0 = squeeze( zeros(size(RatioImageStack(:,:,:,1))) );
error = squeeze( zeros(size(RatioImageStack(:,:,:,1))) );

%% Define mask (for testing)

% Define slice and pixel indices in small ROI in capsule
zs = [1:size(ratioimg, 1)];
ys = [opts.patchsize:size(ratioimg, 2)-opts.patchsize];
xs = [opts.patchsize:size(ratioimg, 3)-opts.patchsize];

% Patch width
pw = ceil((opts.patchsize-1)/2);



%% Define histogram
disp('Calibrating...')
for zindx = zs
    for yindx = ys
        for xindx = xs
            
            patchvals = RatioImageStack(zindx, yindx-pw:yindx+pw, xindx-pw:xindx+pw, :);
            patchvals  = patchvals(:);
            
            % N Define outlier range
            iqr = (prctile(patchvals, 75) - prctile(patchvals, 25));
            up = median(patchvals) + 1.5*iqr;
            down = median(patchvals) - 1.5*iqr;
            
            % Remove outliers
            outlierbools = or(patchvals>up, patchvals<down);
            patchvals(outlierbools) = [];



            % == INITIAL GUESSES

            % T2
            T2guess = -(TE2-TE1)/log(mean(patchvals));

            % Deal with unrealistic values
            T2guess(T2guess<25) = 25;
            T2guess(T2guess>250) = 250;


            % sigma0
            rTE = (2*(TE1^2))/(TE1+TE2);
            sigma0guess = std(patchvals)/(sqrt(2)*exp(rTE/T2guess));
            sigma0guess(sigma0guess > 0.25) = 0.25;            
            sigma0guess(sigma0guess < 0.001) = 0.001;


           
            % % == USING RATIO DISTRIBUTION (SLOW)

            % % Define bin edges and centers
            % binmin = 0;
            % binmax = 2;
            % nbin = opts.nbin;
            % binedges = linspace(binmin, binmax, nbin+1);
            % bincenters = (binedges(:,1:end-1) + binedges(:,2:end)) / 2 ;
            % binspacing = bincenters(2)-bincenters(1);
            % 
            % 
            % counts = histcounts(patchvals, binedges);

            % [coeffs, resnorm] = fitDistToHist( ...
            %     counts, ...
            %     bincenters, ...
            %     fd = 1, ...
            %     TE = [TE1,TE2], ...
            %     disttype = opts.disttype, ...
            %     beta0guess = [sigma0guess, T2guess,1], ...
            %     boundstype = 'fixed',...
            %     lb = [0.001, 30, 0.1],...
            %     ub = [0.25, 400, 1]...
            %     );
            % 
            % sigma0fit = coeffs(1)
            % T2fit = coeffs(2)
            

            % == JUST USING INITIAL GUESSES (fast)
 
            sigma0fit = sigma0guess;
            T2fit = T2guess;
            resnorm = 1;


            
            % Append to arrays
            T2(zindx, yindx, xindx) = T2fit;
            sigma0(zindx, yindx, xindx) = sigma0fit;
            error(zindx, yindx, xindx) = resnorm;
            
            % % 
            % % Define histogram
            % % f = figure('visible','off');
            % f= figure;
            % H = histogram(patchvals, binedges);
            % hold on;
            % counts = H.Values;
            % 
            % % Generate fitted distribution
            % [dist, signals] = RatioDistRician(exp(-TE1/T2fit), exp(-TE2/T2fit), sigma0fit, Nav_ratio=1);
            % 
            % % Add to histogram
            % hold on
            % plot(signals, dist*sum(counts)*binspacing, '-*')
            % % xlim([0.5,1.5])
            % % T2fit
            % % sigma0fit
            % close(f);

        end
    end
end


%% Smoothing

smoothsigma = opts.smoothsigma;
% % Filter slices
for slice = zs
    T2(slice,:,:) = imgaussfilt(squeeze( T2(slice,:,:)), smoothsigma);
    sigma0(slice,:,:) = imgaussfilt(squeeze( sigma0(slice,:,:)), smoothsigma);
end

%% Save maps

if opts.save

    % Save Output images as mat files
    MapOutputFolder = char(opts.OutputFolder); 
    save([MapOutputFolder '/T2.mat'], "T2")
    save([MapOutputFolder '/sigma0.mat'], "sigma0")
    % save([MapOutputFolder '/DTI1.mat'], "img1")
    % save([MapOutputFolder '/DTI2.mat'], "img2")

end
% % % Filter slices
% slice = 5;
% % T2slice = squeeze(T2(slice,:,:));
% % sigma0slice = squeeze(sigma0(slice, :,:));
% T2slice = imgaussfilt(squeeze( T2(slice,:,:)), 1);
% sigma0slice = imgaussfilt(squeeze( sigma0(slice,:,:)), 1);
% 
% % figure;
% imshow(squeeze(img1(slice,:,:)), [])
% colorbar
% figure;
% imshow(T2slice, [0, 200])
% colorbar
% figure;
% imshow(sigma0slice*0.5, [0.0 0.075])
% colorbar

end