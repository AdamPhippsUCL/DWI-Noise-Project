% Script to save DICOM scans as .mat files


% Define DICOM folder
dfolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\Noise-Statistics-Project\New\Phantom Experiments\Imaging Data\DICOM");


% Define MAT folder
matfolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\Noise-Statistics-Project\New\Phantom Experiments\Imaging Data\MAT");


% Read contents of directory
contents = dir(dfolder);
fnames = {contents(3:end).name};


% Loop over fnames and attempt to read as DICOM

for findx = 1:length(fnames)

    fname = [dfolder '/' char(fnames{findx})];

    % Try to parse as DICOM file
    try
        dinfo = dfparse(fname);
        Name = dinfo(1).ProtocolName;
        darray = squeeze(dicomread(fname)); 

        % PV Scaling
        RS = [dinfo.RescaleSlope];
        RS = reshape(RS, [1, size(RS)]);
        RS = repmat(RS, [size(darray, 1), size(darray, 2), 1]);
    
        RI = [dinfo.RescaleIntercept];
        RI = reshape(RI, [1, size(RI)]);
        RI = repmat(RI, [size(darray, 1), size(darray, 2), 1]);
    
        SS = [dinfo.Private_2005_100e];
        SS = reshape(SS, [1, size(SS)]);
        SS = repmat(SS, [size(darray, 1), size(darray, 2), 1]);
    
        % b=0 image
        whereb0 = [dinfo.DiffusionBValue]== 0;
        b0img = (double(darray(:,:, whereb0)).*RS(:,:,whereb0) + RI(:,:,whereb0))./(RS(:,:,whereb0) .* SS(:,:,whereb0));
    
        % b>0 image (trace image)
        whereb = ([dinfo.DiffusionBValue]~= 0).*([dinfo.DiffusionDirectionality]== 2) == 1;
        bimg = (double(darray(:,:, whereb)).*RS(:,:,whereb) + RI(:,:,whereb))./(RS(:,:,whereb) .* SS(:,:,whereb));
    
        % Stack Images
        ImageArray = cat(4, b0img, bimg);
    
        % Array of b values (Matching last two dimensions of image array)
        bMatrix =  squeeze(cat(3, [dinfo(whereb0).DiffusionBValue], [dinfo(whereb).DiffusionBValue]));
    
        % Create directory
        folder = [matfolder '/' Name];
        mkdir(folder)
    
        % Save .mat files
        save([folder '/ImageArray.mat'], 'ImageArray');
        save([folder '/bMatrix.mat'], 'bMatrix');

        disp(['Saved data for ' Name])
        disp('')

    catch
        disp([fnames{findx} ' not a DICOM file'])
        continue
    end


end
