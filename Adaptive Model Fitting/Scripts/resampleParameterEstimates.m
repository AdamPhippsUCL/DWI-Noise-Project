% SCRIPT to resample Ratio noise model parameter estimate maps

%% Load original parameter estimates

outputfolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Code\DWI-Noise-Project\Adaptive Model Fitting\Outputs\Parameter Estimates\20240912_Patient21");

% Original geometry
geomin = 'shV3';

sigma0 = load([outputfolder '/' geomin '/sigma0.mat']).sigma0;
T2 = load([outputfolder  '/' geomin '/T2.mat']).T2;
b0fromhighb = double(load([outputfolder  '/' geomin '/b0fromhighb.mat']).b0fromhighb);
META = load([outputfolder  '/' geomin '/META.mat']).META;

% Example DICOM with same image geometry
dfname1 = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Imaging Data\Patient Volunteers\20240912_Patient21\shV3\Series 401 [MR - b1800shV3]\1.3.6.1.4.1.5962.99.1.3883698028.1870441860.1726165616492.8.0.dcm";
% Get coordinates of DICOM image
COORDS1 = constructVoxelCoordinates(dfname1);


%% Resample parameter maps

% New image geometry
geomout = 'Original';

% Example DICOM with new image geometry
dfname2 = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Imaging Data\Patient Volunteers\20240912_Patient21\Original\Series 801 [MR - SWITCH DB TO YES b3000 80]\0.dcm";

% Get coordinates of new DICOM image
COORDS2 = constructVoxelCoordinates(dfname2);

% Resample each parameter map
sigma0 = resample(sigma0, COORDS1, COORDS2);
T2 = resample(T2, COORDS1, COORDS2);
b0fromhighb = resample(b0fromhighb, COORDS1, COORDS2);

% Save new estimates
mkdir([outputfolder '/' geomout])
save([outputfolder '/' geomout '/sigma0.mat'], 'sigma0');
save([outputfolder '/' geomout '/T2.mat'], 'T2');
save([outputfolder '/' geomout '/b0fromhighb.mat'], 'b0fromhighb');
save([outputfolder '/' geomout '/META.mat'], 'META');
