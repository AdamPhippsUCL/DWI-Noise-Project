% Testing parameter estimation code


% % 20240912 Patient 26
% Echo1fname = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Imaging Data\Patient Volunteers\20241018_Patient26\Series 501 [MR - b1te50P4]\1.3.6.1.4.1.5962.99.1.2688242394.826411277.1729265128154.116.0.dcm";
% Echo2fname = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Imaging Data\Patient Volunteers\20241018_Patient26\Series 601 [MR - b1te125P4]\1.3.6.1.4.1.5962.99.1.2688242394.826411277.1729265128154.255.0.dcm";
% b1800fname = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Imaging Data\Patient Volunteers\20241018_Patient26\Series 301 [MR - b1800shV4]\1.3.6.1.4.1.5962.99.1.2688242394.826411277.1729265128154.8.0.dcm";
% outputfolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Code\DWI-Noise-Project\Adaptive Model Fitting\Outputs\Parameter Estimates\20241018_Patient26\shV4");
% b1800 = dicomread(b1800fname);
% b0fromhighb = squeeze(b1800(:,:,1,1:2:end));

% % 20240912 Patient 22
% Echo1fname = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Imaging Data\Patient Volunteers\20240912_Patient22\Series 601 [MR - b1te50P3]\1.3.6.1.4.1.5962.99.1.3883698028.1870441860.1726165616492.1325.0.dcm";
% Echo2fname = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Imaging Data\Patient Volunteers\20240912_Patient22\Series 701 [MR - b1te125P3]\1.3.6.1.4.1.5962.99.1.3883698028.1870441860.1726165616492.1464.0.dcm";
% b1800fname = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Imaging Data\Patient Volunteers\20240912_Patient22\Series 401 [MR - b1800shV3]\1.3.6.1.4.1.5962.99.1.3883698028.1870441860.1726165616492.1218.0.dcm";
% outputfolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Code\DWI-Noise-Project\Adaptive Model Fitting\Outputs\Parameter Estimates\20240912_Patient22\shV3");
% b1800 = dicomread(b1800fname);
% b0fromhighb = squeeze(b1800(:,:,1,1:2:end));


% 20240911 Patient 21
Echo1fname = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Imaging Data\Patient Volunteers\20240912_Patient21\shV3\Series 601 [MR - b1te50P3]\1.3.6.1.4.1.5962.99.1.3883698028.1870441860.1726165616492.115.0.dcm";
Echo2fname = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Imaging Data\Patient Volunteers\20240912_Patient21\shV3\Series 701 [MR - b1te125P3]\1.3.6.1.4.1.5962.99.1.3883698028.1870441860.1726165616492.254.0.dcm";
b1800fname = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Imaging Data\Patient Volunteers\20240912_Patient21\shV3\Series 401 [MR - b1800shV3]\1.3.6.1.4.1.5962.99.1.3883698028.1870441860.1726165616492.8.0.dcm";
outputfolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Code\DWI-Noise-Project\Adaptive Model Fitting\Outputs\Parameter Estimates\20240912_Patient21\shV3");
b1800 = dicomread(b1800fname);
b0fromhighb = squeeze(b1800(:,:,1,1:2:end));

% 
% % 20240911 Patient 20
% Echo1fname = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Imaging Data\Patient Volunteers\20240911_Patient20\Series 601 [MR - b1te50P3]\1.3.6.1.4.1.5962.99.1.3763965273.1080835530.1726045883737.149.0.dcm";
% Echo2fname = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Imaging Data\Patient Volunteers\20240911_Patient20\Series 701 [MR - b1te125P3]\1.3.6.1.4.1.5962.99.1.3763965273.1080835530.1726045883737.288.0.dcm";
% b1800fname = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Imaging Data\Patient Volunteers\20240911_Patient20\Series 401 [MR - b1800shV3]\1.3.6.1.4.1.5962.99.1.3763965273.1080835530.1726045883737.427.0.dcm";
% outputfolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Code\DWI-Noise-Project\Adaptive Model Fitting\Outputs\Parameter Estimates\20240911_Patient20\shV3");
% b1800 = dicomread(b1800fname);
% b0fromhighb = squeeze(b1800(:,:,1,1:2:end));

% % 20240628 Patient 10
% Echo1fname = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Imaging Data\Patient Volunteers\20240628_Patient10\Series 1201 [MR - b1 NSA1LOW6te50P]\1.3.6.1.4.1.5962.99.1.1967159687.1975438338.1719954110855.349.0.dcm";
% Echo2fname = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Imaging Data\Patient Volunteers\20240628_Patient10\Series 1401 [MR - b1 NSA1LOW6te125P]\1.3.6.1.4.1.5962.99.1.1967159687.1975438338.1719954110855.99.0.dcm";
% b1800fname = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Imaging Data\Patient Volunteers\20240628_Patient10\Series 1001 [MR - b1800shVP]\1.3.6.1.4.1.5962.99.1.1967159687.1975438338.1719954110855.1335.0.dcm";
% outputfolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Code\DWI-Noise-Project\Adaptive Model Fitting\Outputs\Parameter Estimates\20240628_Patient10\shV3");
% b1800 = dicomread(b1800fname);
% b0fromhighb = squeeze(b1800(:,:,1,1:5:end));


mkdir(outputfolder);


patchsize = 3;

[sigma0, T2, resnorm] = EstimateNoiseParameters( ...
    Echo1fname, ...
    Echo2fname, ...
    b0fromhighb=b0fromhighb, ...
    register=true, ...
    patchsize = patchsize);



% Save parameter estimates
save([outputfolder '/sigma0.mat'], 'sigma0');
save([outputfolder '/T2.mat'], 'T2');
save([outputfolder '/b0fromhighb.mat'], 'b0fromhighb');


% Meta
META = struct();
META.patchsize = patchsize;
META.Echo1fname = Echo1fname;
META.Echo2fname = Echo2fname;
save([outputfolder '/META.mat'], "META");

