% Script to run VERDICT with adaptive fitting

%% Folders...

% Image data folder
imagedatafolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Imaging Data\Patient Volunteers");

% Schemes folder
schemesfolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Code\VERDICT-Processing\Schemes");

% Models folder
modelsfolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Code\DWI-Noise-Project\Adaptive Model Fitting\MLP\models");

% Noise parameter etsimates folder
paramestfolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Code\DWI-Noise-Project\Adaptive Model Fitting\Outputs\Parameter Estimates");

% Output folder
outputfolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Code\DWI-Noise-Project\Adaptive Model Fitting\Outputs\VERDICT Outputs");


%% Fitting details

% Patient number
patnum = '20240912_Patient21';

% Model type
modeltype = 'Original VERDICT';

% Scheme name/ geometry
schemename = 'Original';

% Fitting technique (AMICO or MLP)
fittingtechnique = 'MLP';

% Fitting type (normal or adaptive)
fittingtype = 'adaptive';

% noise type
noisetype = 'Ratio';

% (for normal fitting)
sigma0train = 0.04;
T2train = 10000;


%% Preprocessing

% dicom folder
dfolder = [imagedatafolder '/' patnum '/' schemename];

% load scheme
load([schemesfolder '/' schemename]);

[Y, bvinY] = verdict_preprocess(dfolder, scheme);


%% Apply fitting

switch fittingtype

    case 'normal'

        % model folder
        thismodelfolder = [char(modelsfolder) '\' char(modeltype) '\' char(schemename) '\' noisetype '\T2_' num2str(T2train) '/sigma_' num2str(sigma0train) ];

        % Fitting
        [fIC, fEES, fVASC, R, rmse] = verdict_fit( ...
            Y, ...
            scheme, ...
            modeltype=modeltype, ...
            fittingtechnique=fittingtechnique, ...
            modelfolder=thismodelfolder ...
            );



    case 'adaptive'

        thismodelfolder = [char(modelsfolder) '\' char(modeltype) '\' char(schemename) '\' noisetype];

        % Load T2 and sigma0 estimate maps
        load([paramestfolder '/' patnum '/' schemename '/sigma0.mat'])
        load([paramestfolder '/' patnum '/' schemename '/T2.mat'])


        % NEED TO SORT THIS OUT!!!
        T2 = flip(T2, 3);
        sigma0=flip(sigma0, 3);

        % Load possible T2 and sigma0 values
        T2vals = load([modelsfolder '/' modeltype '/' schemename '/' noisetype '/T2s']).T2s;
        sigma0vals = load([modelsfolder '/' modeltype '/' schemename '/' noisetype '/sigma0s']).sigma0s;


        [fIC, fEES, fVASC, R, rmse] = adaptive_fit( ...
            Y,...
            schemename,...
            modeltype = modeltype,...
            modelfolder = modelsfolder,...
            schemesfolder = schemesfolder,...
            noisetype = noisetype,...
            T2 = T2,...
            sigma0 = sigma0,...
            T2vals = T2vals,...
            sigma0vals = sigma0vals...
            );
        

end

