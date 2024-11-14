% Script to create new training data for MLP model training

% Base folder
basef = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Code\VERDICT-Processing");

% Number of training data samples
Nvoxel = 20000;

% Noise 
noisetype = 'Rice';
sigma0s = [ 0.02, 0.04, 0.06];
T2s = [ 10000];

% Training data folder
savedata = true;
TrainingDataFolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Code\DWI-Noise-Project\Adaptive Model Fitting\MLP\training data");

% % Save T2 and sigma0 values
% save([TrainingDataFolder '/sigma0s.mat'], "sigma0s")
% save([TrainingDataFolder '/T2s.mat'], "T2s")

%% Protocols

% If multiple, define as a cell array of char arrays
modeltypes = {'Original VERDICT'};
schemenames = {'Original'};

% Saving sheme as mat file
schemesfolder = [basef '/Schemes'];
savescheme = true;


%% Cell sizes

% Cell Radius distribution R~Normal( muR, sigmaR)
muRmin = 6;
muRmax = 9;
sigmaRmin = 0.99;
sigmaRmax = 1.01;


%% Create training data

for sigma0 = sigma0s
    for T2 = T2s
        for indx = 1:length(modeltypes)
        
            modeltype = modeltypes{indx};
            schemename = schemenames{indx};
        
            disp([modeltype ' ' schemename])
        
            outputfolder = [TrainingDataFolder '/' modeltype '/' schemename '/' noisetype '/T2_' num2str(T2) '/sigma_' num2str(sigma0)];
            
            createVERDICTdata( ...
                modeltype,...
                schemename,...
                Nvoxel = Nvoxel,...
                noisetype = noisetype,...
                sigma0 = sigma0,...
                T2 = T2,...
                randmuRs=[muRmin, muRmax],...
                randsigmaRs=[sigmaRmin, sigmaRmax],...
                schemesfolder=schemesfolder,...
                savedata=savedata,...
                outputfolder=outputfolder...
                );
        
        end
    end
end
