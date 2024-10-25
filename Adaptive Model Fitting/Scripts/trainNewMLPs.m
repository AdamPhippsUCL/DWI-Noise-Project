% Script to train MLP models

% Language to configure MLP in (KEEP AS MATLAB)
configure = 'MATLAB';

% Noise 
noisetype = 'Ratio';
sigma0s = [0.01 0.02 0.04];
T2s = [40 50 60 80];

% Specify protocols 
modeltypes = {'No VASC VERDICT'} ;
schemenames = {'Short Scheme v3'}; 

pythonfolder = '';
TrainingDataFolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\DWI-Noise-Project\Adaptive Model Fitting\MLP\training data");
ModelFolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\DWI-Noise-Project\Adaptive Model Fitting\MLP\models");;

for sigma0 = sigma0s
    for T2 = T2s

        for indx = 1:length(modeltypes)
            
            modeltype = modeltypes{indx};
            schemename = schemenames{indx};
        
            datafolder = [TrainingDataFolder '/' modeltype '/' schemename '/' noisetype '/T2_' num2str(T2) '/sigma_' num2str(sigma0)];
            modelfolder = [ModelFolder '/' modeltype '/' schemename '/' noisetype '/T2_' num2str(T2) '/sigma_' num2str(sigma0)];
        
            trainMLP( ...
                datafolder,...
                modelfolder...
                );
            
        end

    end
end