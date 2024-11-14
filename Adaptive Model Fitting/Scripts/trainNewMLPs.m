% Script to train MLP models

% Language to configure MLP in (KEEP AS MATLAB)
configure = 'MATLAB';

% Noise 
noisetype = 'Rice';
sigma0s = [0.02 0.04 0.06];
T2s = [10000];

% Specify protocols 
modeltypes = {'Original VERDICT'} ;
schemenames = {'Original'}; 

pythonfolder = '';
TrainingDataFolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Code\DWI-Noise-Project\Adaptive Model Fitting\MLP\training data");
ModelFolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Code\DWI-Noise-Project\Adaptive Model Fitting\MLP\models");

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


for indx = 1:length(modeltypes)
    
    modeltype = modeltypes{indx};
    schemename = schemenames{indx};
    save([ModelFolder '/' modeltype '/' schemename '/' noisetype '/T2s.mat' ], "T2s");
    save([ModelFolder '/' modeltype '/' schemename '/' noisetype '/sigma0s.mat' ], "sigma0s");

end