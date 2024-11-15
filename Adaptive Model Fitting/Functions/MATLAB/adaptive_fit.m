% Function to implement adaptive fitting
% (Calls verdict_fit.m with different masks)

function [fIC, fEES, fVASC, R, rmse] = adaptive_fit(Y, schemename, opts)

arguments
    Y % Normalised data array [nx, ny, nz, nscheme]
    schemename % scheme structure name

    % == OPTIONS
    
    opts.modeltype 
    opts.modelfolder
    opts.schemesfolder
    opts.noisetype
    opts.T2 % T2 map [nx, ny, nz]
    opts.sigma0 % sigma0 map [nx, ny, nz]

    opts.T2vals % possible network T2 values
    opts.sigma0vals % possible network sigma0 values

end


% Load scheme
load([opts.schemesfolder '/' schemename]);

% Round T2 and sigma0 maps to possible values
T2map = roundtowardvec(opts.T2, opts.T2vals);
sigma0map = roundtowardvec(opts.sigma0, opts.sigma0vals);


% Initialise arrays
fIC = zeros(size(T2map));
fEES = zeros(size(T2map));
fVASC = zeros(size(T2map));
R = zeros(size(T2map));
rmse = zeros(size(T2map));

for T2 = opts.T2vals
    for sigma0 = opts.sigma0vals

        thismask = and(T2map == T2, sigma0map == sigma0);
        sum(thismask(:));


        thismodelfolder = [opts.modelfolder '/' opts.modeltype '/' schemename '/' opts.noisetype '/T2_' num2str(T2) '/sigma_' num2str(sigma0)];

        % Call verdict_fit function with mask
        [thisfIC, thisfEES, thisfVASC, thisR, thisrmse] = verdict_fit( ...
            Y, ...
            scheme,...
            modeltype=opts.modeltype,...
            fittingtechnique = 'MLP',...
            mask = thismask,...
            modelfolder = thismodelfolder);


        fIC(thismask) = thisfIC(thismask);
        fEES(thismask) = thisfEES(thismask);
        fVASC(thismask) = thisfVASC(thismask);
        R(thismask) = thisR(thismask);
        rmse(thismask) = thisrmse(thismask);
    end
end

end