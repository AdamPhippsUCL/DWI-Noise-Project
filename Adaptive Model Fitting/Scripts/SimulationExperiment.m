% MATLAB code to set up and run simulation experiment

% Two grids with varying T2 and sigma0 values across the grid

%% Define grid 

% Number of voxels per grid section
ngrid = 30;

% fIC values
fICmin = 0.05;
fICmax = 0.95;
NfIC = 19;
fICs = linspace(fICmin, fICmax, NfIC);

% T2 values
T2s = [50, 60, 80];
NT2 = length(T2s);

% sigma0
sigma0s = 0.02*ones(size(T2s)); 
Nsigma0 = length(sigma0s);

% fIC grid
fICgrid = zeros(NT2*ngrid, NfIC*ngrid);
for fICindx = 1:NfIC
    for gridindx = 1:ngrid
        fICgrid( :, (fICindx-1)*ngrid + gridindx) = fICs(fICindx);
    end
end


% T2 grid
T2grid = zeros(NT2*ngrid, NfIC*ngrid);
for T2indx = 1:NT2
    for gridindx = 1:ngrid
        T2grid( (T2indx-1)*ngrid + gridindx, :) = T2s(T2indx);
    end
end


% sigma0 grid
% sigma0grid = sigma0*ones(size(T2grid));
sigma0grid = zeros(Nsigma0*ngrid, NfIC*ngrid);
for sigma0indx = 1:Nsigma0
    for gridindx = 1:ngrid
        sigma0grid( (sigma0indx-1)*ngrid + gridindx, :) = sigma0s(sigma0indx);
    end
end

% Rs
muRs = [7 8];
sigmaRs = [1 1];

Rs = linspace(0.1, 15.1, 17);



%% SIMULATE SIGNALS OVER GRID

schemesfolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Code\VERDICT-Processing\Schemes");

modeltype = 'Original VERDICT';
schemename = 'Original';

load([schemesfolder '/' schemename '.mat'])
nscheme = length(scheme);

% Signal grid
signalgrid = zeros(NT2*ngrid, NfIC*ngrid, nscheme);
normsignalgrid = zeros(NT2*ngrid, NfIC*ngrid, nscheme);


% Simulate signals over grid
for jndx = 1:size(T2grid,1)
    for indx = 1:size(fICgrid, 2)

        disp([num2str(jndx) ', ' num2str(indx) ])
        
        sigma0 = sigma0grid(jndx, indx);
        T2 = T2grid(jndx, indx);

        % Rs
        muR = muRs(1) + (muRs(2)-muRs(1))*rand();
        sigmaR = sigmaRs(1) + (sigmaRs(2)-sigmaRs(1))*rand();
        fRs = normpdf(Rs, muR, sigmaR);
        fRs = (1/sum(fRs))*fRs;

        % fIC
        fIC = fICgrid(jndx, indx);
        
        % fVASC
        switch modeltype
            case 'Original VERDICT'
                fVASC = min(0.1, 1-fIC)*rand()
            case 'No VASC VERDICT'
                fVASC = 0;
        end
        % fEES
        fEES = 1-fIC-fVASC;


        for schemeIndx = 1:length(scheme)

            scan_params = [scheme(schemeIndx).delta, scheme(schemeIndx).DELTA, scheme(schemeIndx).bval];

            if scan_params(3)==0;
                normsignalgrid(jndx,indx, schemeIndx) = 1;
                continue;
            end

            TE = scheme(schemeIndx).TE;
            NSA = scheme(schemeIndx).NSA;
            Rav = scheme(schemeIndx).Rav;

            % Tissue parameter vector
            switch modeltype
                case 'Original VERDICT'
                    tps = [fIC*fRs, fEES ,fVASC];
                case 'No VASC VERDICT'
                    tps = [fIC*fRs, fEES];
            end
            
            % Diffusion weighting fraction
            fd = simulateSignal( ...
                tps, ...
                scan_params, ...
                modeltype,...
                Rs = Rs...
                );


            % ==== b=0 signal distribution

            % Noiseless signal
            b0signal = exp(-TE/T2);

            % Define Rice distribution
            b0dist = makedist('Rician', s=b0signal, sigma = sigma0 );
            
            signals = linspace(0,2,201);
            b0pdf = b0dist.pdf(signals);



            % ==== b>0 signal distribution

            % Noiseless signal
            bsignal = b0signal*fd;


            % Define Rice distribution
            bdist = makedist('Rician', s=bsignal, sigma = sigma0 );
            
            signals = linspace(0,2,201);
            bpdf = bdist.pdf(signals);


            % == b=0 signal

            % Take NSA samples
            sample = 0;
            for s = 1:NSA
                sample = sample + sampleDistribution(b0pdf, signals);
            end
            b0sample = sample/NSA;  

            signalgrid(jndx, indx, schemeIndx-1) = b0sample;

            % == b>0 signal

            % Take Nb samples
            Nb = NSA*Rav;
            sample = 0;
            for s = 1:Nb
                sample = sample + sampleDistribution(bpdf, signals);
            end
            bsample = sample/Nb;    

            signalgrid(jndx, indx, schemeIndx) = bsample;

            % == Normalised signal
            normsignalgrid(jndx, indx, schemeIndx) = bsample/b0sample;



        end


    end
end



%% MP-PCA DENOISING

usedenoise = false;



if usedenoise

    normsignalgridDN = ones(size(signalgrid));

    windowsize = 3;
    denoisedsignalgrid = zeros(size(signalgrid));
    
    for ischeme = 1:round(nscheme/2)
    
        [denoised,S2,P]= denoise(signalgrid(:,:,2*ischeme-1:2*ischeme),[windowsize windowsize]);
        signalgridDN(:,:,2*ischeme-1:2*ischeme) = denoised;

        normsignalgridDN(:,:,2*ischeme) = signalgridDN(:,:,2*ischeme)./signalgridDN(:,:,2*ischeme-1);

    
    end


    Y = normsignalgridDN;
    
else

    Y = normsignalgrid;

end




%% APPLY FITTING

% == Settings
fittingtype = 'adaptive';  % 'normal' using Rice network; 'adaptive' using Ratio networks

% for normal fitting
fittingtechnique = 'MLP';
noisetype = 'Rice';
T2train = 10000;
sigma0train = 0.04;

modelfolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Code\DWI-Noise-Project\Adaptive Model Fitting\MLP\Simulation Experiment\models");

switch fittingtype


    case 'normal'

        thismodelfolder = [modelfolder '/' modeltype '/' schemename '/' noisetype '/T2_' num2str(T2train) '/sigma_' num2str(sigma0train) ];
        
        [fIC, fEES, fVASC, R, rmse] = verdict_fit( ...
            Y, ...
            scheme, ...
            modeltype = modeltype,...
            fittingtechnique=fittingtechnique,...    
            modelfolder=thismodelfolder...
            );


    case 'adaptive'

        % Load possible T2 and sigma0 values
        T2vals = load([modelfolder '/' modeltype '/' schemename '/Ratio/T2s']).T2s;
        sigma0vals = load([modelfolder '/' modeltype '/' schemename '/Ratio/sigma0s']).sigma0s;

        [fIC, fEES, fVASC, R, rmse] = adaptive_fit( ...
            Y,...
            schemename,...
            modeltype = modeltype,...
            modelfolder = modelfolder,...
            schemesfolder = schemesfolder,...
            noisetype = 'Ratio',...
            T2 = T2grid,...
            sigma0 = sigma0grid,...
            T2vals = T2vals,...
            sigma0vals = sigma0vals...
            );


end


% % reshape
% RicefIC = reshape(RicefIC, size(fICgrid));
% save([char(outputfolder) '/RicefIC.mat'], 'RicefIC');
% 
% 
% 
% % == Adaptive Ratio Networks
% 
% T2train = T2grid;
% sigma0train = sigma0grid;
% 
% [AdaptivefIC, AdaptivefEES, AdaptivefVASC, AdaptiveR, Adaptivermse] = verdict_MLP_fit( ...
%     schemename, ...
%     modeltype, ...
%     Y, ...
%     noisetype='Adaptive', ...
%     T2train = T2train, ...
%     sigma0train=sigma0train, ...
%     schemesfolder=schemesfolder,...
%     modelsfolder=modelsfolder,...
%     outputfolder = outputfolder...
%     );
% 
% % reshape
% AdaptivefIC = reshape(AdaptivefIC, size(fICgrid));
% save([char(outputfolder) '/AdaptivefIC.mat'], 'AdaptivefIC');
% 


%% Display fitting results
% 
% f=figure;
% imshow(fICgrid, [0, 1] );
% c = colorbar;
% ylabel(c, 'f_{IC}')
% h = gca;
% h.Visible = 'On';
% yticks([round(ngrid/2) + ngrid*(0:NT2-1) ])
% yticklabels(T2s)
% ylabel('T2 (ms)')
% xticks([round(ngrid/2) + ngrid*(0:NfIC-1) ])
% xticklabels([])
% % xlabel('f_{IC} (ms)')
% title('Ground Truth')
% ax = gca();
% ax.FontSize=12;
% 
% 
% figure;
% imshow(fIC, [0, 1]);
% c = colorbar;
% ylabel(c, 'f_{IC}')
% h = gca;
% h.Visible = 'On';
% yticks([round(ngrid/2) + ngrid*(0:NT2-1) ])
% yticklabels(T2s)
% ylabel('T2 (ms)')
% xticks([round(ngrid/2) + ngrid*(0:NfIC-1) ])
% xticklabels([])
% % xlabel('f_{IC} (ms)')
% switch fittingtype
%     case 'normal'
%         switch fittingtechnique           
%             case 'AMICO'            
%                 title('AMICO')
%             case 'MLP'
%                 title('MLP (Rice)')
%         end
%     case 'adaptive'
%         title('adaptive')
% end
% 
% 
% ax = gca();
% ax.FontSize=12;


% saveas(f, [outputfolder '/fIC grids.fig'])
% saveas(f, [outputfolder '/fIC grids.png'])

%% Bias and std results

outputfolder = char("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD\Projects\DWI Noise Project\Code\DWI-Noise-Project\Adaptive Model Fitting\Outputs\Simulation");

% Difference maps
Diffs = fIC-fICgrid;

% Biases
Biases = zeros(NT2, NfIC);

% Variances
Variances = zeros(NT2, NfIC);

for fICindx = 1:NfIC
    for T2indx = 1:NT2

        vals = Diffs((T2indx-1)*ngrid+1:T2indx*ngrid, (fICindx-1)*ngrid+1:fICindx*ngrid);
        bias = mean(vals(:));
        vari = var(vals(:));
        Biases(T2indx, fICindx) = bias;
        Variances(T2indx, fICindx) = vari;

    end
end


% BIAS FIGURE

f1=figure;
for T2indx = 1:NT2
    T2 = T2s(T2indx);
    plot(fICs, Biases(T2indx,:), '-*', DisplayName = ['T2 = ' num2str(T2) ' ms'])
    hold on
end
ylim([-0.05 0.5])
xlim([min(fICs)-0.05, max(fICs)+0.05])
xlabel('f_{IC}')
ylabel('Bias')
switch fittingtype
    case 'normal'
        title(['Rice (\sigma_0 = ' num2str(sigma0s(1)) ')' ])
    case 'adaptive'
        title(['Ratio (\sigma_0 = ' num2str(sigma0s(1)) ')' ])
end
grid on
ax = gca();
ax.FontSize=12;
legend;
f1.Position = [400, 200, 350, 400];



% VARIANCE FIGURE

f2=figure;
for T2indx = 1:NT2
    T2 = T2s(T2indx);
    plot(fICs, sqrt(Variances(T2indx,:)), '-*', DisplayName = ['T2 = ' num2str(T2) ' ms'])
    hold on
end
ylim([-0.015 0.15])
xlim([min(fICs)-0.05, max(fICs)+0.05])
xlabel('f_{IC}')
ylabel('Standard deviation')
switch fittingtype
    case 'normal'
        title(['Rice (\sigma_0 = ' num2str(sigma0s(1)) ')' ])
    case 'adaptive'
        title(['Ratio (\sigma_0 = ' num2str(sigma0s(1)) ')' ])
end
grid on
ax = gca();
ax.FontSize=12;
% legend;
f2.Position = [400, 200, 350, 400];


%% Save figures

switch fittingtype
    case 'normal'
        mkdir([outputfolder '/' schemename '/Rice/sigma0_' num2str(sigma0s(1))])
        saveas(f1, [outputfolder '/' schemename '/Rice/sigma0_' num2str(sigma0s(1)) '/Bias.fig'])
    case 'adaptive'
        mkdir([outputfolder '/' schemename '/Ratio/sigma0_' num2str(sigma0s(1))])
        saveas(f1, [outputfolder '/' schemename '/Ratio/sigma0_' num2str(sigma0s(1)) '/Bias.fig'])
end

switch fittingtype
    case 'normal'
        saveas(f2, [outputfolder '/' schemename '/Rice/sigma0_' num2str(sigma0s(1)) '/Variance.fig'])
    case 'adaptive'
        saveas(f2, [outputfolder '/' schemename '/Ratio/sigma0_' num2str(sigma0s(1)) '/Variance.fig'])
end


% 
% %% Fitting residual error
% 
% % Resimulate signals using model parameter estimate and look at mean
% % squared error between 
% 
% resimsignalgrid = zeros(size(signalgrid));
% 
% for indx = 1:size(signalgrid,1)
%     for jindx = 1:size(signalgrid, 2)
% 
% 
%         % simulate signal with fitted values 
% 
% 
%     end
% end
% 
