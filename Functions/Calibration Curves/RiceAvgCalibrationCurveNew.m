% MATLAB Script to generate calibration curves for averaging Rice
% distributions at different SNR and NSA values


% Save folder
folder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\General VERDICT Code\General-VERDICT-Code\Signal Simulation\Calibration Curves";

% Number of samples
Nsample = 100000;

NSAs = [48, 60];

% Define SNR
SNRs = linspace(0.2,5,49);

save([char(folder) '/SNRs.mat'], 'SNRs')


% Loop over NSAs
for NSA = NSAs

    FitVs = zeros(size(SNRs));
    FitSigmas = zeros(size(SNRs));

    % Loop over SNRs
    for SNRIndx = 1:length(SNRs)

        SNR = SNRs(SNRIndx);

        disp([num2str(NSA) ', ' num2str(SNR) ]);
        v = 1;
        sigma = 1/SNR;


        % == Make Rice distribution

        ricedist = makedist('Rician','s',v,'sigma',sigma);
        
        maxsignal = max([4*(1/SNR), 2]);
        signals = linspace(0,maxsignal,1001);
        
        pdf = ricedist.pdf(signals);


        % == Take multi-average samples of Rice distribution
        
        samples = zeros(Nsample, 1);
        
        for indx = 1:Nsample
        
            % Average some signal measurements
            sample = 0;
            for AVindx = 1:NSA
                sample = sample + sampleDistribution(pdf, signals);
            end

            sample = sample/NSA;
        
            samples(indx, 1) = sample;
        
        end

        % f = figure;
        % H=histogram(samples);
        % hold on
        % plot(signals, pdf*H.BinWidth*Nsample)
        % 
        % 

        %% Fit Rician distribution to multi-average samples 
        x = (samples) + eps;
        pd = fitdist(x,'Rician');
        
        fitv = pd.s;
        FitVs(SNRIndx) = fitv;
    
        fitsigma = pd.sigma;
        FitSigmas(SNRIndx) = fitsigma;

        % 
        % fitpdf = pd.pdf(signals);
        % plot(signals, fitpdf*H.BinWidth*Nsample)
        % close(f)
        % 
        % if fitsigma/sigma > 1/sqrt(NSA)
        %     fitsigma = sigma/sqrt(NSA);
        % end



    end


    % Plot results
    f = figure;
    plot(SNRs, FitVs);
    hold on
    plot(SNRs, (FitSigmas.*SNRs)*sqrt(NSA))
    close(f);



    % Make dictionarys

    try
        mkdir([char(folder) '/NSA ' num2str(NSA)])
    catch
        disp('')
    end

    MeanDict = dictionary(SNRs, FitVs);
    save([char(folder) '/NSA ' num2str(NSA) '/MeanCalibration.mat'], 'MeanDict' )

    SigmaDict = dictionary(SNRs, FitSigmas.*SNRs);
    save([char(folder) '/NSA ' num2str(NSA) '/SigmaCalibration.mat'], 'SigmaDict' )










end
