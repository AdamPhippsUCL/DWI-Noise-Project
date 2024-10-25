

FittingResults = load("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\DWI-Noise-Project\Phantom Experiments\Outputs\24-Oct-2024 11-54-03\FittingResults.mat").FittingResults;
FittingResultsDN = load("C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\DWI-Noise-Project\Phantom Experiments\Outputs\24-Oct-2024 12-04-04\FittingResults.mat").FittingResults;

% figure;
% hold on
% scatter(ones(1,8), [FittingResults(1:8).sigma0fit]./[FittingResultsDN(1:8).sigma0fit], '*');
% boxplot( [FittingResults(1:8).sigma0fit]./[FittingResultsDN(1:8).sigma0fit], Positions = [1]);
% scatter(ones(1,8), [FittingResults(9:16).sigma0fit]./[FittingResultsDN(1:8).sigma0fit], '*');
% boxplot( [FittingResults(9:16).sigma0fit]./[FittingResultsDN(9:16).sigma0fit], Positions=[2]);
% boxplot( [FittingResults(17:24).sigma0fit]./[FittingResultsDN(17:24).sigma0fit], Positions=[3]);
% boxplot( [FittingResults(25:32).sigma0fit]./[FittingResultsDN(25:32).sigma0fit], Positions=[4]);
% boxplot([FittingResults(33:40).sigma0fit]./[FittingResultsDN(33:40).sigma0fit], Positions=[5]);
% boxplot( [FittingResults(41:48).sigma0fit]./[FittingResultsDN(41:48).sigma0fit], Positions=[6]);
% ylim([0 5])


% Number of ROIs
ROINames = unique([FittingResults(:).ROIName]);
NROI = length(ROINames);

figure;
hold on

for ROIindx = 1:NROI
    
    ROIName = ROINames(ROIindx);
    whereROI = (string({FittingResults.ROIName}) == ROIName);

    Nimg = sum(whereROI);

    scatter(ROIindx*ones(1,Nimg), [FittingResults(1 + Nimg*(ROIindx-1):Nimg*ROIindx).sigma0fit]./[FittingResultsDN(1 + Nimg*(ROIindx-1):Nimg*ROIindx).sigma0fit], '*');
    boxplot( [FittingResults(1 + Nimg*(ROIindx-1):Nimg*ROIindx).sigma0fit]./[FittingResultsDN(1 + Nimg*(ROIindx-1):Nimg*ROIindx).sigma0fit], Positions = [ROIindx]);

end

ylim([0 5]);
