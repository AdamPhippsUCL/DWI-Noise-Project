% Create MATLAB structure containing image details

ImageDetails = struct();

% Protocol names
ProtocolNames = {
    "b1500_Ex1",...
    "b2000_Ex1",...
    "b2100_Ex1",...
    "b2400_Ex1",...
    "b2500_Ex1",...
    "b2750_Ex1",...
    "b3000_Ex1",...
    "b3500_Ex1",...
    "b4000_Ex1",...
    ...
    "b2500_NSA2_Ex2",...
    "b2500_NSA4_Ex2",...
    "b2500_NSA6_Ex2",...
    "b2500_NSA8_Ex2",...
    ...
    "b2500_Avg1_Ex3",...
    "b2500_Avg2_Ex3",...
    "b2500_Avg4_Ex3",...
    "b2500_Avg6_Ex3",...
    "b2500_Avg8_Ex3",...    
    };


% b values
bvalues = {
    1500,...
    2000,...
    2100,...
    2400,...
    2500,...
    2750,...
    3000,...
    3500,...
    4000,...
    ...
    2500,...
    2500,...
    2500,...
    2500,...
    ...
    2500,...
    2500,...
    2500,...
    2500,...
    2500,...
};


% TEs
TEs = {
    450,...
    200,...
    100,...
    234,...
    300,...
    450,...
    400,...
    150,...
    250,...
    ...
    300,...
    300,...
    300,...
    300,...
    ...
    300,...
    300,...
    300,...
    300,...
    300,...
    };

% TRs
TRs = {
    4000,...
    4000,...
    4000,...
    4000,...
    4000,...
    4000,...
    4000,...
    4000,...
    4000,...
    ...
    4000,...
    4000,...
    4000,...
    4000,...
    ...
    4000,...
    4000,...
    4000,...
    4000,...
    4000,...
    };


% NSA
NSAs = {
    1,...
    1,...
    1,...
    1,...
    1,...
    1,...
    1,...
    1,...
    1,...
    ...
    2,...
    4,...
    6,...
    8,...
    ...
    2,...
    2,...
    2,...
    2,...
    2,...
};


% Rav
Ravs = {
    3,...
    3,...
    3,...
    3,...
    3,...
    3,...
    3,...
    3,...
    3,...
    ...
    3,...
    3,...
    3,...
    3,...
    ...
    3*1,...
    3*2,...
    3*4,...
    3*6,...
    3*8,...
    };


% Append all to structure
ImageDetails.ProtocolName = ProtocolNames;
ImageDetails.bvalue = bvalues;
ImageDetails.TE = TEs;
ImageDetails.TR = TRs;
ImageDetails.NSA = NSAs;
ImageDetails.Rav = Ravs;


% Save structure as mat file
folder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\Noise Statistics Project\Code\DWI-Noise-Project\Phantom Experiments\Imaging Data";
save([char(folder) '/ImageDetails.mat' ], 'ImageDetails')
