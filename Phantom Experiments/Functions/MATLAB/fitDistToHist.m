function [coeffs, resnorm] = fitDistToHist(counts, bincenters, params, opts)

arguments

    counts % Array of counts
    bincenters % Array of bin centers

    % Distribution parameters 
    params.sigma0 = NaN
    params.T2 = NaN
    params.fd = NaN
    params.TE = 50 % Specify at [TE1, TE2] for DTIRatio
    params.N0 = 1
    params.Nb = 1

    opts.disttype = 'Ratio'
    opts.beta0guess = [0.05, 250, 0.8]

    opts.boundstype = 'fixed';
    opts.lb = [0.001, 50, 0.001]
    opts.ub = [0.5, 1500, 1]
    opts.boundsfrac = 0.25

end


%% Steps

% 1. Normalise histogram counts
% 2. Define unknown parameter vector for fitting
% 3. Apply fitting

%% 1. Normalise histogram counts
counts = counts/sum(counts);


%% 2. Define unknown parameters for fitting

Nparam = 3;
if isfield(params, 'fd')
    if isnan(params.fd)
        params = rmfield(params, 'fd');
    else
        Nparam = Nparam-1; % Just fitting for sigma0 and T2; x = [sigma0, T2]
    end
end

if isfield(params, 'T2')
    if isnan(params.T2)
        params = rmfield(params, 'T2');
    else
        Nparam = Nparam-1; % Just fitting for sigma0; x = sigma0
    end
end



%% Fit

% Use nlm function: fitnlm

% Assign variable to base (for access in functions)
assignin('base', 'current_params', params)


% Define bounds
switch opts.boundstype

    case 'fixed'
        % Bounds set by user
        disp('');

    case 'dependent'
        opts.lb = (1-opts.boundsfrac)*opts.beta0guess;
        opts.ub = (1-opts.boundsfrac)*opts.beta0guess;


end


options = optimoptions('lsqcurvefit', 'OptimalityTolerance', 1e-12);


switch opts.disttype

    case 'Ratio'
        if Nparam == 1
            [coeffs, resnorm] = lsqcurvefit(@RatioDist1, opts.beta0guess(1), bincenters, counts, opts.lb(1), opts.ub(1), options);   
        elseif Nparam == 2
            [coeffs, resnorm] = lsqcurvefit(@RatioDist2, opts.beta0guess(1:2), bincenters, counts, opts.lb(1:2), opts.ub(1:2), options );
        elseif Nparam == 3
            [coeffs, resnorm] = lsqcurvefit(@RatioDist3, opts.beta0guess(1:3), bincenters, counts, opts.lb(1:3), opts.ub(1:3), options);
        end

    case 'Rice'

        if Nparam == 1
            [coeffs, resnorm] = lsqcurvefit(@RiceDist1, opts.beta0guess(1), bincenters, counts, opts.lb(1), opts.ub(1:3) );  
        elseif Nparam == 2
            [coeffs, resnorm] = lsqcurvefit(@RiceDist2, opts.beta0guess(1:2), bincenters, counts, opts.lb(1:2), opts.ub(1:2) );   
        elseif Nparam == 3
            [coeffs, resnorm] = lsqcurvefit(@RiceDist3, opts.beta0guess(1:3), bincenters, counts, opts.lb(1:3), opts.ub(1:3) ); 
        end


end

if Nparam == 1
    coeffs = [coeffs, params.T2, params.fd];
elseif Nparam == 2
    coeffs = [coeffs, params.fd];
elseif Nparam == 3
    coeffs = [coeffs];
end


end


% Function to fit Ratio distribution

function binfreqs = RatioDist1(b, x)

% Generate ratio distribution given x
% Outputs bin frequencies at bin centers predicted
% for ratio distribution with sigma0

arguments
    b % sigma0 value
    x % bin centers
end

% Get parameters
params = evalin('base', 'current_params');
T2 = params.T2;
TE = params.TE;
fd = params.fd;
N0 = params.N0;
Nb = params.Nb;
sigma0 = b;

% b0 signal
b0signal = exp(-TE/T2);

% b signal
bsignal = fd*b0signal;

% Bin spacings
binmin = min(x);
binmax = max(x);
binspacing = x(2)-x(1);


% Generate pdf over bin centers
pdfvals = RatioDistRician( ...
    b0signal, ...
    bsignal, ...
    sigma0, ...
    N0=N0, ...
    Nb=Nb, ...
    zmin = binmin, ...
    zmax = binmax, ...
    dz = binspacing,...
    ymin = binmin,...
    ymax = binmax,...
    dy = binspacing);

% Evaluate bin frequencies
binfreqs = pdfvals*binspacing;

binfreqs = reshape(binfreqs, size(x));

end

function binfreqs = RatioDist2(b, x)

% Generate ratio distribution given x
% Outputs bin frequencies at bin centers predicted
% for ratio distribution with sigma0

arguments
    b % [sigma0, T2] value
    x % bin centers
end

% Get parameters
params = evalin('base', 'current_params');
TE = params.TE;
fd = params.fd;
N0 = params.N0;
Nb = params.Nb;
sigma0 = b(1);
T2 = b(2);

% b0 signal
b0signal = exp(-TE/T2);

% b signal
bsignal = fd*b0signal;

% Bin spacings
binmin = min(x);
binmax = max(x);
binspacing = x(2)-x(1);

% Generate pdf over bin centers
pdfvals = RatioDistRician( ...
    b0signal, ...
    bsignal, ...
    sigma0, ...
    N0=N0, ...
    Nb=Nb, ...
    zmin = binmin, ...
    zmax = binmax, ...
    dz = binspacing,...
    ymin = binmin,...
    ymax = binmax,...
    dy = binspacing);

% Evaluate bin frequencies
binfreqs = pdfvals*binspacing;

binfreqs = reshape(binfreqs, size(x));


end

function binfreqs = RatioDist3(b, x)

% Generate ratio distribution given x
% Outputs bin frequencies at bin centers predicted
% for ratio distribution with sigma0

arguments
    b % [sigma0, T2] value
    x % bin centers
end

% Get parameters
params = evalin('base', 'current_params');
TE = params.TE;
N0 = params.N0;
Nb = params.Nb;
sigma0 = b(1);
T2 = b(2);
fd = b(3);

% b0 signal
b0signal = exp(-TE/T2);

% b signal
bsignal = fd*b0signal;

% Bin center
bincentres = x;

% Bin spacings
binmin = min(x);
binmax = max(x);
binspacing = x(2)-x(1);

% Generate pdf over bin centers
pdfvals = RatioDistRician( ...
    b0signal, ...
    bsignal, ...
    sigma0, ...
    N0=N0, ...
    Nb=Nb, ...
    zs = bincentres,...
    ys = linspace(0,2,250),...
    ymin = binmin,...
    ymax = binmax,...
    dy = binspacing);

% Evaluate bin frequencies
binfreqs = pdfvals*binspacing;

binfreqs = reshape(binfreqs, size(x));


end


% Function to fit Rice distribution

function binfreqs = RiceDist1(b, x)

% Generate Rice distribution given x
% Outputs bin frequencies at bin centers predicted
% for ratio distribution with sigma0

arguments
    b % sigma0 value
    x % bin centers
end

% Get parameters
params = evalin('base', 'current_params');
fd = params.fd;
sigma0 = b;

% Bin spacings
binmin = min(x);
binmax = max(x);
binspacing = x(2)-x(1);

% Generate pdf over bin centers
pdfvals = RiceDist(1, fd, sigma0, zmin = binmin, zmax = binmax, dz = binspacing);

% Evaluate bin frequencies
binfreqs = pdfvals*binspacing;

binfreqs = reshape(binfreqs, size(x));

end

function binfreqs = RiceDist2(b, x)

% Generate ratio distribution given x
% Outputs bin frequencies at bin centers predicted
% for ratio distribution with sigma0

arguments
    b % [sigma0, T2] values
    x % bin centers
end

% Get parameters
params = evalin('base', 'current_params');
fd = params.fd;
sigma0 = b(1);

% Bin spacings
binmin = min(x);
binmax = max(x);
binspacing = x(2)-x(1);

% Generate pdf over bin centers
pdfvals = RiceDist(1, fd, sigma0, zmin = binmin, zmax = binmax, dz = binspacing);

% Evaluate bin frequencies
binfreqs = pdfvals*binspacing;

binfreqs = reshape(binfreqs, size(x));

end

function binfreqs = RiceDist3(b, x)

% Generate ratio distribution given x
% Outputs bin frequencies at bin centers predicted
% for ratio distribution with sigma0

arguments
    b % [sigma0, T2, fd] values
    x % bin centers
end

% Get parameters
params = evalin('base', 'current_params');
sigma0 = b(1);
fd = b(3);

% Bin center
bincentres = x;

% Bin spacings
binmin = min(x);
binmax = max(x);
binspacing = x(2)-x(1);

% Generate pdf over bin centers
pdfvals = RiceDist(1, fd, sigma0, zs = bincentres);

% Evaluate bin frequencies
binfreqs = pdfvals*binspacing;

binfreqs = reshape(binfreqs, size(x));

end
