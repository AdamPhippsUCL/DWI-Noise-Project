% Define function
function out = func(x)

    yprime = evalin('base', 'yprime');
    BetaAlphaRatios = evalin('base', 'BetaAlphaRatios');
    SNRs = evalin('base', 'SNRs');

    [~,I] = min( abs(x-SNRs));
    y = BetaAlphaRatios(I)*SNRs(I);

    out = y-yprime;

end