function [P,Sw,alpha] = makeSpectrum(P)
%MAKESPECTRUM Makes Gaussian spectrum.
%==========================================================================

%% Load Parameters
lam_bw = P.lam_bw;
lam0 = P.lambda;
NW = P.NW;

%% Create Spectrum
% Create Gaussian spectrum and normalize so sum is 1:
lams = linspace(-lam_bw,lam_bw,NW) + lam0;
Sw = exp(-4*log(2)*((lams-lam0)/lam_bw).^2 );
Sw = Sw/sum(Sw);

% Calculate scaling factor:
alpha = lams/lam0;

% Save:
P.Sw = Sw;
P.alpha = alpha;
P.lams = lams;

end

