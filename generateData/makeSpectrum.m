function [P,Sw,alpha] = makeSpectrum(P)
%MAKESPECTRUM Makes spectrum and generates C for monochromatization.
%==========================================================================

%% Load Parameters
lam_bw = P.lam_bw;
lam0 = P.lambda;
NW = P.NW;
N = P.N;
rel_bw = P.rel_bw;

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

%% Generate C:
% If C for this bandwidth, number of wavelengths, and number of pixels
% doesn't exist, generate it:
name = ['C_bw' num2str(rel_bw*100) '_pix' num2str(N) ...
    '_lam' num2str(NW) '.mat'];

% if ~exist(name,'file')
%     % Start a parallel pool with local CPU cores (accelerates generation of
%     % C by 2-3 times), generate C, and terminate parallel pool:
%     parpool
%     C = generateC(P);
%     pp = gcp('nocreate');
%     delete(pp);
% 
%     % Save to unique name:
%     save(name,'C')
% else
%     load(name,'C');
% end
% P.C = C;

end

