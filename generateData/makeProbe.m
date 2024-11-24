function [P,prb] = makeProbe(P)
%MAKEPROBE Makes a Gaussian probe with a hard aperture.
%==========================================================================

%% Load Parameters
N = P.N;
ds_sam = P.ds_sam;
overSam = P.overSam;

%% Make Probes:
% Make position grid and find probe radius [m]:
y = ds_sam*(-(N/2) : (N/2-1));
x = y';
r = sqrt(x.^2 + y.^2);
prbR = ds_sam * N/(overSam*2);

%%
% Create Gaussian probe with hard aperture:
prb = exp(-(r/prbR).^2) .* exp(-1i*(r.^2 / ds_sam^2) /1e3);
prb( abs(prb) < exp(-1) ) = 0;

% Update P:
P.prbR = prbR;
end

