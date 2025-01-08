%COMPAREALGOS Compares ePIE to PaCMAN for a single probe overlap,
%bandwidth, and flux.
%==========================================================================

close all; clear; clc
P.name = 'single';

%% Section 1: Initialize Data
% Switch parameters (1=yes,0=no):
P.useGPU = 1;               % use GPU
P.useSinglePrec = 0;        % use single precision
P.useShotNoise = 1;         % apply shot (Poisson) noise
P.useDetectorNoise = 1;     % apply detector (Gaussian) noise
P.useBackground = 1;        % apply parasitic scattering background
P.makeESWAIm = 0;           % display generated surface waves

% Detector noise parameters:
P.QE = 0.8;                 % quantum efficiency (0 to 1)
P.sigmaReadout = 2;         % standard deviation of noise (in counts)
P.bitDepth = 16;            % bit depth of camera

% Illumination parameters:
P.prbFlux = 1e8;                % number of photons incident on object
P.lambda = 4.13e-9;             % 300 eV, center wavelength [m]
P.rel_bw = 0.2;                 % relative bandwidth (delta_lambda/lambda)
P.lam_bw = P.lambda*P.rel_bw;   % actual FWHM intensity bandwidth [m]
P.NW = 33;                     % number of wavelength samples in bandwidth

% Geometrical parameters:
P.overSam = 4;              % oversampling (larger value, smaller probe)
P.N = 256;                  % num pixels in diff pat (NxN)
P.z = 5e-2;                 % distance from sample to detector [m]
P.ds_det = (1.04e-5)*P.overSam;             % eff. detector pix. size [m]
P.ds_sam = (P.lambda*P.z)/(P.ds_det*P.N);   % object pixel size [m]

% Square scan grid parameters:
P.probeOverlap = 0.7;       % Overlap between the probe radii (~0.4-0.9)
P.randOffsetFrac = 0.2;     % Scan pos offset as fraction of step (~0.2):
P.NP = 13^2;                % Total number of scan positions
                            % (enabling makeESWAIm helps to pick NP)

% Transmission bounds for object:
P.tMin = 0.5;               % min target amplitude transmission (max is 1)
P.pMax = 1;                 % max phase shift in radians (min is 0)

%% Section 2: Generate Data
P = generateData(P);

%% Section 3: ePIE Reconstructions
rP.NIT = 200;       % number of iterations
rP.NPRB = 10;       % iteration at which probe starts updating
rP.plotObjPrb = 1;  % plot object/probe at each iteration (1=yes,0=no)

% Object/probe gain:
rP.beta_obj = 0.5;
rP.beta_prb = 0.5;

% Load data:
load(['data_' P.name '.mat'],'P')

% Run ePIE reconstruction with measured data:
ESWA = P.ESWA_measured;

% Set parasitic scattering region to zero:
load('mask.mat','mask');
mask(mask == 1) = 2; mask(mask == 0) = 1; mask(mask == 2) = 0;
ESWA = ESWA.*mask;

% Initialize object and probe (use [] for default initializations):
obj = [];
prb = [];

% Reconstruct:
rP = ePIE(rP,P,obj,prb,ESWA);

% Save variables:
save(['reconstruction_' P.name '_epie.mat'],'rP')
fprintf('ePIE reconstruction is finished.\n')

%% Section 4: PaCMAN Reconstructions
% Iteration parameters:
rP.NIT = 200;       % number of iterations
rP.NPRB = 10;       % iteration at which probe starts updating
rP.NN = 20;         % iteration at which noise is corrected
rP.NS = 1;          % total number of spatial modes to use (default=1)
rP.plotObjPrb = 1;  % plot object/probe at each iteration (1=yes,0=no)

% Regularizing parameter:
rP.r = 1.5;

% Select monochromatization algorithm ("CGLS" or "BiCGSTAB") and maximum
% number of iterations to run:
rP.monoAlgo = "BiCGSTAB";
rP.kmax = 50;

% Load data:
load(['data_' P.name '.mat'],'P')

% Set parasitic scattering region to zero:
ESWA = P.ESWA_measured;
load('mask.mat','mask')
mask(mask == 1) = 2; mask(mask == 0) = 1; mask(mask == 2) = 0;
ESWA = ESWA.*mask;

% ===================== No Monochromatization =============================
% Initialize object and probe (use [] for default initializations):
obj = [];
prb = [];

% Reconstruct:
rP = PaCMAN(rP,P,obj,prb,ESWA);

% Save variables:
save(['reconstruction_' P.name '_PaCMAN.mat'],'rP')
fprintf('PaCMAN (no monochromatization) reconstruction is finished.\n')

% ==================== With Monochromatization ============================
% Initialize object and probe (use [] for default initializations):
obj = [];
prb = [];

% Monochromatize:
ESWA_mono = monochromatize(ESWA,P,rP);

% Reconstruct:
rP = PaCMAN(rP,P,obj,prb,ESWA_mono);

% Save variables:
save(['reconstruction_' P.name '_PaCMAN_' ...
    convertStringsToChars(rP.monoAlgo) '.mat'],'rP')
fprintf('PaCMAN (BiCGSTAB) reconstruction is finished.\n')

%% Section 5: Analyze Reconstructions
% Identify central region of object:
cen = P.cen;
cen = cen(5:(end-5));
NC = length(cen);

% Select camera region for NRMSE calculation:
camx = 30:85;
camy = 90:145;

% Initialize normalized root-mean-square error (NRMSE)/object matrices:
objects = zeros(NC,NC,3);
errors = zeros(1,3);

% Exact object (for error calculations):
exa_obj = P.exa_obj(cen,cen);

% ePIE:
load(['reconstruction_' P.name '_epie.mat'],'rP')
obj = rP.rec_obj(cen,cen);
objects(:,:,1) = obj;
errors(1) = NRMSE(obj(camx,camy),exa_obj(camx,camy));

% PaCMAN without monochromatization:
load(['reconstruction_' P.name '_PaCMAN.mat'],'rP')
obj = rP.rec_obj(cen,cen);
objects(:,:,2) = obj;
errors(2) = NRMSE(obj(camx,camy),exa_obj(camx,camy));

% PaCMAN with monochromatization (BiCGSTAB):
load(['reconstruction_' P.name '_PaCMAN_BiCGSTAB.mat'],'rP')
obj = rP.rec_obj(cen,cen);
objects(:,:,3) = obj;
errors(3) = NRMSE(obj(camx,camy),exa_obj(camx,camy));

% Save the analyses:
save([P.name '_analysis.mat'],'errors','objects','exa_obj');

%% Section 6: Plot
% Make a 1x3 plot that is zoomed in on the camera, with the normalized
% root-mean-square error (NRMSE) plotted above.

% Object tranmission limits for colorbar:
clb = min(abs(P.exa_obj),[],'all');
cub = max(abs(P.exa_obj),[],'all');

fig = figure('windowstate','maximized');
subplot(1,3,1)
imagesc(abs(objects(camx,camy,1)))
axis off; axis square; colorbar; caxis([clb cub])
title(['ePIE: NRMSE = ' num2str(errors(1))])
set(gca,'fontsize',16)

subplot(1,3,2)
imagesc(abs(objects(camx,camy,2)))
axis off; axis square; colorbar; caxis([clb cub])
title(['PaCMAN (no mono): NRMSE = ' num2str(errors(2))])
set(gca,'fontsize',16)

subplot(1,3,3)
imagesc(abs(objects(camx,camy,3)))
axis off; axis square; colorbar; caxis([clb cub])
title(['PaCMAN (BiCGSTAB): NRMSE = ' num2str(errors(3))])
set(gca,'fontsize',16)
