function rP = ePIE(rP,P,obj,prb,ESWA)
%EPIE Runs the extended ptychographic iterative engine (ePIE), returning
%the reconstructed object and probe.
%==========================================================================

%% Initialize
% Reconstruction parameters:
NIT = rP.NIT;
NPRB = rP.NPRB;
beta_prb = rP.beta_prb;
beta_obj = rP.beta_obj;

% Data parameters:
prbFlux = P.prbFlux;
NP = P.NP;
xposn = P.xposn;
yposn = P.yposn;
N = P.N;

% Random number generator seed:
rng(1)

%% Object/Probe Guesses
% Make reconstruction the same size as exact object:
[NY,NX] = size(P.exa_obj);

% If no object guess is given, initialize object as ones:
if isempty(obj)
    obj = ones(NY,NX);
end

% If no probe guess is given, initialize probe as uniform intensity within
% the support:
if isempty(prb)
    prb = P.exa_prb;
    prb(prb>0) = max(prb,[],'all');
end

%% Fourier Transforms
% Shift ESWA to save time:
ESWA = fftshift(fftshift(ESWA,1),2);

% Scale FTs so that Parseval's theorem is obeyed:
FT2 = @(x) fft2(x)/N;
IFT2 = @(x) ifft2(x)*N;       

%% Format Precision and Processing
if (P.useGPU == 1) && (P.useSinglePrec == 1)
    form = @(x) single(gpuArray(x));
elseif (P.useGPU == 1) && (P.useSinglePrec == 0)
    form = @(x) double(gpuArray(x));
elseif (P.useGPU == 0) && (P.useSinglePrec == 1)
    form = @(x) single(x);
elseif (P.useGPU == 0) && (P.useSinglePrec == 0)
    form = @(x) double(x);
end

% Apply desired format:
obj = form(obj);
prb = form(prb);
xposn = form(xposn);
yposn = form(yposn);
ESWA = form(ESWA);

%% Reconstruct
% Define 1xN vectors that represent the pixel values of the indexed object
% at the center of the full object grid:
xcen = (NX/2) + (-N/2+1 : N/2);
ycen = (NY/2) + (-N/2+1 : N/2);

% Initialize figure if desired:
if rP.plotAll == 1; fig = figure('windowstate','maximized'); end

% For each iteration...
for n = 1:NIT

    % Randomize scan position order:
    order = randperm(NP);

    % For each scan position...
    for np = 1:NP

        % Select random position:
        pos = order(np);
        
        % Use position to create subobject with same dim as prb:
        xp = xcen + xposn(pos);
        yp = ycen - yposn(pos);
        obj_p = obj(yp,xp);
        
        % Calculate exit surface wave and propagate to detector:
        esw = obj_p.*prb;
        ESW = FT2(esw);

        % Enforce modulus constraint in detector plane:
        esw_p = IFT2(ESWA(:,:,pos) .* exp(1i*angle(ESW)));

        % Update object:
        esw_diff = esw_p - esw;
        prb_maxSquared = max(abs(prb(:))).^2;
        obj(yp,xp) = obj_p + beta_obj*conj(prb).*esw_diff./ prb_maxSquared;

        % Update probe:
        if n >= NPRB
            obj_maxSquared = max(abs(obj_p(:))).^2;
            prb = prb + beta_prb*conj(obj_p).*esw_diff/obj_maxSquared;

            % Apply constraints to probe (optional):
            prb(P.exa_prb == 0) = 0;                        % support
            prb = prb*sqrt(prbFlux/sum(abs(prb).^2,'all')); % total flux
        end
    end

    % Update plot of obj, prb, and errors:
    if rP.plotAll == 1
        plotIterations(fig,NIT,n,obj,prb,P);
    end
end

% Save reconstructed object and probe:
rP.rec_obj = gather(obj);
rP.rec_prb = gather(prb);

end