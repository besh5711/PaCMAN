function rP = Ms_PaCMAN(rP,P,obj,prb,ESWA)
%PACMAN Runs the Partial Coherence and Monochromatization Algorithm with
%Noise (PaCMAN), returning the reconstructed object and probe.
%==========================================================================

%% Initialize
% Reconstruction parameters:
NIT = rP.NIT;
NPRB = rP.NPRB;
NN = rP.NN;
NS = rP.NS;
r = rP.r;
NW = rP.NW;
cenw = ceil(NW/2);
cenw_gen = ceil(P.NW/2);

% Data parameters:
prbFlux = P.prbFlux;
N = P.N;
NP = P.NP;
lam_bw = P.lam_bw;
lam0 = P.lambda;
xpos = P.xpos;
ypos = P.ypos;
ds_sam = P.ds_sam;

% Random number generator seed:
rng(1)

%% Object/Probe Guesses
% Make reconstruction the same size as exact object:
[NY,NX,~] = size(P.exa_obj);

% If no object guess is given, initialize object as ones:
if isempty(obj)
    obj = ones(NY,NX,NS,NW);
end

% % If no probe guess is given, initialize probe as uniform intensity within
% % the support:
% if isempty(prb)
%     prb = P.exa_prb;
%     prb(prb>0) = max(prb,[],'all');
% end

% Calculate which wavelengths to reconstruct and the corresponding probe
% positions:
lams = linspace(-lam_bw,lam_bw,NW) + lam0;
Sw = exp(-4*log(2)*((lams-lam0)/lam_bw).^2 );
Sw = Sw/sum(Sw);
alpha = lams/lam0;
ds = ds_sam*alpha;
xposn = xpos./ds.';
yposn = ypos./ds.';

% If no probe guess is given, initialize probe as uniform intensity within
% the support and scaled for each reconstructed wavelength:
if isempty(prb)

    prb = P.exa_prb;
    prb(prb>0) = max(prb,[],'all');

    % Initialize probe at each wavelength as a copy of this guess:
    prb = repmat(prb(:,:,ceil(P.NW/2)),[1 1 NW]);
    
    yd = P.ds_sam*(-(N/2) : (N/2-1));
    xd = yd.';
    
    % Calculate probes by scaling:
    prb0 = prb(:,:,cenw);
    for nw = 1:NW
        xinterp = xd*alpha(nw);
        yinterp = yd*alpha(nw);
    
        prb_interp = interp2(xd,yd,prb0,xinterp,yinterp);
        prb_interp(isnan(prb_interp)) = 0;
        prb(:,:,nw) = sqrt(Sw(nw)/Sw(cenw))*prb_interp*alpha(nw);
    end
    mask = (prb > 0);
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

% NWxNPxN vectors containining probe position:
xp = reshape(xcen,[1 1 N]) + round(xposn);
yp = reshape(ycen,[1 1 N]) - round(yposn);

% Initialize PaCMAN variables:
[Lambda,z,PSI,subObj] = deal(form(zeros(N,N,NP,NS,NW)));
Lambda_hat = form(zeros(N,N,NP));
mubar = form(zeros(N,N));

% Initialize subobject matrix including each position and mode:
for np = 1:NP
    for ns = 1:NS
        for nw = 1:NW
            subObj(:,:,np,ns,nw) = obj(yp(nw,np,:),xp(nw,np,:),ns,nw);
        end
    end
end

% Initialize z (these are fftshifted):
for ns = 1:NS
    for nw = 1:NW
        z(:,:,:,ns,nw) = FT2(subObj(:,:,:,ns,nw) .* prb(:,:,ns,nw));
    end
end

% Initialize figure if desired:
if rP.plotObjPrb == 1; fig = figure('windowstate','maximized'); end

% For each iteration...
for n = 1:NIT
    
    % Initialize (or reset) object numerators to 0:
    [obj_den,obj_num] = deal(form(zeros(size(obj))));
        
    % Calculate zbar (these are NOT fftshifted):
    zbar = IFT2(z + Lambda);

    %========================= Probe Update ===============================  
    % Update probe:
    if n >= NPRB
    
        % Calculate numerator of probe (NxN):
        prb_num = conj(subObj).*zbar;
        prb_num = squeeze(sum(prb_num,3));
    
        % Precalculate sum of subobjects over all positions:
        prb_den = squeeze(sum(abs(subObj).^2,3));
    
        % Update probe:
        prb = (prb + prb_num) ./ (1 + prb_den + eps);
    
        % Apply support constraint to primary spatial mode only (optional):
        prb(:,:,1,:) = prb(:,:,1,:).*mask;
    
        % Apply flux total flux constraint (optional):
        prb = prb*sqrt(prbFlux/sum(abs(prb).^2,'all'));    
    end
    
    %========================= Object Update ==============================  
    obj_den_sum = abs(prb).^2;
    
    % For each position...
    for np = 1:NP

        % For each wavelength...
        for nw = 1:NW
        
            % Update numerator:
            obj_num(yp(nw,np,:),xp(nw,np,:),:,nw) = ...
                obj_num(yp(nw,np,:),xp(nw,np,:),:,nw) ...
                + conj(prb(:,:,:,nw)).*squeeze(zbar(:,:,np,:,nw));
        
            % Update denominator:
            obj_den(yp(nw,np,:),xp(nw,np,:),:,nw) = ...
                obj_den(yp(nw,np,:),xp(nw,np,:),:,nw)+obj_den_sum(:,:,nw);
        end
    end
    
    % Update object:
    obj = (obj + obj_num) ./ (1 + obj_den + eps);
    
    %======================= Update Field =================================  
    % Update sub-object matrix:
    for np = 1:NP
        for ns = 1:NS
            for nw = 1:NW
                subObj(:,:,np,ns,nw) = obj(yp(nw,np,:),xp(nw,np,:),ns,nw);
            end
        end
    end
    
    % Calculate propagated updated exit waves (fftshifted) for this probe:
    for ns = 1:NS
        for nw = 1:NW
            PSI(:,:,:,ns,nw) = FT2(subObj(:,:,:,ns,nw).*prb(:,:,ns,nw));
        end
    end
    
    %====================== Update Regularizers ===========================  
    % Temporary variables:
    modX = sqrt(sum((abs(PSI - Lambda)).^2,[4 5]) ...
        + abs(mubar-Lambda_hat).^2 + eps);
    rho = (r*modX +sqrt((r*modX).^2 + 4*(Cp+r).*Cp.*ESWA.^2))./(2*(Cp+r));
    coeff = rho./modX;
    
    % Update z, mu, and the average of mu:
    z = coeff.*(PSI - Lambda);
    mu = coeff.*(mubar - Lambda_hat);
    mubar = (1/NP)*sum(mu,3);

    %================= Update Background/Multipliers ======================      
    % Update Lambdas:
    Lambda = Lambda + z - PSI;
    Lambda_hat = Lambda_hat + mu - mubar;
    
    % Enable structured noise correction at iteration NN:
    if n == NN
        temp_avg = (1/NP) * sum(ESWA.^2 - sum(abs(PSI).^2,[4 5]),3);
        temp_avg(temp_avg < 0) = 0;
        mubar = sqrt(temp_avg);
    end

    % Update plot of obj, prb, and errors:
    if rP.plotObjPrb == 1
        plotIterations_multi(fig,NIT,n,obj,prb,P,cenw);
    end
end

% Save reconstructed object and probe:
rP.rec_obj = gather(obj);
rP.rec_prb = gather(prb);

end