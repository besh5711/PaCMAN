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

% Data parameters:
prbFlux = P.prbFlux;
N = P.N;
NP = P.NP;

% Random number generator seed:
rng(1)

%% Object/Probe Guesses
% Make reconstruction the same size as exact object:
[NY,NX] = size(P.exa_obj);

% If no object guess is given, initialize object as ones:
if isempty(obj)
    obj = ones(NY,NX,NS,NW);
end

% If no probe guess is given, initialize probe as uniform intensity within
% the support:
if isempty(prb)
    prb = P.exa_prb;
    prb(prb>0) = max(prb,[],'all');
end

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

    % Initialize spatial modes, and throw error if we don't initialize all
    % of them:
    prb0 = prb;
    prb = zeros(N,N,NS);
    prb(:,:,1) = prb0;
    % prb(:,:,2) = (1/10)*max(abs(prb0),[],'all')*exp(1i*randn(N,N));
    % prb(:,:,3) = (1/100)*max(abs(prb0),[],'all')*exp(1i*randn(N,N));
    if min(sum(abs(prb(:,:,NS)),[1 2])) == 0
        error(['The number of probes must match the number of ' ...
            'reconstructed modes! Define more modes in PaCMAN.'])
    end

    % Initialize wavelength-dependent probes:
    prb = repmat(prb,[1 1 1 NW]);
    
    % Grid for center wavelength:
    yd = ds_sam*(-(N/2) : (N/2-1));
    xd = yd.';
    
    % Scale each wavelength- and mode-dependent probe:
    for ns = 1:NS
        prb0 = squeeze(prb(:,:,ns,:));
        for nw = 1:NW
            xinterp = xd*alpha(nw);
            yinterp = yd*alpha(nw);
        
            prb_interp = interp2(xd,yd,prb0,xinterp,yinterp);
            prb_interp(isnan(prb_interp)) = 0;
            prb(:,:,ns,nw) = sqrt(Sw(nw)/Sw(cenw))*prb_interp*alpha(nw);
        end
    end

    % Mask to enforce support on primary spatial mode (optional):
    mask = (prb(:,:,1,:) > 0);
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
if rP.plotAll == 1; fig = figure('windowstate','maximized'); end

% For each iteration...
for n = 1:NIT
    
    % Initialize (or reset) object numerators to 0:
    [obj_den,obj_num] = deal(form(zeros(size(obj))));
        
    % Calculate zbar (these are NOT fftshifted) [Step 1]:
    zbar = IFT2(z + Lambda);

    %========================= Probe Update ===============================  
    % Update probe:
    if n >= NPRB
    
        % Calculate numerator of probe (NxN):
        prb_num = conj(subObj).*zbar;
        prb_num = squeeze(sum(prb_num,3));
    
        % Precalculate sum of subobjects over all positions:
        prb_den = squeeze(sum(abs(subObj).^2,3));
    
        % Update probe [Step 2]:
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
    
    % Update object [Step 3]:
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
    % [Step 4]
    for ns = 1:NS
        for nw = 1:NW
            PSI(:,:,:,ns,nw) = FT2(subObj(:,:,:,ns,nw).*prb(:,:,ns,nw));
        end
    end
    
    %====================== Update Regularizers ===========================  
    % Temporary variables [Step 4]:
    modX = sqrt(sum((abs(PSI - Lambda)).^2,[4 5]) ...
        + abs(mubar-Lambda_hat).^2 + eps);
    rho = (r*modX +sqrt((r*modX).^2 + 4*(Cp+r).*Cp.*ESWA.^2))./(2*(Cp+r));
    coeff = rho./modX;
    
    % Update z, mu, and the average of mu [Step 5]:
    z = coeff.*(PSI - Lambda);
    mu = coeff.*(mubar - Lambda_hat);
    mubar = (1/NP)*sum(mu,3);

    %================= Update Background/Multipliers ======================      
    % Update Lambdas [Step 6]:
    Lambda = Lambda + z - PSI;
    Lambda_hat = Lambda_hat + mu - mubar;
    
    % Enable structured noise correction at iteration NN:
    if n == NN
        temp_avg = (1/NP) * sum(ESWA.^2 - sum(abs(PSI).^2,[4 5]),3);
        temp_avg(temp_avg < 0) = 0;
        mubar = sqrt(temp_avg);
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