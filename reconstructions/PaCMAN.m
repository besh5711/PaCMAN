function rP = PaCMAN(rP,P,obj,prb,ESWA)
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

% Data parameters:
xposn = P.xposn;
yposn = P.yposn;
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
    obj = ones(NY,NX,NS);
end

% If no probe guess is given, initialize probe as uniform intensity within
% the support:
if isempty(prb)
    prb = P.exa_prb;
    prb(prb>0) = max(prb,[],'all');
end

% Initialize spatial modes, and throw error if we don't initialize all of
% them:
prb0 = prb;
prb = zeros(N,N,NS);
prb(:,:,1) = prb0;
% prb(:,:,2) = (1/10)*max(abs(prb0),[],'all')*exp(1i*randn(size(prb0)));
% prb(:,:,3) = (1/100)*max(abs(prb0),[],'all')*exp(1i*randn(size(prb0)));
if min(sum(abs(prb(:,:,NS)),[1 2])) == 0
    error(['The number of probes must match the number of ' ...
        'reconstructed modes! Define more modes in PaCMAN.'])
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

% 1xN vectors containining probe position:
xp = xcen + xposn(1:NP).';
yp = ycen - yposn(1:NP).';

% Initialize PaCMAN variables:
[Lambda,z,PSI,subObj] = deal(form(zeros(N,N,NP,NS)));
Lambda_hat = form(zeros(N,N,NP));
mubar = form(zeros(N,N));

% Initialize subobject matrix including each position and mode:
for np = 1:NP
    for ns = 1:NS
        subObj(:,:,np,ns) = obj(yp(np,:),xp(np,:),ns);
    end
end

% Initialize z (these are fftshifted):
for ns = 1:NS
    z(:,:,:,ns) = FT2(subObj(:,:,:,ns) .* prb(:,:,ns));
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
    
        % Apply support constraint to primary mode only (optional):
        temp = prb(:,:,1);
        temp(P.exa_prb == 0) = 0;
        prb(:,:,1) = temp;
    
        % Apply flux total flux constraint (optional):
        prb = prb*sqrt(prbFlux/sum(abs(prb).^2,'all'));    
    end
    
    %========================= Object Update ==============================  
    obj_den_sum = abs(prb).^2;
    
    % For each position...
    for np = 1:NP
        
        % Update numerator:
        obj_num(yp(np,:),xp(np,:),:) = obj_num(yp(np,:),xp(np,:),:) ...
            + conj(prb).*squeeze(zbar(:,:,np,:));
    
        % Update denominator:
        obj_den(yp(np,:),xp(np,:),:) = obj_den(yp(np,:),xp(np,:),:) ...
            + obj_den_sum;
    end
    
    % Update object:
    obj = (obj + obj_num) ./ (1 + obj_den + eps);
    
    %======================= Update Field =================================  
    % Update sub-object matrix:
    for np = 1:NP
        for ns = 1:NS
            subObj(:,:,np,ns) = obj(yp(np,:),xp(np,:),ns);
        end
    end
    
    % Calculate propagated updated exit waves (fftshifted) for this probe:
    for ns = 1:NS
        PSI(:,:,:,ns) = FT2(subObj(:,:,:,ns).*prb(:,:,ns));
    end
    
    %====================== Update Regularizers ===========================  
    % Temporary variables:
    modX = sqrt(sum((abs(PSI - Lambda)).^2,4) ...
        + abs(mubar-Lambda_hat).^2 + eps);
    rho = (r*modX +sqrt((r*modX).^2 + 4*(1+r).*ESWA.^2))./(2*(1+r));
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
        temp_avg = (1/NP) * sum(ESWA.^2 - sum(abs(PSI).^2,4),3);
        temp_avg(temp_avg < 0) = 0;
        mubar = sqrt(temp_avg);
    end

    % Update plot of obj, prb, and errors:
    if rP.plotObjPrb == 1
        plotIterations(fig,NIT,n,obj,prb,P);
    end
end

% Save reconstructed object and probe:
rP.rec_obj = gather(obj);
rP.rec_prb = gather(prb);

end