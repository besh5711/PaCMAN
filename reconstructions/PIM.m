function rP = PIM(rP,P,obj,prb,ESWA)
%PIM Runs Ptychographic Information Multiplexing (PIM), returning the
%reconstructed object and probe.
%==========================================================================

%% Initialize
% Reconstruction parameters:
NIT = rP.NIT;
NPRB = rP.NPRB;
beta_prb = rP.beta_prb;
beta_obj = rP.beta_obj;
NW = rP.NW;
cenw = ceil(NW/2);

% Data parameters:
prbFlux = P.prbFlux;
NP = P.NP;
N = P.N;
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
    obj = ones(NY,NX,NW);
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

    % Initialize probe at each wavelength as a copy of this guess:
    prb = repmat(prb,[1 1 NW]);
    
    % Grid for center wavelength:
    yd = ds_sam*(-(N/2) : (N/2-1));
    xd = yd.';
    
    % Scale each wavelength-dependent probe:
    prb0 = prb(:,:,cenw);
    for nw = 1:NW
        xinterp = xd*alpha(nw);
        yinterp = yd*alpha(nw);
    
        prb_interp = interp2(xd,yd,prb0,xinterp,yinterp);
        prb_interp(isnan(prb_interp)) = 0;
        prb(:,:,nw) = sqrt(Sw(nw)/Sw(cenw))*prb_interp*alpha(nw);
    end

    % Mask to enforce support (optional):
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
[esw,ESW,subObj] = deal(form(zeros(N,N,NW)));

%% Reconstruct
% Define 1xN vectors that represent the pixel values of the indexed object
% at the center of the full object grid:
xcen = (NX/2) + (-N/2+1 : N/2);
ycen = (NY/2) + (-N/2+1 : N/2);

% Initialize figure if desired:
if rP.plotObjPrb == 1; fig = figure('windowstate','maximized'); end

% For each iteration...
for n = 1:NIT

    % Randomize scan position order:
    order_p = randperm(NP);

    % Randomize wavelength order:
    order_w = randperm(NW);

    % For each scan position...
    for np = 1:NP
        pos = order_p(np);

        % For each wavelength... (Loop 1)
        for nw = 1:NW
            wl = order_w(nw);

            % Use position to create subobject with same dim as prb:
            xp = xcen + round(xposn(wl,pos));
            yp = ycen - round(yposn(wl,pos));
            subObj(:,:,wl) = obj(yp,xp,wl);
            
            % Calculate exit surface wave and propagate to detector:
            esw(:,:,wl) = subObj(:,:,wl) .* prb(:,:,wl);
            ESW(:,:,wl) = FT2(esw(:,:,wl));
        end

        % Take the incoherent sum of the diffraction patterns:
        Itot = sum(abs(ESW).^2,3);

        % For each wavelength... (Loop 2)
        for nw = 1:NW
            wl = order_w(nw);

            % Normalize this wavelength:
            ESWPrime = ESW(:,:,wl).*ESWA(:,:,pos)./sqrt(Itot);
            
            % Backpropagate and take difference in esw:
            esw_new = IFT2(ESWPrime);
            esw_diff = esw_new - esw(:,:,wl);

            % Update subobject:
            prb_max2 = max(abs(prb(:,:,wl)),[],'all').^2;
            subObj_new = subObj(:,:,wl) + beta_obj ...
                * conj(prb(:,:,wl)) .* esw_diff / prb_max2;

            % Update full object:
            xp = xcen + round(xposn(wl,pos));
            yp = ycen - round(yposn(wl,pos));
            obj(yp,xp,wl) = subObj_new;

            % Update probe:
            if n >= NPRB
                obj_max2 = max(abs(subObj(:,:,wl)),[],'all').^2;
                prb(:,:,wl) = prb(:,:,wl) + beta_prb*conj(subObj(:,:,wl)) ...
                    .* esw_diff / obj_max2;

                % Apply support constraint (optional):
                prb(:,:,wl) = prb(:,:,wl).*mask(:,:,wl);

                % Apply spectral weight constraint (optional):
                % prb(:,:,wl) = prb(:,:,wl)*sqrt(prbFlux*Sw(wl) ...
                %     /sum(abs(prb(:,:,wl)).^2,'all'));
            end
        end

        % Apply total probe flux constraint (optional):
        if n >= NPRB
            prb = prb * sqrt(prbFlux/sum(abs(prb).^2,'all'));
        end
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