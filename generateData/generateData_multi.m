function P = generateData_multi(P)
%GENERATEDATA This function takes the inputs in structure P and computes
%all the data required to perform the ptychographic reconstructions.
%==========================================================================

%% Load Parameters
N = P.N;
NP = P.NP;
ds_sam = P.ds_sam;
ds_det = P.ds_det;
NW = P.NW;
cenw = ceil(NW/2);

% Format precision and GPU:
if (P.useGPU == 1) && (P.useSinglePrec == 1)
    form = @(x) single(gpuArray(x));
elseif (P.useGPU == 1) && (P.useSinglePrec == 0)
    form = @(x) double(gpuArray(x));
elseif (P.useGPU == 0) && (P.useSinglePrec == 1)
    form = @(x) single(x);
elseif (P.useGPU == 0) && (P.useSinglePrec == 0)
    form = @(x) double(x);
end

%% Make Spectrum
[P,Sw,alpha] = makeSpectrum(P);

%% Make Probes
[P,prb0] = makeProbe(P);

% Use the same probe for each wavelength, but interpolated to the proper
% size:
y = ds_sam*(-(N/2) : (N/2-1));
x = y';

prb = zeros([size(prb0) NW]);
for nw = 1:NW
    xinterp = x*alpha(nw);
    yinterp = y*alpha(nw);

    prb_interp = interp2(x,y,prb0,xinterp,yinterp);
    prb_interp(isnan(prb_interp)) = 0;
    prb(:,:,nw) = sqrt(Sw(nw))*prb_interp*alpha(nw);
end

%% Make Scan Grid
[xpos,ypos] = makeScanGrid(P);
ds = ds_sam*alpha;
xposn = round(xpos./ds.'); P.xposn = xposn;
yposn = round(ypos./ds.'); P.yposn = yposn;

P.xpos = xpos;
P.ypos = ypos;

%% Make Objects
% Find number of pixels in object reconstruction based on probe pixel
% dimension and positions:
Nobj = N + (max(xposn,[],'all')-min(xposn,[],'all'));
Nobj = round(Nobj*1.1);

% Round to even number and save:
if rem(Nobj,2) == 1; Nobj = Nobj + 1; end
P.Nobj = Nobj;

% Build object:
[obj0,P] = makeObject(P);
obj = zeros([size(obj0) NW]);

% Use the same probe for each wavelength, but interpolated to the proper
% size:
y = ds_sam*(-(Nobj/2) : (Nobj/2-1));
x = y';

for nw = 1:NW
    xinterp = x*alpha(nw);
    yinterp = y*alpha(nw);

    obj_interp = interp2(x,y,obj0,xinterp,yinterp);
    obj_interp(isnan(obj_interp)) = 1;
    obj(:,:,nw) = obj_interp;
end

%% Calculate Diffraction Patterns
% Format the variables:
obj = form(obj);
prb = form(prb);

% Compute max dose and filter transmission:
cen = (P.Nobj/2) + (-N/2+1 : N/2);

% Define properly normalized Fourier transform:
FT2 = @(x) 1/N * fftshift(fftshift(fft2(x),1),2);

% Normalize probe intensity using flux:
for nw = 1:NW
    intensity = sum(abs(prb(:,:,nw)).^2,'all');
    prb(:,:,nw) = prb(:,:,nw) .* sqrt(P.prbFlux*Sw(nw)/intensity);
end
intensity = sum(abs(prb0).^2,'all');
prb0 = prb0 .* sqrt(P.prbFlux/intensity);

% Initialize diffraction data for each scan position:
[ESWA_ideal,ESWA_broad] = deal(form(zeros(N,N,NP)));

% Define central region of object:
cen = (Nobj/2) + (-N/2+1 : N/2);

% Make spatial frequency grid:
[FX,FY] = meshgrid( -N/2:N/2-1,-N/2:N/2-1 );
FX = repmat(FX/N,[1 1 NP]);
FY = repmat(FY/N,[1 1 NP]);

% Set up subpixel shifts in x and y:
xShift = xposn - round(xposn); xShift = reshape(xShift,[1 NW NP]);
yShift = yposn - round(yposn); yShift = reshape(yShift,[1 NW NP]);

% Pre-calculate subpixel shift kernels (NxNxNW):
fsh = @(x) fftshift(fftshift(x,1),2);
subPixKernel = zeros(N,N,NP,NW);
for nw = 1:NW

%     subPixKernel(:,:,:,nw) = fsh(exp(-2i*pi*(FX.*xShift(1,nw,:)/Nobj ...
%         + FY.*yShift(1,nw,:)/Nobj)));
    subPixKernel(:,:,:,nw) = exp(-2i*pi*(FX.*xShift(1,nw,:)/Nobj ...
        + FY.*yShift(1,nw,:)/Nobj));
end
conjSubPixKernel = conj(subPixKernel);

% For each position...
for np = 1:NP

    % Define sub-region of object:
    ynp = cen - round(yposn(:,np));
    xnp = cen + round(xposn(:,np));
    subObj = obj(ynp(cenw,:),xnp(cenw,:),cenw);

    % Apply subpixel shift:
    subObj = subPixelShift(subObj,conjSubPixKernel(:,:,np,cenw));

    % Propagate exit surface wave to get diffraction pattern:
    esw = subObj.*prb0;
    ESWA_ideal(:,:,np) = abs(FT2(esw));

    % If desired, make images of ESWA and eswa at each position:
    if P.makeESWAIm == 1
        plotESWA(ESWA_ideal(:,:,np),abs(esw),np);
    end

    % For each wavelength in the bandwidth...
    for nw = 1:NW

        % Create sub-object:
        subObj = obj(ynp(nw,:),xnp(nw,:),nw);

        % Apply subpixel shift:
        subObj = subPixelShift(subObj,conjSubPixKernel(:,:,np,nw));

        % Calculate broadband pattern as incoherent addition:
        I_nw = abs(FT2(subObj.*prb(:,:,nw))).^2;
        ESWA_broad(:,:,np) = ESWA_broad(:,:,np) + I_nw;
    end
end

% Ideal monochromatic and broadband diffraction patterns:
P.ESWA_ideal = gather(ESWA_ideal);
P.ESWA_broad = gather(sqrt(ESWA_broad));

% Apply noise to ideal broadband patterns and save:
P.ESWA_measured = addNoise(P.ESWA_broad,P);

% Save variables:
P.exa_prb = gather(prb);
P.exa_obj = gather(obj);

% Save P:
name = ['data_' P.name '.mat'];
save(name,'P','-v7.3');

end