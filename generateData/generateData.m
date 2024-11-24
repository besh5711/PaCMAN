function P = generateData(P)
%GENERATEDATA This function takes the inputs in structure P and computes
%all the data required to perform the ptychographic reconstructions.
%==========================================================================

%% Load Parameters
N = P.N;
NP = P.NP;
ds_sam = P.ds_sam;
ds_det = P.ds_det;
NW = P.NW;

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

%% Make Probe
[P,prb] = makeProbe(P);

%% Make Scan Grid
[xpos,ypos] = makeScanGrid(P);
xposn = round(xpos/ds_sam); P.xposn = xposn;
yposn = round(ypos/ds_sam); P.yposn = yposn;

%% Make Object
[obj,P] = makeObject(P);

%% Calculate Diffraction Patterns
% Format the variables:
obj = form(obj);
prb = form(prb);

% Compute max dose and filter transmission:
cen = (P.Nobj/2) + (-N/2+1 : N/2);

% Define properly normalized Fourier transform:
FT2 = @(x) 1/N * fftshift(fftshift(fft2(x),1),2);

% Normalize probe intensity using flux:
intensity = sum(abs(prb).^2,'all');
prb = prb .* sqrt(P.flux/intensity);

% Initialize diffraction data for each scan position:
[ESWA_ideal,ESWA_broad] = deal(form(zeros(N,N,NP)));

% Make detector grid for wavelength interpolation:
yd = ds_det*(-(N/2) : (N/2-1));
xd = yd;
[xd,yd] = meshgrid(xd,yd);

% For each position...
for np = 1:NP

    % Define sub-region of object:
    ynp = cen - yposn(np);
    xnp = cen + xposn(np);

    % Propagate exit surface wave to get diffraction pattern:
    esw = obj(ynp,xnp).*prb;
    ESWA_ideal(:,:,np) = abs(FT2(obj(ynp,xnp).*prb));

    % If desired, make images of ESWA and eswa at each position:
    if P.makeESWAIm == 1
        plotESWA(ESWA_ideal(:,:,np),abs(esw),np);
    end

    % For each wavelength in the bandwidth...
    for nw = 1:NW

        % Find x,y values to interpolate:
        xinterp = xd/alpha(nw);
        yinterp = yd/alpha(nw);

        % Interpolate and scale diffraction pattern intensities:
        I_interp = interp2(xd,yd,ESWA_ideal(:,:,np).^2,xinterp,yinterp);
        I_interp(isnan(I_interp)) = 0;
        I_interp = Sw(nw)*I_interp/(alpha(nw)^2);

        % Add intensities incoherently to form broadband pattern:
        ESWA_broad(:,:,np) = ESWA_broad(:,:,np) + I_interp;
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