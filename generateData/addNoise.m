function ESWA_noise = addNoise(ESWA,P)
%ADDNOISE Adds detector (Gaussian) and shot (Poisson) noise to ESWA.
%==========================================================================

%% Parasitic Scattering
% Apply constant background in specified mask to diffraction patterns:
if P.useBackground == 1

    % If we don't have a region defined, draw it manually and save
    % (for example, try a triangle in the upper right corner):
    if ~exist('mask.mat','file')
        mask = double(roipoly(ESWA(:,:,1)));
        save('mask.mat','mask')
    else
        load('mask.mat','mask')
    end
    
    % Apply background to ESWA:
    ESWA = ESWA + (0.02*max(ESWA,[],'all'))*mask;
end

%% Shot Noise
% Apply Poisson noise to intensity (the number of photons at each pixel):
if P.useShotNoise == 1
    ESWA = sqrt(poissrnd(ESWA.^2));
end

%% Detector Noise
% Apply Gaussian noise to intensity:
if P.useDetectorNoise == 1
    electronNoise = ((P.QE).*abs(ESWA).^2) ...
        + (P.sigmaReadout).*randn(size(ESWA));
    adu = (P.sensitivity.*electronNoise) + (P.ADoffset);
    
    % Apply pixel saturation:
    if (P.useSaturation) == 1
        maxAdu = (2^(P.bitDepth)) - 1;
        adu(adu>maxAdu) = maxAdu;
    end
    
    adu = adu - (P.ADoffset);
    adu(adu<0) = 0;
    adu = double(adu);
    ESWA = sqrt( (adu./(P.sensitivity)) ./ (P.QE) );
end

ESWA_noise = ESWA;

end