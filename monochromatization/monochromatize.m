function ESWA_mono = monochromatize(ESWA,P,rP)
%MONOCHROMATIZE Monochromatizes ESWA.
%==========================================================================

%% Load parameters:
N = P.N;
C = gpuArray(P.C);
type = rP.monoAlgo;
kmax = rP.kmax;
[~,~,NP] = size(ESWA);

% Initialize monochromatization:
ESWA_mono = ESWA;
b_all = ESWA;

%% Monochromatize
% For each scan position...
for np = 1:NP
    b = b_all(:,:,np);

    % Soft thresholding to reduce noise (optional):
    b(b<3) = 0;

    % Run selected algorithm to monochromatize:
    switch type
        case "BiCGSTAB"
            M = sqrt( BiCGSTAB(C,b(:).^2,kmax) );

        case "CGLS"
            reorth = 0;
            M = sqrt( cgls(C,b(:).^2,kmax,reorth) );
    end

    M = reshape(M,[N,N]);
    ESWA_mono(:,:,np) = M;
end

end