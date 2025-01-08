function C = generateC(P)
%GENERATEC Generates sparse matrix C satisfying b=Cm.
%==========================================================================
% Here, b is the broadband diffraction pattern with size N^2 x 1, with N
% being the number of pixels along x or y. m is the monochromatized
% diffraction pattern with size N^2 x 1.

%% Load Parameters
Sw = P.Sw;
ds_det = P.ds_det;
N = P.N;
Np = N^2;
alpha = P.alpha;
NW = P.NW;

%% Initialize
% Find x and y (2D):
xvec = ds_det*( -N/2 : (N/2-1) );
yvec = xvec;
[x,y] = ndgrid(xvec,yvec);

xn = [];
yn = [];
cnj = [];

%% Calculate C
% For each wavelength...
parfor nw = 1:NW

    fprintf('Generating C: Wavelength %d of %d\n',nw,NW)

    % Create 3D interpolation matrix:
    m_int = zeros(N,N,Np);
    for nxy = 1:Np
        [nx,ny] = ind2sub(size(x),nxy);
        m_int(nx,ny,nxy) = 1;
    end

    % Interpolate to find weights (NxNxNp):
    weights = interpn(x,y,m_int,x/alpha(nw),y/alpha(nw))*Sw(nw)/alpha(nw)^2;
    weights(isnan(weights)) = 0;

    % For each position...
    for nxy = 1:Np

        % Find non-zero weights as linear indices:
        [row,col,v] = find(weights(:,:,nxy));
        lin = sub2ind(size(weights(:,:,nxy)),row,col);

        % Convert to indices of column of C:
        xn = [xn lin'];
        yn = [yn nxy*ones(size(lin))'];
        cnj = [cnj v'];
    end
end

% Create sparse matrix:
C = sparse(xn,yn,cnj,Np,Np);
end