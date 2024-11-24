function [xpos,ypos] = makeScanGrid(P)
%MAKESCANGRID_SQUARE Makes a square scan grid with random offset and
%returns the physical coordinates of each position.
%==========================================================================

%% Load Parameters
NP = P.NP;

%% Make Scan Grid

% Calculate the distance between scan points:
ds = 2*P.prbR*(1-P.probeOverlap);

% Number of scans in X and Y:
nScansXY = sqrt(NP);

% Initialize position vectors [pixels]:
[xpos,ypos] = deal(zeros(1,NP));

% Calculate positions:
posInd = 1;
for nx = 1:nScansXY

    % X position [m]:
    x = ds*nx;

    for ny = 1:nScansXY

        % Y position [m]:
        y = ds*ny;

        % Calculate total offset:
        offset = (ds*P.randOffsetFrac)*exp(1i*rand);

        % Find position and plot circle around it (for probe size):
        xpos(posInd) = x + real(offset);
        ypos(posInd) = y + imag(offset);

        % Increment index:
        posInd = posInd+1;
    end
end

% Put center at the origin:
xpos = xpos - mean(xpos);
ypos = ypos - mean(ypos);

P.xpos = xpos;
P.ypos = ypos;

end
