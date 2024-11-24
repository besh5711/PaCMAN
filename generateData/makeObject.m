function [obj,P] = makeObject(P)
%MAKEOBJECT Makes object with a "cameraman" transmission and "peppers"
%phase. The object size is 200x200 pixels.
%==========================================================================

%% Load Parameters
N = P.N;
xposn = P.xposn;
tMin = P.tMin;
pMax = P.pMax;

%% Calculate Full Object Size
% Find number of pixels in object reconstruction based on probe pixel
% dimension and positions:
Nobj = N + (max(xposn)-min(xposn));
Nobj = round(Nobj*1.1);

% Round to even number:
if rem(Nobj,2) == 1; Nobj = Nobj + 1; end

% Save parameters:
P.Nobj = Nobj;

%% Load Amplitude
% Make the amplitude the normalized image of "cameraman":
target_amp = imread('cameraman.tif');
target_amp = double(target_amp(21:220,21:220));     % crop
target_amp = target_amp - min(target_amp,[],'all'); % set min to 0
target_amp = target_amp/max(target_amp(:));         % set max to 1
target_amp = tMin + (1-tMin)*target_amp;

%% Load Phase
% Make the phase the normalized image of "peppers":
target_phase = imread('peppers.png');
target_phase = double(target_phase(121:320,101:300));
target_phase = target_phase - min(target_phase,[],'all');   % set min to 0
target_phase = pMax*target_phase/max(target_phase(:));

%% Create Object
% Find size and coordinates of target:
[tarN,~] = size(target_amp);
cen = (-tarN/2 : tarN/2-1) + Nobj/2;
P.cen = cen;

% Create object with a transmission that varies from tMin to tBack, where
% tBack is the transmission of the background:
obj = ones(Nobj,Nobj);
obj(cen,cen) = target_amp;

% Apply the phase:
obj(cen,cen) = obj(cen,cen).*exp(1i*target_phase);

end