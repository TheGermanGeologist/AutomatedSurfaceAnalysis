function [windowedData,cutoffSize] = ASAfourierWindowing(zData)
%ASAfourierWindowing applies windowing to an input topography dataset
%in order to avoid creating vertical and horizontal artifacts in the DFT

% get size of fade out range (cutoff size)
shortSide = min(size(zData));
cutoffSize = round(shortSide / 10);

% set up grid and vectors
xWindow = ones(1,size(zData,2),class(zData));
yWindow = ones(1,size(zData,1),class(zData));
if isa(zData,'gpuArray')
    fade = cosd(gpuArray.linspace(-90,0,cutoffSize));
else
    fade = cosd(linspace(-90,0,cutoffSize));
end

% create window
xWindow(1,1:cutoffSize)= fade;
xWindow(1,(end-(cutoffSize-1)):end) = flip(fade);
yWindow(1,1:cutoffSize)= fade;
yWindow(1,(end-(cutoffSize-1)):end) = flip(fade);
[winXgrid,winYgrid] = meshgrid(xWindow,yWindow);
windowMask = winXgrid .* winYgrid;

% apply window
windowedData = zData .* windowMask;

end

