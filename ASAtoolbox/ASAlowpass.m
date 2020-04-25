function [lowpassed,dampened] = ASAlowpass(spectralMap,fx,fy,threshold)
%ASAlowpass returns a low-passed topography, given a fourier transform of the topography
%all frequencies above threshold will be dampened

% build array containing the spatial frequencies at the respective
% locations in the spectral map
[fxGrid,fyGrid] = meshgrid(fx,fy);
freqGrid = sqrt(fxGrid.^2 + fyGrid.^2);

% create dampening mask
% fade out higher frequencies exponentially
dampening = exp(1 - freqGrid./threshold);
% set multiplicator of frequencies lower than the threshold to 1
dampening(freqGrid < threshold) = 1;

% apply low-pass
dampened = spectralMap .* dampening;
% inverse fourier transform, creating the low-frequency topography
lowpassed = real(ifft2(ifftshift(dampened)));


end

