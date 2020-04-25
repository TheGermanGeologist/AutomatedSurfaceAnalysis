function [spectralMap,fx,fy] = ASAfft2(fixedData,spacing)
%ASAfft2 calculates the 2D fourier transform for a given topography
%dataset. Topography can't contain NaN, must be interpolated.
%The input data should ideally be square.
%also returns fx and fy, vectors representing the spatial frequency of
%the DC-centered fourier transfrom


% Padding for faster FFT performance
% actually,the slight performance boost is not worth all the trouble
% padding makes dampening unnessecaryly complicated
% and the massively bigger arraysize is likely to exceed the memory

% xPadding = 2^nextpow2(size(fixedData,2));
% yPadding = 2^nextpow2(size(fixedData,1));
% spectralMap = fftshift(fft2(fixedData,xPadding,yPadding));

% 2D FFT and shifting of zero frequency into the middle
% there probably should be a factor of 1/(m*n) here so the autocorrelation,
% PSD etc are correct. Check with people who know!
spectralMap = fftshift(fft2(fixedData));

% the x and y axes of the spectral map as spatial frequency [cycles/mm]
if isa(fixedData,'gpuArray')
    fx = gpuArray(1/(spacing*size(fixedData,2))).*((-size(fixedData,2)/2:size(fixedData,2)/2-1));
    fy = gpuArray(1/(spacing*size(fixedData,1))).*((-size(fixedData,1)/2:size(fixedData,1)/2-1));
else
    fx = (1/(spacing*size(fixedData,2))).*((-size(fixedData,2)/2:size(fixedData,2)/2-1));
    fy = (1/(spacing*size(fixedData,1))).*((-size(fixedData,1)/2:size(fixedData,1)/2-1));
end

end
