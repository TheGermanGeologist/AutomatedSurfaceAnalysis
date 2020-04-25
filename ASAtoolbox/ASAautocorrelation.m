function [autocorrelation,betaMin,betaMinAngle] = ASAautocorrelation(shiftedTransformed,xGrid,yGrid)
%ASAautocorrelation calculates the ACF for a given fourier transformation
%DC component shifted to the center

% Calculation of the Autocorrelation (ACF) by using the Wiener Khinchin Theorem
% for Fourier transforms (Cohen, 1992: Convolution, Filtering, Linear System,
% the Wiener Khinchin Theorem: Generalizations)

% MATLAB applys the factor 1/(m*n) in the inverse fourier transform, but
% the fourier transorm is defined as applying the factor in the transform
% itself and not the inverse in most cases. Apply factor here so it doesn't
% get blown up by multiplying itself with the complex conjugate, then
% remoev it again below.
shiftedTransformed = shiftedTransformed ./ numel(shiftedTransformed);

% shift Fourier transform back to original position
transformed = ifftshift(shiftedTransformed);
% calculate ACF as (inverse) transform of the fourier transform multiplied by
% its complex conjugate

if ~isa(transformed,'gpuArray')
    autocorrelation = fftshift(ifft2(numel(transformed).*(transformed .* conj(transformed))));
else
    autocorrelation = fftshift(ifft2(numel(transformed).*(transformed .* conj(transformed)),'symmetric'));
end

% Autocorrelation parameterization
% find betaMin and the corresponding angle
[autoGridX,autoGridY] = meshgrid(xGrid(1,:)-max(xGrid(1,:))/2,...
    yGrid(:,1)-max(yGrid(:,1))/2);
x = autoGridX(autocorrelation < (max(max(autocorrelation)*0.1+max(max(autocorrelation)*0.005)))...
    & autocorrelation > (max(max(autocorrelation)*0.1-max(max(autocorrelation)*0.005))));
y = autoGridY(autocorrelation < (max(max(autocorrelation)*0.1+max(max(autocorrelation)*0.005)))...
    & autocorrelation > (max(max(autocorrelation)*0.1-max(max(autocorrelation)*0.005))));
distances = sqrt(x.^2+y.^2);
minimumpos = find(distances == min(distances));
betaMin = distances(minimumpos);
betaMinAngle = abs(atand(y(minimumpos) / x(minimumpos)));
if isa(autocorrelation,'gpuArray')
    betaMin = gather(betaMin);
    betaMinAngle = gather(betaMinAngle);
    % no point in returning ACF as gpuArray since it is used in subsequent
    % caculations
    autocorrelation = gather(autocorrelation);
else
end


end

