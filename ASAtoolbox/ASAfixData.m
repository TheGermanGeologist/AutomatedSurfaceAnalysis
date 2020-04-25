function [fixedData,percentOutl,percentNaN] = ASAfixData(rawData,options,varargin)
%ASAfixData removes outliers from the raw data, then fills NaN using interpolation
%takes optional argument "windowSize", containing custom value (integers)
%for the outlier detection window size

if nargin == 3
    outlWin = varargin{2};
elseif nargin > 3
    error('Max. number of arguments is 2.')
else
    % default window sizes
    outlWin = 20;
end

% outlier detection doesn't profit from GPU and 3rd party interpolation fcn
% is incompatible, so bring everything to CPU for the moment

if isa(rawData,'gpuArray')
    rawData = gather(rawData);
    gpuSwitcheroo = 1;
else
    gpuSwitcheroo = 0;
end

if strcmp(options.outlierDetection,'yes')
    % Remove outliers using a moving median window as detection and nearest neighbour as interpolation
    % columns
    outlierPosC = isoutlier(rawData,'movmedian',outlWin,1);
    % rows
    outlierPosR = isoutlier(rawData,'movmedian',outlWin,2);
    % combined
    outliersAdded = outlierPosC + outlierPosR;
    outliersCombined = outliersAdded ~= 0;
    % calculate percentage of outliers
    percentOutl = (sum(sum(outliersCombined)) / numel(rawData)) * 100;
    % remove Outliers
    outliersRem = rawData;
    outliersRem(outliersCombined) = NaN;
else
    percentOutl = 0;
    outliersRem = rawData;
end



% caculate percentage of bad data
percentNaN = (sum(sum(isnan(rawData))) / numel(rawData)) * 100;

% new version using 3rd party function for NaN inpainting by John D'Errico
fixedData = inpaint_nans(outliersRem);

if gpuSwitcheroo
    fixedData = gpuArray(fixedData);
else
end

