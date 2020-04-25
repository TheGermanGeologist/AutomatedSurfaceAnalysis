function [topoCropped,lowpassedCropped,fixedCropped] = ASAcrop(lowpassed,topography,fixedData,fftThreshold,spacing)
%ASAcrop removes the borders of several topography arrays,
%depending on the threshold used to apply the low-pass to the data

% This cropping is done because at the threshold frequency, the borders
% will be incorrectly offset due to the periodicity of the DFT, which is
% not actually realized in the finite data.

% Calculate size of cropping area, this seems to yield reasonable results
cutoff = round((1/(8*fftThreshold))/spacing);

% apply cutoff to all three arrays
topoCropped = topography(cutoff:(end-cutoff),cutoff:(end-cutoff));
lowpassedCropped = lowpassed(cutoff:(end-cutoff),cutoff:(end-cutoff));
fixedCropped = fixedData(cutoff:(end-cutoff),cutoff:(end-cutoff));

end
