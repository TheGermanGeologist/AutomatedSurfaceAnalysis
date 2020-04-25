function [apsd] = ASAapsd(transformedshifted,angleBins,spacing)
%ASAapsd calculates the Angular Power Spectral Density for given
%angle bins: 0 -> (180-stepsize), must be uniformly increasing
%e.g. 0:10:170 or 0:5:175
%for the input fourier transform (DC shifted to the center)

% calculate Power Spectral Density
% as defined in Cohen, 1992 (Convolution, Filtering, Linear
% Systems, the Wiener Khinchin Theorem: Generalizations
PSD = abs(transformedshifted).^2;

% The Power Spectral Density is scaled by dividing by the area of the
% domain.
A = numel(transformedshifted) * spacing^2;
PSD = PSD ./ A;


% Line Map containing all lines, used to calculate weight of points falling
% on a line separating two areas
lineMap = ones(size(PSD),class(PSD));

% set up borders of the angle bins
binNo = numel(angleBins);
angleRange = (angleBins(2) - angleBins(1)) / 2;
lineAngles = angleBins + angleRange;
lineAnglesShifted = [lineAngles(end), lineAngles(1:end-1)];

% create line map
for ii = 1:binNo
    lineAngles(ii);
    tempLines = zeros(size(PSD),class(PSD));
    [xEnd,yEnd,xCenter,yCenter] = getLineCorner(lineMap,lineAngles(ii));
    [yline,xline] = bresenham(xCenter,yCenter,xEnd(1),yEnd(1));
    idx = sub2ind(size(PSD),xline,yline);
    tempLines(idx) = 1;
    [yline,xline] = bresenham(xCenter,yCenter,xEnd(2),yEnd(2));
    idx = sub2ind(size(PSD),xline,yline);
    tempLines(idx) = 1;
    lineMap = lineMap + tempLines; 
end

% apply weight to PSD
PSDweighed = PSD ./ lineMap;

clearvars tempLines PSD

% get angle slices of PSD, calculate APSD
apsd = zeros(1,numel(angleBins),class(PSDweighed));

for ii = 1:binNo
    % create logical array of the current slice
    currSlice = zeros(size(PSDweighed),class(PSDweighed));
    % get two lines
    % line 1
    [xEnd1,yEnd1,xCenter,yCenter] = getLineCorner(lineMap,lineAngles(ii));
    [yline1,xline1] = bresenham(xCenter,yCenter,xEnd1(1),yEnd1(1));
    idx1 = sub2ind(size(PSDweighed),xline1,yline1);
    currSlice(idx1) = 1;
    [yline1,xline1] = bresenham(xCenter,yCenter,xEnd1(2),yEnd1(2));
    idx1 = sub2ind(size(PSDweighed),xline1,yline1);
    currSlice(idx1) = 1;
    % line 2
    [xEnd2,yEnd2,xCenter,yCenter] = getLineCorner(lineMap,lineAnglesShifted(ii));
    [yline2,xline2] = bresenham(xCenter,yCenter,xEnd2(1),yEnd2(1));
    idx2 = sub2ind(size(PSDweighed),xline2,yline2);
    currSlice(idx2) = 1;
    [yline2,xline2] = bresenham(xCenter,yCenter,xEnd2(2),yEnd2(2));
    idx2 = sub2ind(size(PSDweighed),xline2,yline2);
    currSlice(idx2) = 1;
    % fill area
    % get seed point
    [xSeed,ySeed,~,~] = getLineCorner(lineMap,angleBins(ii));
    filled = imfill(logical(currSlice),[ySeed(1) xSeed(1)],4);
    filled = imfill(filled,'holes');
    filled = imfill(logical(filled),[ySeed(2) xSeed(2)],4);
    filled = imfill(filled,'holes');
    
    
    apsd(ii) = real(sum(sum(PSDweighed(logical(filled)))));
    
end
if isa(apsd,'gpuArray')
   apsd = gather(apsd); 
else
end
end

