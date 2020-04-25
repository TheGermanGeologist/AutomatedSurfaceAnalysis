function [xGridFixed,yGridFixed,topographyFixed] = ASAfixDimensions(rawData,metaData)
%ASAfixDimensions converts raw data to height values in micrometers
%shifts the lowest point of the data to zero
%returns X- and Y-grids, containing the spatial dimensions of the data in mm

% get factors
spacing = metaData.pixelsize;
heightfactor = metaData.wavelength / 1e3;

% build X Y grid with appropriate dimensions [mm]
[xGridFixed,yGridFixed] = ASAbuildGrid(rawData,spacing);

% convert raw data to height in µm
topographyScaled = rawData .* heightfactor;
% shift lowest point to zero
minimum = min(min(topographyScaled));
topographyFixed = topographyScaled - minimum;

end

