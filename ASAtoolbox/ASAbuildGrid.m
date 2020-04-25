function [xGrid,yGrid] = ASAbuildGrid(rawZdata,spacing)
%ASAbuilGrid creates X and Y meshgrids from a given raw data and its
%spacing. used in multiple instances.
if isa(rawZdata,'gpuArray')
    x = gpuArray(( 0:(size(rawZdata,2)-1) ) .* spacing);
    y = gpuArray(( 0:(size(rawZdata,1)-1) ) .* spacing);
else
    x = ( 0:(size(rawZdata,2)-1) ) .* spacing;
    y = ( 0:(size(rawZdata,1)-1) ) .* spacing;
end

[xGrid,yGrid] = meshgrid(x,y);

end

