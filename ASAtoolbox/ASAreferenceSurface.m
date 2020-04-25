function [referenceSurface,coeffs,rSquared] = ASAreferenceSurface(xGrid,yGrid,zData,mode)
%ASAreferenceSurface returns a reference Surface for the roughness parameter calculation
%modes: 'mean': uses the mean elevation as reference. 'surf1' - 'surf5': uses a polynomial
%fit surface as reference. 'zero' uses 0 elevation as reference at all
%points.
%This is essentially a wrapper for ASAfitSurface and the special cases
%'mean' & 'zero'

switch mode
    case 'zero'
        if isa(zData,'gpuArray')
            referenceSurface = zeros(size(zData),'gpuArray');
        else
            referenceSurface = zeros(size(zData));
        end
        coeffs = [];
        rSquared = [];
    case 'mean'
        zVec = reshape(zData(~isnan(zData)),[numel(zData(~isnan(zData))),1]);
        if isa(zData,'gpuArray')
            referenceSurface = ones(size(zData),'gpuArray') .* mean(zVec);
        else
            referenceSurface = ones(size(zData)) .* mean(zVec);
        end
        coeffs = [];
        rSquared = [];
    case 'surf1'
        [referenceSurface,coeffs,rSquared] = ASAfitSurface(xGrid,yGrid,zData,1);
    case 'surf2'
        [referenceSurface,coeffs,rSquared] = ASAfitSurface(xGrid,yGrid,zData,2);
    case 'surf3'
        [referenceSurface,coeffs,rSquared] = ASAfitSurface(xGrid,yGrid,zData,3);
    case 'surf4'
        [referenceSurface,coeffs,rSquared] = ASAfitSurface(xGrid,yGrid,zData,4);
    case 'surf5'
        [referenceSurface,coeffs,rSquared] = ASAfitSurface(xGrid,yGrid,zData,5);
    otherwise
        error('invalid argument for reference value')
end
        

end

