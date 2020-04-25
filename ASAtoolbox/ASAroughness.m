function [roughness] = ASAroughness(data,reference,spacing)
%ASAroughness calculates the roughness parameters
%SA, SQ, SSK and SKU using the reference surface provided

zDeviation = data - reference;

nancount = sum(sum(isnan(data)));
A = size(data,1) * size(data,2) * spacing^2;

% check if this is correct or if it skews the results...
% zDeviation(isnan(data)) = 0;
% actually there is no reason to preserve spatial information past this
% point so just kick out all NaN values
zDeviation = zDeviation(~isnan(zDeviation));

% arithmetical mean height SA
SA = sum(sum((abs(zDeviation)) .* spacing^2)) / (A - nancount * spacing^2);

% root mean square height distribution SQ
SQ = sqrt( sum(sum( (zDeviation .^2) .* spacing^2 )) / (A - nancount * spacing^2) );

% skewness SSK
SSK = (1/SQ^3) * (1/(A - nancount * spacing^2)) * sum(sum( (zDeviation .^3) .*spacing^2 ));

% kurtosis SKU
SKU = (1/SQ^4) * (1/(A - nancount * spacing^2)) * sum(sum( (zDeviation .^4) .*spacing^2 ));

if isa(data,'gpuArray')
    roughness.SA = gather(SA);
    roughness.SQ = gather(SQ);
    roughness.SSK = gather(SSK);
    roughness.SKU = gather(SKU);
else
    roughness.SA = SA;
    roughness.SQ = SQ;
    roughness.SSK = SSK;
    roughness.SKU = SKU;
end

end

