function [gamma,sill] = ASAsemivariogram2D(rawData,range)

%Developed by Hannes Krietsch (2014)

range = round(range);
if mod(range,2) ~= 0
    warning('Range of semivariogram is not even, rounding down to next even number')
    range = range-1;
end
sill = var(rawData(~isnan(rawData)));

gamma = zeros(2*range+1,class(rawData));

for i = 0:range
    tic
    for j = 0:range
        rawDatadiff = rawData(i+1:end,j+1:end) - rawData(1:end-i,1:end-j);
        gamma(i+range+1,j+range+1) = 0.5 * mean(rawDatadiff(~isnan(rawDatadiff)).^2);
        gamma(-i+range+1,-j+range+1) = gamma(i+range+1,j+range+1);
    end

    for j = 1:range
        rawDatadiff = rawData(i+1:end,1:end-j) - rawData(1:end-i,j+1:end);
        gamma(i+range+1,-j+range+1) = 0.5 * mean(rawDatadiff(~isnan(rawDatadiff)).^2);
        gamma(-i+range+1,j+range+1) = gamma(i+range+1,-j+range+1);
    end
    t = toc
    tRemaining = (range-i) * t;
    disp(['Estimated remaining time: ',num2str(tRemaining),' s'])
end


end

