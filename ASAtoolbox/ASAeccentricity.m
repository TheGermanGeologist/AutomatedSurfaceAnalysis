function [Xpoints,Ypoints,ellipses,eccentr,truedir] = ASAeccentricity(gamma)
%ASAeccentricity fit three ellipses to isolines at equidistant positions
%in the input 2D-semivariogram. The distances are fixed to 25, 50 and 75%
%of the distance from edge to center in x-direction 

% check if the range is correct (must be even)
% the algorithm could be adapted for uneven ranges, but it's
% extra work and using the next even number for range doesn't hurt anyone
range = length(gamma);
if mod(range,2) == 0
    error('Range of semivariogram must be even')
end

% get the positions at which we want to sample for isolines
pos1 = round(0.375 * range);
pos2 = round(0.25 * range);
pos3 = round(0.125 * range);

midline = round(0.5 * range);

% get the values of these isolines
value(1) = gamma(midline,pos1);
value(2) = gamma(midline,pos2);
value(3) = gamma(midline,pos3);

% define range in which we want to look for points corresponding to the isoline
% this would be even better if it was dynamically calculated depending on the spacing of
% the values obtained above
padding = 10;
midline = midline -1;

[xGrid,yGrid] = meshgrid(-midline:midline,-midline:midline);

for i = 1:3
	% get coordinates for points associated with our isoline
    x = xGrid(gamma < (value(i)+padding) & gamma > (value(i)-padding));
    y = yGrid(gamma < (value(i)+padding) & gamma > (value(i)-padding));
    
	% store coordinates in cell for output
    Xpoints{i} = x;
    Ypoints{i} = y;
    
	% fit ellipse, using third party function by Ohad Gal
    ellipses{i} = fit_ellipse(x,y);
    
	% calculate the eccentricity and store it for output
	% epsilon = e/a, where e = (a^2 + b^2)^(1/2)
%     eccentr(i) = sqrt(ellipses{i}.short_axis^2 + ellipses{i}.long_axis^2)/...
%         ellipses{i}.long_axis;
    % new approach:
    eccentr(i) = sqrt(1 - ellipses{i}.short_axis.^2 / ellipses{i}.long_axis.^2);
    % figure out in which direction the ellipse points
    if ellipses{i}.xRadius < ellipses{i}.yRadius && ellipses{i}.phi > 0
        % direction left
        truedir(i) = rad2deg(ellipses{i}.phi) + 90;
    elseif ellipses{i}.xRadius < ellipses{i}.yRadius && ellipses{i}.phi < 0
        % direction right
        truedir(i) = rad2deg(ellipses{i}.phi) + 90;
    elseif ellipses{i}.xRadius > ellipses{i}.yRadius && ellipses{i}.phi > 0
        % direction right
        truedir(i) = rad2deg(ellipses{i}.phi);
    elseif ellipses{i}.xRadius > ellipses{i}.yRadius && ellipses{i}.phi < 0
        % direction left
        truedir(i) = rad2deg(ellipses{i}.phi) + 180;
    elseif ellipses{i}.xRadius > ellipses{i}.yRadius && ellipses{i}.phi == 0
        truedir(i) = 0;
    elseif ellipses{i}.xRadius < ellipses{i}.yRadius && ellipses{i}.phi == 0
        truedir(i) = 90;
    elseif rad2deg(ellipses{i}.phi) == 45 || rad2deg(ellipses{i}.phi) == -45
        warning(['Ellipse no. ',num2str(i),' returned an angle of exactly 45 degree'])
        warning('Ellipse algorithm is unstable at this edge case, might be off by 90 degree')
        truedir(i) = 45;
    elseif ellipses{i}.long_axis == ellipses{i}.shortaxis
        warning(['Ellipse no. ',num2str(i),' appears to be a circle.'])
        truedir(i) = 0;
    else
        warning(['Unable to determine ellipse direction for ellipse no. ',num2str(i)])
        truedir(i) = NaN;
    end
end


end