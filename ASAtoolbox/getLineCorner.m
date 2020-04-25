function [xEnd,yEnd,xCenter,yCenter] = getLineCorner(grid,angle)
%getLineCorner gets the x and y positions where a line of a given angle
%intersects the borders of a given grid

% This is probably a really convoluted and silly way to do this since there are
% a lot of edge cases which are caught in if statements
% possibly there are even more edge cases which might break the algorithm, but
% so far it seems to work

% get center point
% this is where the DC component is at
if mod(size(grid,2),2) == 0
    xCenter = (size(grid,2)/2) + 1;
else
    xCenter =(size(grid,2)+1)/2;
end
if mod(size(grid,1),2) == 0
    yCenter = (size(grid,1)/2) + 1;
else
    yCenter =(size(grid,1)+1)/2;
end
grid(yCenter,xCenter) = 1;

% Equations calculating the corresponding y or x coordinates
% for given x or y depending on the angle
lineEqX = @(angleVar,x) tand(angleVar) * x;
lineEqY = @(angleVar,y) cotd(angleVar) * y;

if angle >= 0 && angle <= 90
    % first try at getting the coordinates from the equations
	% but one of them is always out of bounds or inf
    yEnd = yCenter - lineEqX(angle,size(grid,2) -xCenter);
    xEnd = lineEqY(angle,size(grid,1) -yCenter);
    if yEnd > (size(grid,1)/2)+1 || yEnd < 1
		% intersecting at the top border or perfectly horizontal
        yEnd = 1;
        xEnd = xCenter + round(xEnd);
        
    elseif xEnd > size(grid,2)/2
		% intersecting at the right border 
        xEnd = size(grid,2);
        yEnd = round(yEnd);
        
    else
		% this was some sort of annyoing edge case
		% possibly when it hits the corner of the grid perfectly
        xEnd = xCenter + round(xEnd);
        yEnd = round(yEnd);
    end
    grid(yEnd,xEnd) = 1;
    xEnd(1) = xEnd;
    yEnd(1) = yEnd;
    
	% now the same thing for the lower / left corner intersection
    yEnd2 = round(abs(lineEqX(-angle,xCenter-1))) +yCenter;
    xEnd2 = round(lineEqY(-angle,yCenter)) + xCenter;
    if yEnd2 > size(grid,1)
		% intersecting at bottom border
        yEnd2 = size(grid,1);
        if xEnd2 < 1
			% for some reason we're now out of bonds to the left as well
			% some sort of edge case again
            xEnd2 = 1;
        else
        end
    elseif xEnd2 < 1
		% intersecting at left border
        if yEnd2 < 1
            yEnd2 = 1;
        else
        end
        xEnd2 = 1;
    else
    end
    
    grid(yEnd2,xEnd2) = 1;
    xEnd(2) = xEnd2;
    yEnd(2) = yEnd2;
    
    
elseif angle > 90 && angle <= 180
    
    yEnd = yCenter - lineEqX(angle,size(grid,2) -xCenter);
    xEnd = lineEqY(angle,size(grid,1) -yCenter);
    if yEnd > (size(grid,1)/2)+1
        xEnd = xCenter + round(xEnd);
        if xEnd < 1
            xEnd = 1;
            yEnd = size(grid,1)+1-round(yEnd);
        else
            yEnd = 1;
        end
        
    else
         if xEnd < 1
            xEnd = 1;
            yEnd = size(grid,1)+2-round(yEnd);
        else
            yEnd = 1;
        end
    end
    grid(yEnd,xEnd) = 1;
    xEnd(1) = xEnd;
    yEnd(1) = yEnd;
    
    yEnd2 = round(abs(lineEqX(-angle,xCenter-1))) +yCenter;
    xEnd2 = round(lineEqY(-angle,yCenter)) + xCenter;
    if xEnd2 > size(grid,2)
        xEnd2 = size(grid,2);

    else
    end
    if yEnd2 > size(grid,1)
        yEnd2 = size(grid,1);
        
    else
    end
    grid(yEnd2,xEnd2) = 1;
    xEnd(2) = xEnd2;
    yEnd(2) = yEnd2;
end

end

% honestly I have no idea what exactly I did here,
% it works for now, but it's probably a good idea
% to come up with a more transparent, less confusing solution