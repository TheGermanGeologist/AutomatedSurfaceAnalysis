function [untilted,varargout] = ASAremoveTilt(rawZdata,spacing,mode)
%ASAremoveTilt removes the (apparent) tilt from raw topography data
%using the mode specified. 'substract': substracts a 1st order fit plane
%from the data. 'rotate' applies a rotation matrix to the data, depending
%on the dip angles obtained from the 1st order fit surface

% build X & Y grid (disposable)
[xGrid,yGrid] = ASAbuildGrid(rawZdata,spacing);

% get first order plane
[fittedPlane,coeffs,~] = ASAfitSurface(xGrid,yGrid,rawZdata,1);
          
         
% actually removing the tilt now

switch mode
    case 'substract'
        % just substract the fitted plane from the raw data
        % this is biased against values at the higher end of the fitted
        % plane, but the coordinate system remains unchanged
        untilted = rawZdata - fittedPlane;
        % return spacing if required
        if nargout == 2
            newspacing.dx = spacing;
            newspacing.dy = spacing;
            varargout{1} = newspacing;
        else
        end
        
    case 'rotate'
        % this requires much more work
        % possible problems of this method: topography is not preserved
        % should the sample be at relatively steep angles, but actually
        % this is more of a theoretical concern. More important: the grid
        % spacing changes when rotating and the xy-axis aren't orthogonal
        % anymore. not really recommended at this point
		% it should be possible to come up with some sort of distortion correction
		% but it's complicated, confusing and might go wrong
		% In practice the change of the coordinate system is likely very small,
		% but it is kind of unpredictable
        warning('Using rotate for tilt correction is currently not recommended.')
        
        % first step: create map of good datapoints (not NaN)
        mapping = find(~isnan(rawZdata));
        
        % second step: shift origin to middle of xy plane
        yshift = -0.5 * max(max(yGrid));
        xshift = -0.5 * max(max(xGrid));
        yVec = yVec + yshift;
        xVec = xVec + xshift;
        
        % third step: get angles of plane
        xSlope = coeffs(1);
        yAngle = -atan(xSlope);

        ySlope = coeffs(2);
        xAngle = -atan(ySlope);
        
        % calculate roatation axis and roatation matrix
        M = [ 0, 0, -yAngle;...
              0, 0, -xAngle;...
             yAngle, xAngle, 0 ];
        R = expm(M);
        
        % finally: perform the roation
        planePoints = [xVec';yVec';zVec'];
        rotatedData = R * planePoints;
        
        % put everything back into arrays
        rotatedXgrid = zeros(size(rawZdata),class(rawZdata));
        rotatedXgrid(:,:) = NaN;
        rotatedXgrid(mapping) = rotatedData(1,:);
        rotatedYgrid = zeros(size(rawZdata),class(rawZdata));
        rotatedYgrid(:,:) = NaN;
        rotatedYgrid(mapping) = rotatedData(2,:);
        rotatedZgrid = zeros(size(rawZdata),class(rawZdata));
        rotatedZgrid(:,:) = NaN;
        rotatedZgrid(mapping) = rotatedData(3,:);
        % fixing the X Y Grids with linear interpolation should work fine
        % replace with inpaint_nans at some point, but the whole rotation
        % thing is broken anyways atm.
        rotatedYgrid = fillmissing(rotatedYgrid,'linear');
        rotatedXgrid = fillmissing(rotatedXgrid,'linear');
        
        % push everything out, let user deal with the rest outside
        untilted = rotatedZgrid;
        varargout{1} = spacing;         % FIX THIS!!!
        varargout{2} = rotatedXgrid;
        varargout{3} = rotatedYgrid;
       
        
    otherwise
        error('mode must either be ''substract'' or ''rotate''')
end


          
end