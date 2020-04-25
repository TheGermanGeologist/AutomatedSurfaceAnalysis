function [fittedSurface,coeffs,rSquared] = ASAfitSurface(xGrid,yGrid,zData,order)

% reshape to column vectors
zVec = reshape(zData,[numel(zData), 1]);
xVec = reshape(xGrid,[numel(xGrid), 1]);
yVec = reshape(yGrid,[numel(yGrid), 1]);

% remove NaN values
x = xVec(~isnan(zVec));
y = yVec(~isnan(zVec));
z = zVec(~isnan(zVec));
clearvars xVec yVec zVec

switch order
    case 1
        % ax + by + c = 0
        equationSystem = [x y ones(size(z))];
        clearvars x y
        coeffs = equationSystem \ z;
        clearvars equationSystem
        fittedSurface = xGrid .* coeffs(1) + yGrid .*coeffs(2) + ...
              ones(size(xGrid),class(zData)) .* coeffs(3);
          
    case 2
        % ax^2 + 2bxy + cy^2 + d = 0
        equationSystem = [(x.^2) (2.*y.*x) (y.^2) ones(size(z))];
        clearvars x y
        coeffs = equationSystem \ z;
        clearvars equationSystem
        fittedSurface = xGrid.^2 .* coeffs(1) + yGrid.*xGrid .* coeffs(2) +...
            yGrid.^2 .* coeffs(3) + ones(size(xGrid),class(zData)) .* coeffs(4);
        
    case 3
        % ax^3 + 3bx^2y + 3cxy^2 + dy^3 + e = 0
        equationSystem = [(x.^3) (3.*x.^2.*y) (3.*x.*y.^2) (y.^3) ones(size(z))];
        clearvars x y
        coeffs = equationSystem \ z;
        clearvars equationSystem
        fittedSurface = xGrid.^3 .* coeffs(1) + xGrid.^2.*yGrid .* coeffs(2) + ...
            xGrid.*yGrid.^2 .* coeffs(3) + yGrid.^3 .* coeffs(4) + ...
            ones(size(xGrid),class(zData)) .* coeffs(5);
        
    case 4
        % ax^4 + 4bx^3y + 6cx^2y^2 + 4dxy^3 + ey^4 + f = 0
        equationSystem = [(x.^4) (4.*x.^3.*y) (6.*x.^2.*y.^2) (4.*x.*y.^3) (y.^4) ones(size(z))];
        clearvars x y
        coeffs = equationSystem \ z;
        clearvars equationSystem
        fittedSurface = xGrid.^4 .* coeffs(1) + 4.*xGrid.^3.*yGrid .* coeffs(2) + ...
            6.*xGrid.^2.*yGrid.^2 .* coeffs(3) + 4.*xGrid.*yGrid.^3 .* coeffs(4) + ...
            yGrid.^4 .* coeffs(5) + ones(size(xGrid),class(zData)) .* coeffs(6);
    case 5
        % ax^5 + 5bx^4y + 10cx^3y^2 + 10dx^2y^3 + 5exy^4 + fy^5 + g = 0
        equationSystem = [(x.^5) (5.*x.^4.*y) (10.*x.^3.*y.^2) (10.*x.^2.*y.^3) (5.*x.*y.^4) (y.^5) ones(size(z))];
        clearvars x y
        coeffs = equationSystem \ z;
        clearvars equationSystem
        fittedSurface = xGrid.^5 .* coeffs(1) + 5.*xGrid.^4.*yGrid .*coeffs(2) + ...
            10.*xGrid.^3.*yGrid.^2 .* coeffs(3) + 10.*xGrid.^2.*yGrid.^3 .* coeffs(4) + ...
            5.*xGrid.*yGrid.^4 .*coeffs(5) + yGrid.^5 .* coeffs(6) + ...
            ones(size(xGrid),class(zData)) .* coeffs(7);
        
    otherwise
        error('Order of fitted surface must be between 1 and 5')
    
end

if isa(coeffs,'gpuArray')
    coeffs = gather(coeffs);
else
end

% calculate R squared
% build fitted vector without NaN
fittedExNaN = fittedSurface(~isnan(zData));
clearvars zData
fittedVec = reshape(fittedExNaN,[numel(fittedExNaN) 1]);
clearvars fittedExNaN
% Chi squared
chiSquared = sum((fittedVec - z).^2);
clearvars fittedVec
% Chi zero squared
chi0squared = sum( (ones(size(z),class(fittedSurface)) .* mean(z) - z).^2 );

if isa(chi0squared,'gpuArray')
    rSquared = gather(1 - chiSquared / chi0squared);
else
    rSquared = 1 - chiSquared / chi0squared;
end

end
