% /* 
% * @Author: Athul Vijayan
% * @Date:   2014-08-17 22:17:21
% * @Last Modified by:   Athul Vijayan
% * @Last Modified time: 2014-08-17 23:50:34
% */

clear all
numGrids = 10;                       % number of different grids over which solution is sought
numFirstGridElems = 4; xVal = 1/2;  % number of elements in the first (coarsest) grid
xArray = cell(numGrids,1);          % array storing the x values at grid locations
phiFD = cell(numGrids,1);           % Forward difference solution vectors
midValsFD = zeros(numGrids,1);      % an array storing all phi values at x = 0.5, FD
stepVals = zeros(numGrids,1);       % array storing delta x values for different grids
phi_0 = 0;                          % phi b.c. at x = 0
phi_1 = 1;
k = -4*9*0.1/3^4*0.08;
%
%% Solve using Forward Differences
for i = 1:numGrids
    numElem = numFirstGridElems*2^(i-1);      % number of elements
    numPoints = numElem-1;                      % number of  grid points, N in class
    xArray{i} = linspace(0,1,numPoints);   % {} references Cell elements, () references array or matrix elements
    % xArray{i} = xArray{i}(2:numPoints+1);     % [0,1] divided into numPoints interior grid points
    deltaX = 1/numElem;                       % step size
    %
    diagVec = -2*ones(numPoints,1);              % main diagonal of length n
    myMat = gallery('tridiag',ones(numPoints-1,1) ,diagVec, ones(numPoints-1,1)); % matrix of coefficients
    for j=1:numPoints
        myVec(j) = k*deltaX^2*(3-j*deltaX)^2;
    end
    myVec(1) = myVec(1) - phi_0;
    myVec(numPoints) = myVec(numPoints) - phi_1;
    myVec = [phi_0 zeros(1,numPoints-2) phi_1]';
    phiFD{i} = myMat\myVec;                   % solution obtained by matrix inversion
    midValsFD(i) = phiFD{i}((numPoints+1)/2);           % pick the computed value at x = xVal; 
    stepVals(i) = deltaX;                     % store step size for plotting later

    trueVal = (-5/108)*(0.25^4-12*0.25^3 + 54*0.25^2 -126*0.25);
    errorAtX1(i) = 4/3*(trueVal - midValsFD(i));  % calculate the error in every iteration of delta x=1/3
    
    xVal = 1/2;
    trueVal = (-5/108)*(0.5^4-12*0.5^3 + 54*0.5^2 -126*0.5);
    index = find(xArray{i}==xVal);            % index returns the location of Xval in the array xArray{i}
    midValsFD(i) = phiFD{i}(index);           % pick the computed value at x = xVal; 
    errorAtX2(i) = 4/3*(trueVal - midValsFD(i));  % error in every iteration of delta x2/3

    xVal = 3/4;
    trueVal = (-5/108)*(0.75^4-12*0.75^3 + 54*0.75^2 -126*0.75);
    index = find(xArray{i}==xVal);            % index returns the location of Xval in the array xArray{i}
    midValsFD(i) = phiFD{i}(index);           % pick the computed value at x = xVal; 
    errorAtX3(i) = 4/3*(trueVal - midValsFD(i));  % error in every iteration of delta x2/3
end
 % To get error as a function of $$\Delta x$$ x we use polyfit.
plot(log(stepVals), log(errorAtX1))
p1 = polyfit(log(stepVals'), log(errorAtX1), 1)
p2 = polyfit(log(stepVals'), log(errorAtX2), 1)
p3 = polyfit(log(stepVals'), log(errorAtX3), 1)
