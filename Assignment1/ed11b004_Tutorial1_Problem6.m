% /* 
% * @Author: Athul Vijayan
% * @Date:   2014-08-17 20:17:08
% * @Last Modified by:   Athul Vijayan
% * @Last Modified time: 2014-08-17 22:16:30
% */

%% Initialize arrays and allocate space
clear all
numGrids = 10;                       % number of different grids over which solution is sought
numFirstGridElems = 4; xVal = 1/2;  % number of elements in the first (coarsest) grid
xArray = cell(numGrids,1);          % array storing the x values at grid locations
phiFD = cell(numGrids,1);           % Forward difference solution vectors
midValsFD = zeros(numGrids,1);      % an array storing all phi values at x = 0.5, FD
stepVals = zeros(numGrids,1);       % array storing delta x values for different grids
phi_0 = 0;                          % phi b.c. at x = 0
phi_1 = 1;
%
%% Solve using Forward Differences
for i = 1:numGrids
    numElem = numFirstGridElems*2^(i-1);      % number of elements
    numPoints = numElem-1;                      % number of  grid points, N in class
    xArray{i} = linspace(0,1,numPoints);   % {} references Cell elements, () references array or matrix elements
    % xArray{i} = xArray{i}(2:numPoints+1);     % [0,1] divided into numPoints interior grid points
    deltaX = 1/numElem;                       % step size
    %
    diagVec = (2 + deltaX^2)*ones(numPoints,1);              % main diagonal of length n
    myMat = gallery('tridiag',-1*ones(numPoints-1,1) ,diagVec, -1*ones(numPoints-1,1)); % matrix of coefficients
    myVec = [phi_0 zeros(1,numPoints-2) phi_1]';
    phiFD{i} = myMat\myVec;                   % solution obtained by matrix inversion
    midValsFD(i) = phiFD{i}((numPoints+1)/2);           % pick the computed value at x = xVal; 
    stepVals(i) = deltaX;                     % store step size for plotting later

end
x=0.5;
trueVal = (exp(1)/(exp(1)^2 -1 ))*(exp(x)-exp(-x));
errorAtX1 = midValsFD - trueVal;
%% plot FD fit and true solution

power = polyfit(log(stepVals), log(errorAtX1), 1)
figure;
plot(log(stepVals), log(abs(errorAtX1')), 'bo','MarkerSize',12, ...
                                            'MarkerFaceColor',[0.8 0.8 0.8], ...
                                            'MarkerEdgeColor','r',...
                                            'LineWidth',2)
annotation('textbox', [.2 .7 .1 .1], 'String', ['slope using polyfit ',num2str(power(1))]);

fitVal = polyval(power,log(stepVals)); % obtain straight line fit to data
hold on % set hold = on to make next plot appear in the same figure
plot(log(stepVals),fitVal,'k:','LineWidth',2) % plot the straight line fit, 'k' for black, ':' for dotted line
%
% set various plot options
set(gca,'XTick',-8:2:0) % tick marks on x-axis
set(gca,'YTick',-24:6:-6) % tick marks on y-axis
set(gca,'FontSize',14) % font size for both axes
xlim([-8,0]); % x-range for plot
ylim([-24 -6]); % y-range for plot
xlabel('$$\log\Delta x$$',... % Label for x-axis
 'Interpreter','LaTeX', ...
 'FontSize',18)
ylabel('$$\log E$$',... % Label for y-axis
 'Interpreter','LaTeX', ...
 'FontSize',18)
grid('on') % grid lines to read off values easier
legend('Finite Difference','Linear Fit','Location','SouthEast') % plot legend
title('Error in $$\phi_{x=0.5}$$ decreases as the square of step size $$\Delta x$$', ...
 'Interpreter','LaTeX','FontSize',18)
