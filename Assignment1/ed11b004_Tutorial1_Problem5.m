% @Author: Athul Vijayan
% @Date:   2014-08-16 11:45:06
% @Last Modified by:   Athul Vijayan
% @Last Modified time: 2014-08-16 11:46:13

% We will solve the following equation using finite-differences:
%
% $$ \frac{d \phi}{d x} = 3 \phi + 4 $$, 
% 
% subject to the boundary conditon
%
% $$\phi = 0 \ \mathrm{at} \ x=0$$
%
% The solution is 
%
% $$\phi = \frac{4}{3}(e^(3 x) -1)$$

% we use 20 grid refinements to study error at x=\frac{1}{2} and x=\frac{2}{3}.
% $\Delta x = [\frac{1}{3}, \frac{1}{6}, \frac{1}{9},......, \frac{1}{60}]$
clear all
numGrids = 10;                       % number of different grids over which solution is sought
numFirstGridElems = 3; xVal = 1/3;  % number of elements in the first (coarsest) grid
xArray = cell(numGrids,1);          % array storing the x values at grid locations
phiFD = cell(numGrids,1);           % solution vectors
midValsFD = zeros(numGrids,1);      % an array storing all phi values at x , FD
stepVals = zeros(numGrids,1);       % array storing delta x values for different grids
phi_0 = 0;                          % phi b.c. at x = 0

for i = 1:numGrids
    numElem = numFirstGridElems*2^(i-1);      % number of elements
    numPoints = numElem;                      % number of  grid points, N in class
    xArray{i} = linspace(0,1,numPoints+1)';   % {} references Cell elements, () references array or matrix elements
    xArray{i} = xArray{i}(2:numPoints+1);     % [0,1] divided into numPoints interior grid points
    deltaX = 1/numElem;                       % step size
    %
    diagVec = ones(numPoints,1);              % main diagonal of length n
    offdiagVec = -1.0*(3*deltaX + 1)*ones(numPoints-1,1);    % off-diagonal vectors of length (n-1)
    myMat = gallery('tridiag',offdiagVec,diagVec,zeros(numPoints-1,1)); % matrix of coefficients
    myVec = 4*deltaX*ones(numPoints,1);       % rhs vector
    myVec(1) = myVec(1) + phi_0*(3*deltaX + 1);
    phiFD{i} = myMat\myVec;                   % solution obtained by matrix inversion
    index = find(xArray{i}==xVal);            % index returns the location of Xval in the array xArray{i}
    midValsFD(i) = phiFD{i}(index);           % pick the computed value at x = xVal; 
    stepVals(i) = deltaX;                     % store step size for plotting later

    errorAtX1(i) = midValsFD(i)-4/3*(exp(3*xVal)-1);  % calculate the RMS error in every iteration of delta x=1/3
    
    xVal = 2/3;
    index = find(xArray{i}==xVal);            % index returns the location of Xval in the array xArray{i}
    midValsFD(i) = phiFD{i}(index);           % pick the computed value at x = xVal; 
    errorAtX2(i) = midValsFD(i)-4/3*(exp(3*xVal)-1);  % error in every iteration of delta x2/3
end
% Now we have error at x1 and x2 for different delta x. we plot log(error) Vs log(delta x) for  checking the order of error
%% plot FD fit and true solution
power = polyfit(log(abs(errorAtX1')), log(stepVals), 1)
figure;
plot(log(stepVals), log(abs(errorAtX1')), 'g')
title('Relationship betweeen $$log(Error)$$ at x=$$1/3$$ Vs $$log(\Delta x)$$', ...
    'Interpreter','LaTeX','FontSize',18)       % plot title at the top
xlabel('$$log(\Delta x)$$',...                 % Label for x-axis
        'Interpreter','LaTeX', ...
        'FontSize',18)
ylabel('$$log(Error)$$',...           % Label for y-axis
        'Interpreter','LaTeX', ...
        'FontSize',18)
annotation('textbox', [.2 .7 .1 .1], 'String', ['slope using polyfit ',num2str(power(1))]);
power = polyfit(log(abs(errorAtX2')), log(stepVals), 1);
figure;
plot(log(stepVals), log(abs(errorAtX2')), 'g')
title('Relationship betweeen $$log(Error)$$ at x=$$2/3$$ Vs $$log(\Delta x)$$', ...
    'Interpreter','LaTeX','FontSize',18)       % plot title at the top
xlabel('$$log(\Delta x)$$',...                 % Label for x-axis
        'Interpreter','LaTeX', ...
        'FontSize',18)
ylabel('$$log(Error)$$',...           % Label for y-axis
        'Interpreter','LaTeX', ...
        'FontSize',18)
annotation('textbox', [.2 .7 .1 .1], 'String', ['slope using polyfit ',num2str(power(1))]);

%% plot FD fit and true solution

figure;
trueVec = 4*(exp(3*xArray{i})-1)/3;
plot(xArray{i},trueVec,'k:','LineWidth',2);
hold on
i=1;
plot(xArray{i},phiFD{i},'bo','MarkerSize',6, ...
                           'MarkerFaceColor',[0.8 0.8 0.8], ...
                           'MarkerEdgeColor','r',...
                           'LineWidth',2);
i=10;
plot(xArray{i},phiFD{i},'bo','MarkerSize',6, ...
                           'MarkerFaceColor',[0.8 0.8 0.8], ...
                           'MarkerEdgeColor','g',...
                           'LineWidth',2);
% set various plot options
set(gca,'XTick',0:1/3:1.0)         % tick marks on x-axis
set(gca,'YTick',0:10:30)           % tick marks on y-axis
set(gca,'FontSize',14)             % font size for both axes
xlim([0 1.0]);                     % x-range for plot
ylim([0 30]);                      % y-range for plot
xlabel('$$x$$',...                 % Label for x-axis
        'Interpreter','LaTeX', ...
        'FontSize',18)
ylabel('$$\phi(x)$$',...           % Label for y-axis
        'Interpreter','LaTeX', ...
        'FontSize',18)
grid('on')                         % grid lines to read off values easier
legend('Exact','Finite Difference at coarsest','Finite Difference at finest','Location','NorthWest') % plot legend
title('Comparison of Finite Difference Solution for coarsest and finest mesh to Exact Solution ', ...
    'Interpreter','LaTeX','FontSize',18)        % plot title at the top

clear all
numGrids = 10;                       % number of different grids over which solution is sought
numFirstGridElems = 3; xVal = 1/3;  % number of elements in the first (coarsest) grid
xArray = cell(numGrids,1);          % array storing the x values at grid locations
phiFD = cell(numGrids,1);           % solution vectors
midValsFD = zeros(numGrids,1);      % an array storing all phi values at x , FD
stepVals = zeros(numGrids,1);       % array storing delta x values for different grids
phi_0 = 0;                          % phi b.c. at x = 0

% Now we do the above solution using backward difference formula for first derivative
for i = 1:numGrids
    numElem = numFirstGridElems*2^(i-1);      % number of elements
    numPoints = numElem;                      % number of  grid points, N in class
    xArray{i} = linspace(0,1,numPoints+1)';   % {} references Cell elements, () references array or matrix elements
    xArray{i} = xArray{i}(2:numPoints+1);     % [0,1] divided into numPoints interior grid points
    deltaX = 1/numElem;                       % step size

    offdiagVec = -1*ones(numPoints-1,1);              % main diagonal of length n
    diagVec = (1 - 3*deltaX)*ones(numPoints,1);    % off-diagonal vectors of length (n-1)
    myMat = gallery('tridiag',offdiagVec,diagVec,zeros(numPoints-1,1)); % matrix of coefficients
    myVec = 4*deltaX*ones(numPoints,1);       % rhs vector
    myVec(1) = myVec(1) + phi_0*(3*deltaX + 1);
    phiFD{i} = myMat\myVec;                   % solution obtained by matrix inversion
    index = find(xArray{i}==xVal);            % index returns the location of Xval in the array xArray{i}
    midValsFD(i) = phiFD{i}(index);           % pick the computed value at x = xVal; 
    stepVals(i) = deltaX;                     % store step size for plotting later

    errorAtX1(i) = midValsFD(i)-4/3*(exp(3*xVal)-1);  % calculate the RMS error in every iteration of delta x=1/3
    
    xVal = 2/3;
    index = find(xArray{i}==xVal);            % index returns the location of Xval in the array xArray{i}
    midValsFD(i) = phiFD{i}(index);           % pick the computed value at x = xVal; 
    errorAtX2(i) = midValsFD(i)-4/3*(exp(3*xVal)-1);  % error in every iteration of delta x2/3
end

% Now we have error at x1 and x2 for different delta x. we plot log(error) Vs log(delta x) for  checking the order of error
%% plot FD fit and true solution
power = polyfit(log(abs(errorAtX1')), log(stepVals), 1)
figure;
plot(log(stepVals), log(abs(errorAtX1')), 'g')
title('Relationship betweeen $$log(Error)$$ at x=$$1/3$$ Vs $$log(\Delta x)$$', ...
    'Interpreter','LaTeX','FontSize',18)       % plot title at the top
xlabel('$$log(\Delta x)$$',...                 % Label for x-axis
        'Interpreter','LaTeX', ...
        'FontSize',18)
ylabel('$$log(Error)$$',...           % Label for y-axis
        'Interpreter','LaTeX', ...
        'FontSize',18)
annotation('textbox', [.2 .7 .1 .1], 'String', ['slope using polyfit ',num2str(power(1))]);
power = polyfit(log(abs(errorAtX2')), log(stepVals), 1);
figure;
plot(log(stepVals), log(abs(errorAtX2')), 'g')
title('Relationship betweeen $$log(Error)$$ at x=$$2/3$$ Vs $$log(\Delta x)$$', ...
    'Interpreter','LaTeX','FontSize',18)       % plot title at the top
xlabel('$$log(\Delta x)$$',...                 % Label for x-axis
        'Interpreter','LaTeX', ...
        'FontSize',18)
ylabel('$$log(Error)$$',...           % Label for y-axis
        'Interpreter','LaTeX', ...
        'FontSize',18)
annotation('textbox', [.2 .7 .1 .1], 'String', ['slope using polyfit ',num2str(power(1))]);





