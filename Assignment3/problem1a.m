% @Author: Athul Vijayan
% @Date:   2014-09-14 17:46:59
% @Last Modified by:   Athul Vijayan
% @Last Modified time: 2014-09-14 19:17:20

clear all
clc

numIter = 5;                            % number of iterations
syms x;
initTerms = 2;                          % Number of terms in the sine expansion included for first iteration
params = cell(numIter,1);               % unknowns
trialFns = cell(numIter, 1);            % trial functions for each iterations
phi = 1 + sin(pi*x/2);
phiHat = cell(numIter, 1);

for i=1:numIter
    % ********************** collocation method ******************.
    collocPts = transpose(linspace(0, 1, initTerms +i +1)); collocPts = collocPts(2:end-1);
    collocFvec = sin(pi*collocPts/2);
    for j=1: initTerms +i -1
        % trialFns{i}{j} = (x^j)*(x-1);
        trialFns{i}{j} = ((-1)^(j-1)*((pi*x)^(2*(j-1) + 1))) / (factorial(2*(j-1) +1)*(2^(2*(j-1)+1)));
    end

    for k=1:size(collocPts, 1)
        for j=1:size(collocPts, 1)
            Kvec(k,j) = subs(trialFns{i}{j}, x, collocPts(k));
        end
    end
    collocParams{i} = inv(Kvec)*collocFvec;
    phiHat{i} = 1 + trialFns{i}*collocParams{i};

    figure
    phiPlot = ezplot(phi, [0, 1]);
    set(phiPlot,'LineWidth',2); set(phiPlot,'color','g'); 
    hold on;
    phiHatPlot = ezplot(phiHat{i}, [0, 1]);
    set(phiHatPlot,'LineWidth',2); set(phiHatPlot,'color','r'); set(phiHatPlot,'LineStyle','--');
    set(get(gca,'XLabel'),'String','x');
    set(get(gca,'YLabel'),'String','$$\phi, \quad \hat{\phi}$$', 'Interpreter','LaTex');
    L = legend('$$\phi = 1 + sin({\pi x \over 2})$$', '$$\hat{\phi}$$', 'Location','southeast');
    set(L,'Interpreter','LaTex');
    title(['Plot showing actual and estimated function using point collocation for iteration ', num2str(i)]);
    hold off;
end