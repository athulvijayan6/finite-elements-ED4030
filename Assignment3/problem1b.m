% @Author: Athul Vijayan
% @Date:   2014-09-13 12:47:10
% @Last Modified by:   Athul Vijayan
% @Last Modified time: 2014-09-14 19:21:43
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

    % ==============galerkin method ========================
    % We start with initTerms number ofterms in sine expansion and increase terms upon iteration.
    % construct trial functions matrix.
    trialFns{i} = cell(i + initTerms -1, 1);
    for j=0: initTerms +i -2
        trialFns{i}{j+1} = ((-1)^j*((pi*x)^(2*j + 1))) / (factorial(2*j +1)*(2^(2*j+1)));
    end
    % The trialFns has trial functions $$N_i$$ now.
    % RHS vector F
    FvecIntegrand = trialFns{i} * sin(pi*x/2);     %RHS integrand
    Fvec = double(int(FvecIntegrand, x, 0, 1));    %integrate
    clear('FvecIntegrand');

    % Construct K matrix
    for j=1:size(trialFns{i},1)
        for k=1:size(trialFns{i},1)
            KvecIntegrand(j,k) = trialFns{i}{j}*trialFns{i}{k}; % Kmatrix integrand
        end
    end
    Kvec = double(int(KvecIntegrand, x, 0, 1));     %integrate
    clear('KvecIntegrand');
    params{i} = inv(Kvec)*Fvec;                     %find unknowns
    phiHat{i} = 1;                                  %Find aproximate phi, phiHat
    collocPhiHat{i} = 1;
    for j=1:size(trialFns{i},1)
        phiHat{i} = phiHat{i} + params{i}(j)*trialFns{i}{j};
    end

    figure;
    phiPlot = ezplot(phi, [0, 1]);
    set(phiPlot,'LineWidth',2); set(phiPlot,'color','g');
    hold on;
    phiHatPlot = ezplot(phiHat{i}, [0, 1]);
    set(phiHatPlot,'LineWidth',2); set(phiHatPlot,'color','r'); set(phiHatPlot,'LineStyle','--');
    set(get(gca,'XLabel'),'String','x');
    set(get(gca,'YLabel'),'String','$$\phi, \quad \hat{\phi}$$', 'Interpreter','LaTex');
    L = legend('$$\phi = 1 + sin({\pi x \over 2})$$', '$$\hat{\phi}$$', 'Location','southeast');
    set(L,'Interpreter','LaTex');
    title(['Plot showing actual and estimated function using galerkin for iteration ', num2str(i)]);
    hold off
    clear('phiPlot', 'phiHatPlot', 'L');
end

