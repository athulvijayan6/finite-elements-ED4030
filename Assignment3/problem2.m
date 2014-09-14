% @Author: Athul Vijayan
% @Date:   2014-09-14 11:27:24
% @Last Modified by:   Athul Vijayan
% @Last Modified time: 2014-09-14 17:31:03

for numTrialFns =4:8;                                        % do for iterations. each iteration adda a trial term
    syms x;   
    distance = transpose(linspace(0, 1, 6));                % given steps in x
    phi = [20 ; 30 ; 50 ; 65 ; 40 ; 30];                    %temperatures
    psi = 20 + 10*x;                                        % $$\psi$$ in estimate

    for i=1:numTrialFns
        trialFns(i) = (x^i)*(x-1);                          % find symbolic trial functions
    end

    for i=1:numTrialFns
        for j=1:numTrialFns
            KvecIntegrand(i,j) = trialFns(i)*trialFns(j); % construct K matrix
        end
        FvecIntegrand = double(times(subs(trialFns(i), x, distance) , phi - (subs(psi, x, distance))));
        Fvec(i) = trapz(distance, FvecIntegrand);           % do numerical integration 
    end
    Fvec = Fvec';
    Kvec = double(int(KvecIntegrand, x, 0, 1));             % integrate to get final K matrix

    params = inv(Kvec)*Fvec;                                % solved parameters
    phiHat = psi + trialFns*params;                         %Recreate symbolic estimate eqn

    % Now we plot both $$\phi$$ and $$\hat{\phi}$$.
    figure
    hold on
    phiPlot = plot(distance, phi);
    set(phiPlot,'LineStyle','o'); set(phiPlot,'color','b'); set(phiPlot,'MarkerFaceColor',[.49 1 .63]);
    phiHatPlot = ezplot(phiHat, [0, 1]);
    set(phiHatPlot,'LineWidth',2); set(phiHatPlot,'color','r');
    set(get(gca,'XLabel'),'String','x');
    set(get(gca,'YLabel'),'String','$$\phi, \quad \hat{\phi}$$', 'Interpreter','LaTex');
    L = legend('$$\phi$$', '$$\hat{\phi}$$', 'Location','southeast');
    set(L,'Interpreter','LaTex');
    title(['Plot showing actual and estimated function using galerkin for iteration ', num2str(numTrialFns - 3)]);
    clear all
    clc
end
