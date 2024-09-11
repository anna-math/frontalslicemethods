close all
clc
clear
%% Script for Figure 7 (right)


addpath('utils/')

%% Parameters
n1 = 100;
n2 = 20;
n3 = 10;
nn = 10;
paras.maxiter = 1000; 
paras.alpha = .0001;
numSystems = 10;

saveBoo = false;

%% Arrays to store results
kappaVals = zeros(numSystems * 3, 1);
muVals = zeros(numSystems * 3, 1);
finalErrors = zeros(numSystems * 3, 1);
systemTypes = {'Gaussian', 'Mixture', 'Uniform'};
markers = {'o', 's', '^'};

%% Problem Set Up
XX = randn(n2,n3,nn);
XX = XX/tnorm(XX);

count = 1;

for st = 1:length(systemTypes)
    systemType = systemTypes{st};
    
    for i = 1:numSystems
        switch systemType
            case 'Gaussian'
                AA = randn(n1, n2, nn);
            case 'Mixture'
                AA = 0.5 * unifrnd(0, 1, [n1, n2, nn]) + 0.5 * randn(n1, n2, nn);
            case 'Uniform'
                AA = unifrnd(0, 1, [n1, n2, nn]);
        end
        
        BB = tprod(AA, XX);
        [err, ~, ~] = frontalDescentAll(AA, XX, BB, 'cyclic', paras);
        
        finalErrors(count) = err(end);
        [kappa, mu] = evalKappaMu(AA, paras.alpha);
        kappaVals(count) = kappa;
        muVals(count) = mu;
        
        fprintf('System: %s, Kappa: %.4f, Mu: %.4f, Final Error: %.4e\n', systemType, kappa, mu, err(end));
        count = count + 1;
    end
end

%% Scatter Plot
figure
hold on
for st = 1:length(systemTypes)
    systemType = systemTypes{st};
    marker = markers{st};
    indices = (st-1)*numSystems+1 : st*numSystems;
    scatter(kappaVals(indices), muVals(indices), 100, finalErrors(indices), marker, 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', systemType)
end
hold off
colorbar
xlabel('\kappa')
ylabel('\mu')
title('Scatter plot of systems with different \kappa and \mu')
legend('Location', 'best')
set(gca, 'FontSize', 14)

if saveBoo
    fname = sprintf('%s_%dn1_%dn2_%dn3_%dnn_scatterplot', mfilename(pwd), n1, n2, n3, nn);
    saveas(gcf, strcat(fname, '.png'))
    saveFigure(strcat(fname, '.fig'))
end