%% Vary number of frontal slices nn
%% Script for Figure 6


close all
clc
clear


addpath('utils/')

rng(2)

%% Parameters
n1 = 100;
n2 = 20;
n3 = 10;
nnVec = 2:1:50;
numTrials = 20;

paras.maxiter = 1000; 
paras.alpha = .001;

saveBoo = false;

methods = {'cyclic', 'randomSampling', 'leverageSampling'};
colors = {'b', 'g', 'r'};
lineStyles = {'-', '--', '-.'};

avgFinalError = zeros(length(nnVec), length(methods));

for tt = 1:numTrials
    tt
    for m = 1:length(methods)
        method = methods{m};
        finalError = zeros(length(nnVec),1);
        for ii = 1:length(nnVec)
            ii
            nn = nnVec(ii);

            %% Problem Set Up
            AA = randn(n1,n2,nn);
            XX = randn(n2,n3,nn);
            XX = XX/tnorm(XX);
            BB = tprod(AA,XX);

            %% Recovery
            [err, time, X] = frontalDescentAll(AA, XX, BB, method, paras);
            
            finalError(ii) = err(end);
        end
        avgFinalError(:, m) = avgFinalError(:, m) + finalError;
    end
end

avgFinalError = avgFinalError/numTrials;

%% Plots
figure
set(0,'defaultAxesFontSize', 16)
set(0,'defaultlinelinewidth', 2)
hold on
for m = 1:length(methods)
    semilogy(nnVec, avgFinalError(:, m), 'Color', colors{m}, 'LineStyle', lineStyles{m}, 'DisplayName', methods{m});
end
hold off
xlabel('Number of Frontal Slices')
ylabel('Approximation Error')
legend('show')

if(saveBoo)
    fname = sprintf('%s_%dn1_%dn2_%dn3_errorplot', mfilename(pwd), n1, n2, n3);
    saveas(gcf, strcat(fname ,'.png'))
    saveFigure(strcat(fname ,'.fig'))
end
