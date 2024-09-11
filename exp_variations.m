%% Compare variations: cyclic, leverage score, random, TRK
%% Script for Figure 5

close all
clc
clear


addpath('utils/')

rng(4)

%% Parameters
n1 = 100;
n2 = 10;
n3 = 10;
nn = 20;

paras.maxiter = 1000; 
paras.alpha = 4;

saveBoo = false;

%% Data Set Up
XX = randn(n2,n3,nn);
XX = XX / tnorm(XX);
AA = randn(n1,n2,nn);

for ii=1:nn
	AA(:,:,ii) = AA(:,:,ii)/tnorm( AA(:,:,ii));
end


%% Blurred image
BB = tprod(AA,XX);

%% Recovery
[err_TRK, time_TRK, X_TRK] = TRK(AA, XX, BB, paras);
[err_FGD_cyclic, time_FGD_cyclic, X_FGD_cyclic] = frontalDescentAll(AA, XX, BB, 'cyclic', paras);
[err_FGD_random, time_FGD_random, X_FGD_random] = frontalDescentAll(AA, XX, BB, 'randomSampling', paras);
[err_FGD_leverage, time_FGD_leverage, X_FGD_leverage] = frontalDescentAll(AA, XX, BB, 'leverageSampling', paras);

%% Plots
set(0,'defaultAxesFontSize', 16)
set(0,'defaultlinelinewidth', 2)

figure
semilogy(err_FGD_cyclic, 'LineWidth', 3, 'DisplayName', 'Cyclic');
hold on 
semilogy(err_FGD_random, 'LineWidth', 3, 'DisplayName', 'Random Sampling');
hold on 
semilogy(err_FGD_leverage, 'LineWidth', 3, 'DisplayName', 'Leverage Score Sampling');
hold on 
semilogy(err_TRK, 'LineWidth', 3, 'DisplayName', 'TRK');
hold on 
xlabel('Iterations')
ylabel('Approximation Error')
if(saveBoo)
	fname = sprintf('%s_errorplot', mfilename(pwd));
	saveas(gcf, strcat(fname ,'.png'))
	saveFigure(strcat(fname ,'.fig'))
end






