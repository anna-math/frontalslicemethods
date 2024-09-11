%% Varying system type: Uniform and Gaussian
%% Script for Figure 7 (left)

close all
clc
clear


addpath('utils/')



rng(3)

%% Parameters
n1 = 100;
n2 = 20;
n3 = 10;
nn = 10;

paras.maxiter = 1000; 
paras.alpha = .00005;

saveBoo = false;

%% Plots
figure
set(0,'defaultAxesFontSize', 16)
set(0,'defaultlinelinewidth', 2)

%% Problem Set Up
XX = randn(n2,n3,nn);
XX = XX/tnorm(XX);


%% Recovery
AA = randn(n1,n2,nn);
BB = tprod(AA,XX);
[err_cyclic, ~, ~] = frontalDescentAll(AA, XX, BB, 'randomSampling', paras);
semilogy(err_cyclic, 'DisplayName', 'Gaussian')
hold on


[kap, mu] = evalKappaMu(AA,paras.alpha)


AA = .5*unifrnd(0,1, [n1,n2,nn]) + .5*randn(n1,n2,nn);
BB = tprod(AA,XX);
[err_cyclic, ~, ~] = frontalDescentAll(AA, XX, BB, 'randomSampling', paras);
semilogy(err_cyclic, 'DisplayName', 'Mixture')

[kap, mu] = evalKappaMu(AA,paras.alpha)

AA = unifrnd(0,1, [n1,n2,nn]);
BB = tprod(AA,XX);
[err_cyclic, ~, ~] = frontalDescentAll(AA, XX, BB, 'randomSampling', paras);
semilogy(err_cyclic, 'DisplayName', 'Uniform')

[kap, mu] = evalKappaMu(AA,paras.alpha)



legend('show', 'location', 'southwest')
xlabel('Iterations')
ylabel('Approximation Error')
if(saveBoo)
	fname = sprintf('%s_%dn1_%dn2_%dn3_%dnn_errorplot', mfilename(pwd), n1, n2, n3, nn);
	saveas(gcf, strcat(fname ,'.png'))
	saveFigure(strcat(fname ,'.fig'))
end



