%% Varying block size experiment
% Script for Figure 8


close all
clc
clear


addpath('utils/')

rng(4)

%% Parameters
n1 = 100;
n2 = 20;
n3 = 10;
nn = 10;

paras.maxiter = 500; 
paras.alpha = .0001;

saveBoo = false;

%% Data Set Up
XX = randn(n2,n3,nn);
XX = XX / tnorm(XX);
AA = randn(n1,n2,nn);


%% Blurred image
BB = tprod(AA,XX);

%% Recovery
[err_TRK, time_TRK, X_TRK] = TRK(AA, XX, BB, paras);
paras.blocksize = 1;
[err_CFGD_b1, time_CFGD_b1, X_CFGD_b1] = frontalDescentAll(AA, XX, BB, 'cyclic', paras);
paras.blocksize = 5;
[err_CFGD_b5, time_CFGD_b5, X_CFGD_b5] = frontalDescentAll(AA, XX, BB, 'cyclic', paras);
paras.blocksize = 10;
[err_CFGD_b10, time_CFGD_b10, X_CFGD_b10] = frontalDescentAll(AA, XX, BB, 'cyclic', paras);


%% Evaluate the difference for a chosen slice
nslice = 5;
original_slice = XX(:,:,nslice);
recovered_CFGD_b1_slice = X_CFGD_b1(:,:,nslice);
recovered_CFGD_b5_slice = X_CFGD_b5(:,:,nslice);
recovered_CFGD_b10_slice = X_CFGD_b10(:,:,nslice);

mse_CFGD_b1 = immse(recovered_CFGD_b1_slice, original_slice);
psnr_CFGD_b1 = psnr(recovered_CFGD_b1_slice, original_slice);
ssim_CFGD_b1 = ssim(recovered_CFGD_b1_slice, original_slice);

mse_CFGD_b5 = immse(recovered_CFGD_b5_slice, original_slice);
psnr_CFGD_b5 = psnr(recovered_CFGD_b5_slice, original_slice);
ssim_CFGD_b5 = ssim(recovered_CFGD_b5_slice, original_slice);

mse_CFGD_b10 = immse(recovered_CFGD_b10_slice, original_slice);
psnr_CFGD_b10 = psnr(recovered_CFGD_b10_slice, original_slice);
ssim_CFGD_b10 = ssim(recovered_CFGD_b10_slice, original_slice);

%% Display metrics
fprintf('CFGD b=1 Recovery Metrics:\n');
fprintf('MSE: %0.3e\n', mse_CFGD_b1);
fprintf('PSNR: %0.3e\n', psnr_CFGD_b1);
fprintf('SSIM: %0.3e\n', ssim_CFGD_b1);
fprintf('\n')

fprintf('CFGD b=5 Recovery Metrics:\n');
fprintf('MSE: %0.3e\n', mse_CFGD_b5);
fprintf('PSNR: %0.3e\n', psnr_CFGD_b5);
fprintf('SSIM: %0.3e\n', ssim_CFGD_b5);
fprintf('\n')

fprintf('CFGD b = 10 Recovery Metrics:\n');
fprintf('MSE: %0.3e\n', mse_CFGD_b10);
fprintf('PSNR: %0.3e\n', psnr_CFGD_b10);
fprintf('SSIM: %0.3e\n', ssim_CFGD_b10);


%% Plots
set(0,'defaultAxesFontSize', 16)
set(0,'defaultlinelinewidth', 2)

figure
semilogy(err_CFGD_b1, 'LineWidth', 3, 'DisplayName', 's=1');
hold on 
semilogy(err_CFGD_b5, 'LineWidth', 3, 'DisplayName', 's=5');
semilogy(err_CFGD_b10, 'LineWidth', 3, 'DisplayName', 's=10');
xlabel('Iterations')
ylabel('Approximation Error')
legend('show', 'location', 'southwest')
if(saveBoo)
	fname = sprintf('%s_errorplot', mfilename(pwd));
	saveas(gcf, strcat(fname ,'.png'))
	saveFigure(strcat(fname ,'.fig'))
end






