%% Deblurring experiment (real world video and images)
%% Script for Figure 9, 10


close all
clc
clear


addpath('utils/')
addpath('data/')


rng(0)

%% Parameters
band = 9;
sigma = 3;
dataset = 'image'; %image or video


paras.maxiter = 200; 
paras.alpha = .01;
paras.bandsize = band+1;

%paras.stepsize = 1;

saveBoo = false;

%% DATA SET UP
switch dataset
	case 'image'
		new_size = [60, 60];
		XX = imread('data/HSTgray.jpg');
		XX = double(XX);
		XX = imresize(XX, new_size);
		XX = permute(XX, [1 3 2]);
		[N,~] = size(XX);
		XX = XX / tnorm(XX);
		nslice = 1;
	case 'video'
		load 'data/trimmedvideo.mat'
		XX = reshape(Trimmed(25:175,75:225,1,100:130), [151, 151, 31]);
		XX = double(XX);
		XX = permute(XX, [1 3 2]);
		[N,~] = size(XX)
		XX = XX/tnorm(XX);
		nslice = 15;
end


%% Deblurring operator
z1 = [exp(-((0 : band - 1).^2)/(2*sigma^2)), zeros(1, N - band)];
A2 = 1/(sigma*sqrt(2*pi)) * toeplitz(z1);
z2 = [z1(1) fliplr(z1(end - length(z1) + 2 : end))];
A1 =  1/(sigma*sqrt(2*pi)) * toeplitz(z1,z2);

AA = zeros(N,N,N);
for ii=1:N
	AA(:,:,ii) = A1(ii,1) * A2;
end

% Normalize row slices for TRK
for ii=1:N
    AA(ii,:,:) =  AA(ii,:,:)/tnorm(AA(ii,:,:));
end

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
original_slice = reshape(XX(:,nslice,:), [N,N]);
recovered_TRK_slice = reshape(abs(X_TRK(:,nslice,:)), [N,N]);
recovered_CFGD_b1_slice = reshape(X_CFGD_b1(:,nslice,:), [N,N]);
recovered_CFGD_b5_slice = reshape(X_CFGD_b5(:,nslice,:), [N,N]);
recovered_CFGD_b10_slice = reshape(X_CFGD_b10(:,nslice,:), [N,N]);


% Compute error metrics
mse_TRK = immse(recovered_TRK_slice, original_slice);
psnr_TRK = psnr(recovered_TRK_slice, original_slice);
ssim_TRK = ssim(recovered_TRK_slice, original_slice);

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
fprintf('TRK Recovery Metrics:\n');
fprintf('MSE: %0.3e\n', mse_TRK);
fprintf('PSNR: %0.3e\n', psnr_TRK);
fprintf('SSIM: %0.3e\n', ssim_TRK);
fprintf('\n')

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
semilogy(err_CFGD_b1, 'LineWidth', 3, 'DisplayName', 'Cyclic, b=1');
hold on 
semilogy(err_CFGD_b5, 'LineWidth', 3, 'DisplayName', 'Cyclic, b=5');
semilogy(err_CFGD_b10, 'LineWidth', 3, 'DisplayName', 'Cyclic, b=10');
semilogy(err_TRK, 'LineWidth', 3, 'DisplayName', 'TRK');
legend('show')
xlabel('Iterations')
ylabel('Approximation Error')
if(saveBoo)
	fname = sprintf('%s_%sDataset_errorplot', mfilename(pwd), dataset);
	saveas(gcf, strcat(fname ,'.png'))
	saveFigure(strcat(fname ,'.fig'))
end

figure
set(gcf, 'PaperPosition', [0 0 17 5])
subplot(1,3,1)
imagesc(reshape(BB(:,nslice,:), [N,N]))
title('Blurred Image')
subplot(1,3,2)
imagesc(recovered_CFGD_b5_slice)
title('Cyclic, b=5')
subplot(1,3,3)
imagesc(recovered_TRK_slice)
title('TRK')

if(saveBoo)
	fname = sprintf('%s_%sDataset_images', mfilename(pwd), dataset);
	saveas(gcf, strcat(fname ,'.png'))
	saveFigure(strcat(fname ,'.fig'))
end




