% Figure 2 and tables


close all
clc 
clear 


addpath('utils/')
addpath('data/')

%% Parameters
band = 9;
sigma = 4;
new_size = [151, 151]; % Adjust as needed

paras.maxiter = 1000; 
paras.alpha = .01;
paras.stepsize = 1;
paras.bandsize = band+1;
paras.blocksize = 15;

saveBoo = false;

%% Data Set Up
load 'data/trimmedvideo.mat'
XX = reshape(Trimmed(25:175,75:225,1,100:130), [151, 151, 31]);
XX = double(XX);
[N,~] = size(XX);
XX = XX / tnorm(XX);

%% Image Resizing
XX_resized = zeros(new_size(1), new_size(2), size(XX, 3));
for i = 1:size(XX, 3)
    XX_resized(:, :, i) = imresize(XX(:, :, i), new_size);
end
XX_resized = permute(XX_resized, [1 3 2]);
[N1, N2, N3] = size(XX_resized);
XX_resized = XX_resized / tnorm(XX_resized);

%% Deblurring Operator
z1 = [exp(-((0 : band - 1).^2)/(2*sigma^2)), zeros(1, N1 - band)];
A2 = 1/(sigma*sqrt(2*pi)) * toeplitz(z1);
z2 = [z1(1) fliplr(z1(end - length(z1) + 2 : end))];
A1 = 1/(sigma*sqrt(2*pi)) * toeplitz(z1, z2);

AA = zeros(N1, N1, N1);
for ii = 1:N1
    AA(:, :, ii) = A1(ii, 1) * A2;
end

% Normalize row slices for TRK
for ii = 1:N1
    AA(ii, :, :) = AA(ii, :, :) / tnorm(AA(ii, :, :));
end

%% Measurements with resized data
BB_resized = tprod(AA, XX_resized);

%% Recovery with timing
tic;
[err_TRK, ~, X_TRK] = TRK(AA, XX_resized, BB_resized, paras);
time_TRK = toc;

paras.blocksize = 1;
tic;
[err_CFGD_b1, ~, X_CFGD_b1] = frontalDescentAll(AA, XX_resized, BB_resized, 'cyclic', paras);
time_CFGD_b1 = toc;

paras.blocksize = 5;
tic;
[err_CFGD_b5, ~, X_CFGD_b5] = frontalDescentAll(AA, XX_resized, BB_resized, 'cyclic', paras);
time_CFGD_b5 = toc;

paras.blocksize = 10;
tic;
[err_CFGD_b10, ~, X_CFGD_b10] = frontalDescentAll(AA, XX_resized, BB_resized, 'cyclic', paras);
time_CFGD_b10 = toc;

tic;
[err_RFGD, ~, X_RFGD] = frontalDescentAll(AA, XX_resized, BB_resized, 'randomSampling', paras);
time_RFGD = toc;

%% Evaluate the difference for a chosen slice
nslice = 15; % Choose a slice to display
original_slice_resized = imresize(squeeze(XX_resized(:, nslice, :)), new_size);
recovered_TRK_slice = permute(X_TRK(:, nslice, :), [1, 3, 2]);
recovered_CFGD_b1_slice = permute(X_CFGD_b1(:, nslice, :), [1, 3, 2]);
recovered_CFGD_b5_slice = permute(X_CFGD_b5(:, nslice, :), [1, 3, 2]);
recovered_CFGD_b10_slice = permute(X_CFGD_b10(:, nslice, :), [1, 3, 2]);
recovered_RFGD_slice = permute(X_RFGD(:, nslice, :), [1, 3, 2]);

% Ensure the sizes match
if ~isequal(size(original_slice_resized), size(recovered_TRK_slice)) || ...
   ~isequal(size(original_slice_resized), size(recovered_CFGD_b1_slice))
    error('The dimensions of the original and recovered slices do not match.');
end

% Compute error metrics
recovered_TRK_slice = real(recovered_TRK_slice);
mse_TRK = immse(recovered_TRK_slice, original_slice_resized);
psnr_TRK = psnr(recovered_TRK_slice, original_slice_resized);
ssim_TRK = ssim(recovered_TRK_slice, original_slice_resized);

mse_CFGD_b1 = immse(recovered_CFGD_b1_slice, original_slice_resized);
psnr_CFGD_b1 = psnr(recovered_CFGD_b1_slice, original_slice_resized);
ssim_CFGD_b1 = ssim(recovered_CFGD_b1_slice, original_slice_resized);

mse_CFGD_b5 = immse(recovered_CFGD_b5_slice, original_slice_resized);
psnr_CFGD_b5 = psnr(recovered_CFGD_b5_slice, original_slice_resized);
ssim_CFGD_b5 = ssim(recovered_CFGD_b5_slice, original_slice_resized);

mse_CFGD_b10 = immse(recovered_CFGD_b10_slice, original_slice_resized);
psnr_CFGD_b10 = psnr(recovered_CFGD_b10_slice, original_slice_resized);
ssim_CFGD_b10 = ssim(recovered_CFGD_b10_slice, original_slice_resized);

mse_RFGD = immse(recovered_RFGD_slice, original_slice_resized);
psnr_RFGD = psnr(recovered_RFGD_slice, original_slice_resized);
ssim_RFGD = ssim(recovered_RFGD_slice, original_slice_resized);

%% Display metrics and timing
fprintf('TRK Recovery Metrics:\n');
fprintf('MSE: %0.3e\n', mse_TRK);
fprintf('PSNR: %0.3e\n', psnr_TRK);
fprintf('SSIM: %0.3e\n', ssim_TRK);
fprintf('Time: %0.3f seconds\n', time_TRK);
fprintf('\n');

fprintf('CFGD b=1 Recovery Metrics:\n');
fprintf('MSE: %0.3e\n', mse_CFGD_b1);
fprintf('PSNR: %0.3e\n', psnr_CFGD_b1);
fprintf('SSIM: %0.3e\n', ssim_CFGD_b1);
fprintf('Time: %0.3f seconds\n', time_CFGD_b1);
fprintf('\n');

fprintf('CFGD b=5 Recovery Metrics:\n');
fprintf('MSE: %0.3e\n', mse_CFGD_b5);
fprintf('PSNR: %0.3e\n', psnr_CFGD_b5);
fprintf('SSIM: %0.3e\n', ssim_CFGD_b5);
fprintf('Time: %0.3f seconds\n', time_CFGD_b5);
fprintf('\n');

fprintf('CFGD b=10 Recovery Metrics:\n');
fprintf('MSE: %0.3e\n', mse_CFGD_b10);
fprintf('PSNR: %0.3e\n', psnr_CFGD_b10);
fprintf('SSIM: %0.3e\n', ssim_CFGD_b10);
fprintf('Time: %0.3f seconds\n', time_CFGD_b10);
fprintf('\n');

fprintf('RFGD Recovery Metrics:\n');
fprintf('MSE: %0.3e\n', mse_RFGD);
fprintf('PSNR: %0.3e\n', psnr_RFGD);
fprintf('SSIM: %0.3e\n', ssim_RFGD);
fprintf('Time: %0.3f seconds\n', time_RFGD);
fprintf('\n');


%% Plots
set(0,'defaultAxesFontSize', 16)
set(0,'defaultlinelinewidth', 2)

figure
semilogy(err_CFGD_b1, 'LineWidth', 3, 'DisplayName', 'Frontal, b=1');
hold on 
semilogy(err_CFGD_b5, 'LineWidth', 3, 'DisplayName', 'Frontal, b=5');
semilogy(err_CFGD_b10, 'LineWidth', 3, 'DisplayName', 'Frontal, b=10');
semilogy(err_RFGD, 'LineWidth', 3, 'DisplayName', 'Random');
semilogy(err_TRK, 'LineWidth', 3, 'DisplayName', 'Rows');

legend('show')
xlabel('Iterations')
ylabel('Approximation Error')
if(saveBoo)
    fname = sprintf('%s_errorplot', mfilename(pwd));
    saveas(gcf, strcat(fname ,'.png'))
    saveFigure(strcat(fname ,'.fig'))
end

figure
set(gcf, 'PaperPosition', [0 0 21 5])
subplot(1,4,1)
imagesc(reshape(XX_resized(:,nslice,:), [N1,N1]))
title('Original Image')
subplot(1,4,2)
imagesc(reshape(BB_resized(:,nslice,:), [N1,N1]))
title('Blurred Image')
subplot(1,4,3)
imagesc(reshape(X_CFGD_b1(:,nslice,:), [N1,N1]))
title('Frontal')
subplot(1,4,4)
imagesc(reshape(abs(X_TRK(:,nslice,:)), [N1,N1]))
title('Row')
if(saveBoo)
    fname = sprintf('%s_images', mfilename(pwd));
    saveas(gcf, strcat(fname ,'.png'))
    saveFigure(strcat(fname ,'.fig'))
end
