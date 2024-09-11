% Vary system types testing kappa mu allowed
% Script for Figure 4

close all
clc
clear

addpath('utils/')

%% Parameters
n1 = 100;
n2 = 10;
n3 = 20;
nn = 2;
num_batches = 1; % Number of batches
alpha_values = logspace(-5, 0, 12); % Define a range of alpha values in log scale

% Initialize storage for the results
kappa_mu_term_Gaussian = zeros(size(alpha_values));
kappa_mu_term_Mixture = zeros(size(alpha_values));
kappa_mu_term_Uniform = zeros(size(alpha_values));
kappa_mu_term_BoundedInner = zeros(size(alpha_values));
kappa_mu_term_BoundedInner2 = zeros(size(alpha_values));

%% Loop over different alpha values
for a = 1:length(alpha_values)
    alpha = alpha_values(a);
    
    %% Gaussian AA
    kappa_accum = 0;
    mu_accum = 0;
    for b = 1:num_batches
        AA = randn(n1,n2,nn); % Gaussian
        [kap, mu] = evalKappaMu(AA, alpha);
        kappa_accum = kappa_accum + kap;
        mu_accum = mu_accum + mu;
    end
    avg_kappa_Gaussian = kappa_accum / num_batches;
    avg_mu_Gaussian = mu_accum / num_batches;
    kappa_mu_term_Gaussian(a) = avg_kappa_Gaussian + alpha * avg_mu_Gaussian * (nn - 1);

    %% Diagonal AA
    kappa_accum = 0;
    mu_accum = 0;
    for b = 1:num_batches
        AA = zeros(n1,n2,nn); 
        for ii=1:min(n1,n2)
            for kk=1:nn
            AA(ii,ii,kk) = randn(1);
            end
        end
        [kap, mu] = evalKappaMu(AA, alpha);
        kappa_accum = kappa_accum + kap;
        mu_accum = mu_accum + mu;
    end
    avg_kappa_diag = kappa_accum / num_batches;
    avg_mu_diag = mu_accum / num_batches;
    kappa_mu_term_diag(a) = avg_kappa_diag + alpha * avg_mu_diag * (nn - 1);

     %% Mutually Orthogonal AA
    kappa_accum = 0;
    mu_accum = 0;
   
    for b = 1:num_batches
        AA = zeros(n1,n2,nn); 
        Amat = randn(n1, n2*nn);
 
        Amat = orth(Amat); size(Amat)
        for ii=1:nn
            AA(:,:,ii) = Amat(:,((ii-1)*n2+1):((ii)*n2));
        end
        [kap, mu] = evalKappaMu(AA, alpha);
        kappa_accum = kappa_accum + kap;
        mu_accum = mu_accum + mu;

        % Display the pairwise tprod for the slices in AA
        disp(['Batch ', num2str(b), ': Pairwise tprod for AA slices']);
        for i = 1:nn
            for j = 1:nn
                tprod_val = tnorm(tprod(tran(AA(:,:,i)), AA(:,:,j)));
                fprintf('tprod(tran(AA(:,:, %d)), AA(:,:, %d)) = %f\n', i, j, tprod_val);
            end
        end
    end
    avg_kappa_mutOrth = kappa_accum / num_batches;
    avg_mu_mutOrth = mu_accum / num_batches;
    kappa_mu_term_mutOrth(a) = avg_kappa_mutOrth + alpha * avg_mu_mutOrth * (nn - 1);

    %% Noisy Mutually Orthogonal AA
    kappa_accum = 0;
    mu_accum = 0;
   
    for b = 1:num_batches
        AA = zeros(n1,n2,nn); 
        Amat = randn(n1,n2*nn);
        Amat = orth(Amat);
        for ii=1:nn
            AA(:,:,ii) = Amat(:,((ii-1)*n2+1):((ii)*n2));
        end
        AA = AA + (1/nn^3)*randn(n1,n2,nn);
        [kap, mu] = evalKappaMu(AA, alpha);
        kappa_accum = kappa_accum + kap;
        mu_accum = mu_accum + mu;
    end
    avg_kappa_mutOrthNoisy = kappa_accum / num_batches;
    avg_mu_mutOrthNoisy = mu_accum / num_batches;
    kappa_mu_term_mutOrthNoisy(a) = avg_kappa_mutOrthNoisy + alpha * avg_mu_mutOrthNoisy * (nn - 1);

   %% Mutually inner product bounded AA
    kappa_accum = 0;
    mu_accum = 0;
   
    for b = 1:num_batches
        AA = zeros(n1,n2,nn); 
        Amat = randn(n1, n2*nn);
 
        Amat = custom_orthonormalize(Amat,0.1); size(Amat)
        for ii=1:nn
            AA(:,:,ii) = Amat(:,((ii-1)*n2+1):((ii)*n2));
        end
        [kap, mu] = evalKappaMu(AA, alpha);
        kappa_accum = kappa_accum + kap;
        mu_accum = mu_accum + mu;

        % Display the pairwise tprod for the slices in AA
        disp(['Batch ', num2str(b), ': Pairwise tprod for AA slices']);
        for i = 1:nn
            for j = 1:nn
                tprod_val = tnorm(tprod(tran(AA(:,:,i)), AA(:,:,j)));
                fprintf('tprod(tran(AA(:,:, %d)), AA(:,:, %d)) = %f\n', i, j, tprod_val);
            end
        end
    end
    avg_kappa_BoundedInner = kappa_accum / num_batches;
    avg_mu_BoundedInner = mu_accum / num_batches;
    kappa_mu_term_BoundedInner(a) = avg_kappa_BoundedInner + alpha * avg_mu_BoundedInner * (nn - 1);

    %% Mutually inner product bounded AA 2
    kappa_accum = 0;
    mu_accum = 0;
   
    for b = 1:num_batches
        AA = zeros(n1,n2,nn); 
        Amat = randn(n1, n2*nn);
 
        Amat = custom_orthonormalize(Amat,0.01); size(Amat)
        for ii=1:nn
            AA(:,:,ii) = Amat(:,((ii-1)*n2+1):((ii)*n2));
        end
        [kap, mu] = evalKappaMu(AA, alpha);
        kappa_accum = kappa_accum + kap;
        mu_accum = mu_accum + mu;

        % Display the pairwise tprod for the slices in AA
        disp(['Batch ', num2str(b), ': Pairwise tprod for AA slices']);
        for i = 1:nn
            for j = 1:nn
                tprod_val = tnorm(tprod(tran(AA(:,:,i)), AA(:,:,j)));
                fprintf('tprod(tran(AA(:,:, %d)), AA(:,:, %d)) = %f\n', i, j, tprod_val);
            end
        end
    end
    avg_kappa_BoundedInner = kappa_accum / num_batches;
    avg_mu_BoundedInner = mu_accum / num_batches;
    kappa_mu_term_BoundedInner2(a) = avg_kappa_BoundedInner + alpha * avg_mu_BoundedInner * (nn - 1);

    %% Mixture AA
    kappa_accum = 0;
    mu_accum = 0;
    for b = 1:num_batches
        AA = 0.5*unifrnd(0,1, [n1,n2,nn]) + 0.5*randn(n1,n2,nn); % Mixture
        [kap, mu] = evalKappaMu(AA, alpha);
        kappa_accum = kappa_accum + kap;
        mu_accum = mu_accum + mu;
    end
    avg_kappa_Mixture = kappa_accum / num_batches;
    avg_mu_Mixture = mu_accum / num_batches;
    kappa_mu_term_Mixture(a) = avg_kappa_Mixture + alpha * avg_mu_Mixture * (nn - 1);

    %% Uniform AA
    kappa_accum = 0;
    mu_accum = 0;
    for b = 1:num_batches
        AA = unifrnd(0,1, [n1,n2,nn]); % Uniform
        [kap, mu] = evalKappaMu(AA, alpha);
        kappa_accum = kappa_accum + kap;
        mu_accum = mu_accum + mu;
    end
    avg_kappa_Uniform = kappa_accum / num_batches;
    avg_mu_Uniform = mu_accum / num_batches;
    kappa_mu_term_Uniform(a) = avg_kappa_Uniform + alpha * avg_mu_Uniform * (nn - 1);
end

%% Plotting the results with log-log scale
figure;
loglog(alpha_values, kappa_mu_term_Gaussian, '-o', 'DisplayName', 'Gaussian', 'LineWidth', 2);
hold on;
loglog(alpha_values, kappa_mu_term_Mixture, '-s', 'DisplayName', 'Mixture', 'LineWidth', 2);
loglog(alpha_values, kappa_mu_term_Uniform, '-^', 'DisplayName', 'Uniform', 'LineWidth', 2);
loglog(alpha_values, kappa_mu_term_diag, '.-', 'DisplayName', 'Diagonal', 'LineWidth', 2);
loglog(alpha_values, kappa_mu_term_mutOrth, '--', 'DisplayName', 'Mutually Orthogonal', 'LineWidth', 2);
loglog(alpha_values, kappa_mu_term_mutOrthNoisy, '-x', 'DisplayName', 'Mutually Orthogonal (Noisy)', 'LineWidth', 2);
loglog(alpha_values, kappa_mu_term_BoundedInner, '-x', 'DisplayName', 'Mutual inner product is bounded by (0.1)', 'LineWidth', 2);
loglog(alpha_values, kappa_mu_term_BoundedInner2, '-x', 'DisplayName', 'Mutual inner product is bounded by (0.01)', 'LineWidth', 2);
yline(1, '--r', 'DisplayName', 'threshold 1', 'LineWidth', 2);
% Set axis labels with larger font size
xlabel('Alpha', 'FontSize', 20);
ylabel('\kappa + \alpha \cdot \mu \cdot (nn - 1)', 'FontSize', 20);
% Increase font size of the axis ticks
set(gca, 'FontSize', 16);
%title('Comparison of \kappa + \alpha \cdot \mu \cdot (nn - 1) for Different AA Types');
legend('Location', 'best', 'NumColumns', 8);
grid on;
hold off;
%% Helper function to generate matrices with desired eigenvalue property
function A = generate_matrix_with_min_eigenvalue(n1, n2, min_eigenvalue)
    % Generate a random matrix of size n1 x n2
    A = randn(n1, n2);
    
    % Ensure the eigenvalues of A'A are greater than the min_eigenvalue
    if n1 >= n2
        [V, D] = eig(A' * A);
    else
        [V, D] = eig(A * A');
    end
    
    % Adjust eigenvalues to ensure they are all greater than min_eigenvalue
    D(D < min_eigenvalue) = min_eigenvalue;
    
    % Reconstruct the matrix with adjusted eigenvalues
    if n1 >= n2
        A = A * V * sqrt(D) * V';
    else
        A = V * sqrt(D) * V' * A;
    end
end