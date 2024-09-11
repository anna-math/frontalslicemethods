function opNorm = computeOperatorNorm(A)
    % computeOperatorNorm computes the operator norm (spectral norm) of a matrix A
    % Input: 
    %   A - an m x n matrix
    % Output: 
    %   opNorm - the operator norm of matrix A

    % Compute the singular values of the matrix A
    singularValues = svd(A);
    
    % The operator norm is the largest singular value
    opNorm = max(singularValues);
end
