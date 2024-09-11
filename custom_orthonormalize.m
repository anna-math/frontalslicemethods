function Amat_modified = custom_orthonormalize(Amat,threshold)
    % Step 1: Normalize columns of Amat to have unit length
    for i = 1:size(Amat, 2)
        Amat(:, i) = Amat(:, i) / norm(Amat(:, i));
    end
    
    % Step 2: Iteratively adjust columns to ensure inner products are bounded
    max_inner_product = threshold;
    tolerance = 1e-10; % To prevent infinite loops due to numerical precision
    max_iterations = 1000; % Prevent infinite loops with a maximum iteration count
    
    for iter = 1:max_iterations
        adjusted = false; % Track if any adjustment is made in this iteration
        
        for i = 1:size(Amat, 2)
            for j = 1:i-1
                inner_product = dot(Amat(:, i), Amat(:, j));
                
                % Adjust if the inner product exceeds the desired threshold
                if abs(inner_product) > max_inner_product
                    adjustment = sign(inner_product) * (abs(inner_product) - max_inner_product) * Amat(:, j);
                    Amat(:, i) = Amat(:, i) - adjustment;
                    Amat(:, i) = Amat(:, i) / norm(Amat(:, i)); % Re-normalize
                    
                    adjusted = true;
                end
            end
        end
        
        % Break loop if no adjustments were made (all inner products are within bounds)
        if ~adjusted
            break;
        end
        
        % If we reached the max iterations, break (prevent infinite loop)
        if iter == max_iterations
            warning('Max iterations reached. The result might not satisfy the inner product constraint fully.');
            break;
        end
    end
    
    Amat_modified = Amat;
end
