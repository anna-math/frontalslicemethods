function [err, timeVec, Xt] = frontalDescentAll(AA, XX, BB, variation, paras)
% Frontal Descent - solve tensor linear system AX = B using frontal slices of A
%	** Implementation only works for 3rd order **

% Inputs
%	AA: 3rd order measurement tensor n1 x n2 x nn
%	XX: 3rd order solution tensor n2 x n3 x nn
%	BB: 3rd order resulting tensor n1 x n3 x nn
% 	variations: 'cyclic', 'leverageSampling', 'randomSampling', 'accumulative', 'full'
%	paras:
%		.maxiter: maximum number of iterations
%		.alpha: learning rate
%		.stepsize: step size for all variations (default 1)
%		.blocksize: number of frontal slices to use at each iteration, default 1
%		.bandsize: effective number of frontal slices, default nn
%		.tolerance: tolerance for early stopping (default 0)
%		
% Outputs
%	err: (paras.maxiter)x1 vector containing approximation error at each iteration
%	timeVec: (paras.maxiter)x1 vector containing elapsed cputime time at each iteration
%	Xt: approximation of solution to system

	%% Standard variable collection and pre-processing
	[n1,n2,nn] = size(AA);
	[~,n3,~] = size(XX);

	maxiter = paras.maxiter;
	alpha = paras.alpha;

	if isfield(paras,'bandsize')
	    bandsize = paras.bandsize;
	else
	    bandsize = nn;
	end

	if isfield(paras,'blocksize')
	    blocksize = paras.blocksize;
	else
	    blocksize = 1;
	end

	if isfield(paras,'stepsize')
	    stepsize = paras.stepsize;
	else
	    stepsize = 1;
	end

	if isfield(paras, 'tolerance')
	    tolerance = paras.tolerance;
	else
	    tolerance = 0;
	end

	numberOfFrontalBlocks = int32(bandsize/blocksize);
    if isfield(paras, 'initX')
        Xt = paras.initX;
    else
        Xt = zeros(n2,n3,nn);
    end
    Xhist = cell(1,nn);
    for jj = 1:numberOfFrontalBlocks
        Xhist{jj} = Xt;
    end
	Rt = BB;

    err = zeros(maxiter,1);
    timeVec = zeros(maxiter,1);

    err(1) = tnorm(XX);

    %% Additional preprocessing for variations
	switch variation
    	case 'cyclic'
        	    ii_vec = repmat(1:numberOfFrontalBlocks, [1, int32(maxiter/numberOfFrontalBlocks)+1]);
	    case 'randomSampling'
    		    ii_vec = randsample(numberOfFrontalBlocks, maxiter, true);
    	case 'leverageSampling'  
		    % Compute the norms of the frontal slices
		    slice_norms = zeros(1, numberOfFrontalBlocks);
		    for jj = 1:numberOfFrontalBlocks
		    	indx=max(1,(jj-1)*blocksize+1):min(jj*blocksize,nn);
		        slice_norms(jj) = tnorm(AA(:,:,indx));
		    end
		    % Normalize the slice norms to get sampling probabilities
		    sampling_probs = slice_norms / sum(slice_norms);
		    ii_vec = randsample(numberOfFrontalBlocks, maxiter, true, sampling_probs);
        case 'full'
            % Full sampling does not need additional preprocessing
            ii_vec = 1:maxiter; % Placeholder to allow the loop to run
	end

	%% Main loop
    for tt = 1:maxiter
        pointer = ii_vec(tt);
    	ss = cputime;

    	% Sample based on variation
        if strcmp(variation, 'full')
            Ai = AA;
        else
            pointer = pointer + stepsize;
            ii = mod(pointer, numberOfFrontalBlocks) + 1;
            indx = max(1, (ii-1) * blocksize + 1) : min(ii * blocksize, nn);
            Ai = zeros(n1, n2, nn);
            Ai(:,:,indx) = AA(:,:,indx);
        end

        % Residual update
        Rt = Rt - tprod(Ai,Xt) + tprod(Ai,Xhist{ii});
        Xhist{ii} = Xt;

        % Iterate update
        if strcmp(variation, 'accumulative')
            % Accumulative update
            Xt = Xt + alpha * tprod(tran(Ai), Rt) / (blocksize * maxiter);
        else
            Xt = Xt + alpha * tprod(tran(Ai), Rt);
        end
        
        % Some metrics
        timeVec(tt+1) = cputime-ss;  % iteration elapsed time
        err(tt+1) = tnorm(Xt - XX);
        
        % Early stopping based on tolerance
        if err(tt+1) < tolerance
            err = err(1:tt+1);
            timeVec = timeVec(1:tt+1);
            break;
        end

    end
end
