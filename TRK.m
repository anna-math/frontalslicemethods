function [err, cpu_time, Xt] = TRK(AA, Xopt, YY, paras)
% Tensor Randomized Kaczmarz - solve tensor linear system Ax = y for x
%	** Implementation only works for 3rd order **

% Inputs
%	AA: 3rd order measurement tensor
%	YY: 3rd order resulting tensor
%	Xopt: desired solution
%	paras:
%		.maxiter: maximum number of iterations
%		
% Outputs
%	Xt: approximation of solution to system
%	error: (paras.maxiter)x1 vector containing approximation error at each iteration
%	cpu_time: (paras.maxiter)x1 vector containing elapsed cpu time at each iteration
%	fps: flop counts


%% Parameters
[mm,ll,nn] = size(AA);
[~,kk,~] = size(YY);

%% Initialization
Xt = zeros(ll,kk,nn);
XX = Xopt;
err = [];
cpu_time = [];
fps = 0;

err(1) = tnorm(Xopt);
cpu_time(1) = 0;
maxiter = paras.maxiter;
err = zeros(maxiter,1);
timeVec = zeros(maxiter,1);
err(1) = tnorm(XX);
%% Main loop
for tt=1:maxiter
    
    ss = cputime

	%% Randomly select measurement
	ik = randsample(mm,1);

	%% Compute (AiAi^T)^{-1}
	Aik = AA(ik,:,:);
	Yik = YY(ik,:,:);
	VV = squeeze(tprod(Aik,tran(Aik))); 
	DDinv = zeros(1,1,nn);
	DDinv(1,1,:) = fft(1./(nn*fft(VV)));
	

	%% Update step
	Xt = Xt - tprod(tprod(tran(Aik), DDinv),tprod(Aik,Xt) - YY(ik,:,:));


	%% Compute metrics
    timeVec(tt+1) = cputime - ss;        
    err(tt+1) = tnorm(Xt - XX);
   
end
    
cpu_time = timeVec;

end