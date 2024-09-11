function [kappa, mu] = evalKappaMu(AA, alpha)

	[n1,n2,nn] = size(AA);

	%% Computing kappa and mu
	kappavec = zeros(1,nn);
	for ii=1:nn
		Ai = zeros(n1,n2,nn);
		Ai(:,:,ii) = AA(:,:,ii);
		kappavec(ii) = computeOperatorNorm(bcirc(teye(n2,nn) - alpha * tprod(tran(Ai), Ai)));
	end

	mumat = zeros(nn,nn);
	for ii=1:nn
		for jj=(ii+1):nn
			Ai = zeros(n1,n2,nn);
			Ai(:,:,ii) = AA(:,:,ii);
			Aj = zeros(n1,n2,nn);
			Aj(:,:,jj) = AA(:,:,jj);
			mumat(ii,jj) = computeOperatorNorm(bcirc(tprod(tran(Ai), Aj)));
		end
	end

	kappa = max(kappavec);
	mu =  max(max(triu(mumat, 1)));
end