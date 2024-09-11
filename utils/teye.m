function tensorI = teye(pp,nn)

tensorI = zeros(pp,pp,nn);

tensorI(:,:,1) = eye(pp);

end