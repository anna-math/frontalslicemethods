function C=tprod(A,B)
% tprod compute the tensor-tensor product using FFTs
% A is n1 x n2 x n3; B is na x n4 x nb; then 
% output tensor is n1 x n4 x n3
% 
% Follows multiplication defined in Kilmer, Martin and Perrone, 2008
% 
% Copyright (C) Misha E. Kilmer, Tufts University,  misha.kilmer@tufts.edu

% modified by Sara Soltani, 
% August 2015, ssol@dtu.dk, sarahsoltani@gmail.com


[n1,n2,n3]=size(A);
[na,n4,nb]=size(B);

if n3~=nb || n2~=na
    disp('WARNING, dimensions are not acceptable')
    return
end

D=fft(A,[],3); 


Bhat=fft(B,[],3); 



C=zeros(n1,n4,n3);
for i=1:n3
    C(:,:,i)=D(:,:,i)*Bhat(:,:,i);
end

C=ifft(C,[],3);