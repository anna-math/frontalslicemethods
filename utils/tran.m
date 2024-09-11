function Z=tran(A)
% return the tensor transpose of A, using definition
% in Kilmer, Martin and Perrone, 2008
%

% Copyright: Misha E. Kilmer, 2009
% Tufts University,  misha.kilmer@tufts.edu

% modified by Sara Soltani, 
% August 2015, ssol@dtu.dk, sarahsoltani@gmail.com


[~,~,n3]=size(A);
idx=[1,n3:-1:2];    

Z=A(:,:,idx);
Z=permute(Z,[2,1,3]);



