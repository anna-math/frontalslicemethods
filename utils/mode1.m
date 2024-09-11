function B = mode1(A,U)

[m,p,n]=size(A);  [r,q] = size(U); 

if q ~= m
    disp('warning, sizes do not match.  Quitting.')
    return
end

%multiply U into mode 1.
A1 = reshape(A,[m,p*n]);  %mode 1 unfolding

B = U*A1;  %left mat-mat product

B = reshape(B,[r,p,n]);