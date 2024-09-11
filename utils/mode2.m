function B = mode2(A,U)

[m,p,n]=size(A); [r,q]=size(U);  

A2 = reshape(permute(A,[2,1,3]),[p,m*n]);

if q~=p
    disp('warning, inappropriate sizes.  Quitting')
    return
end

B = U*A2;

B = reshape(B,[r,m,n]); B = permute(B,[2,1,3]); 