function B = mode3(A,U)

[m,p,n]=size(A);  [r,q]=size(U);

A3=reshape(permute(A,[3,1,2]),[n,p*m]); 

if n~=q
    disp('warning, wrong sizes. Quititng.'); 
    return
end

B = U*A3; B = reshape(B,[r,m,p]); B = permute(B,[2,3,1]);


