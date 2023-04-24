function X=inv_chol(A)  %%function inv of positive square matrix L*L'

[L,flag]=chol(A);
if flag
   [u,s,~]=svd(A);
   L=sqrt(s)*u';
end


tmp=inv(L);
X=tmp*tmp';
