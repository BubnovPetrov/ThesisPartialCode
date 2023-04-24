function C=corr2D_tens(Cx,Cy)

if size(Cx,1)==size(Cx,2) %%% deja covariance sous forme de matrice
    CX=Cx;
else
    CX=toeplitz(Cx);
end
CY=toeplitz(Cy);

C=repmat(CX,[length(Cy),length(Cy)]);
for ii=1:length(Cy)
    for jj=1:length(Cy)
        indx=(ii-1)*length(Cx)+1:ii*length(Cx);
        indy=(jj-1)*length(Cx)+1:jj*length(Cx);
        C(indx,indy)=C(indx,indy)*CY(ii,jj);
    end
end

return

