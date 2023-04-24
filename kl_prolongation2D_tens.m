function LL=kl_prolongation2D_tens(nkl,C2x,C2y,p1,p2,l1,l2,IJ)
% IJ: mode pairs for eigenvalue tensor product
% nkl: number of eigenvalues 

tol=1e-7;
flag_analytic=1; %if 1 the linear systems are solved the analytic solution...if 0 matlab solves them

%% Expansion 1D (2 directions)
%if tensoritation the 1D expansions can be extended
[~,kx,~]=kl_prolongation(l1(1:max(IJ(:,1))),p1(:,1:max(IJ(:,1))),[],toeplitz(C2x),1);
[~,ky,~]=kl_prolongation(l2(1:max(IJ(:,2))),p2(:,1:max(IJ(:,2))),[],toeplitz(C2y),1);

%% Indices to pass in 2D
i1=repmat(IJ(:,1),[1,nkl]);
j1=repmat(IJ(:,1)',[nkl,1]);
i2=repmat(IJ(:,2),[1,nkl]);
j2=repmat(IJ(:,2)',[nkl,1]);
%% 2D Expansions
%x-->
Kx=kx(IJ(:,1),IJ(:,1)) .* (i2==j2);
%y-->
Ky=ky(IJ(:,2),IJ(:,2)) .* (i1==j1);
%x--> y-->
Ku=kx(IJ(:,1),IJ(:,1)) .* ky(IJ(:,2),IJ(:,2)) ;
%x--> y<--
Kd=kx(IJ(:,1),IJ(:,1)) .* ky(IJ(:,2),IJ(:,2))' ;

%% Decomposition
%%% voir notes
I=eye(nkl);O=zeros(nkl);

%% 1 x back-forward;
[L,flag]=chol(I-Kx*Kx'-Kx'*Kx,'lower');
if flag
    warning(['KK non positive: min eigenvalue= ',num2str(min(eig(I-Kx*Kx'-Kx'*Kx)))])
    [u,s,~]=svd(I-Kx*Kx'-Kx'*Kx);
    tmp=u*s*u'; %svd symm
    L=chol(tmp,'lower');
end
LL{1}=[Kx',Kx,L];

if flag_analytic
    %inversion of (I-Kx*Kx'-Kx'*Kx)=L*L'....for after...
    tmp=inv(L');
    Z=tmp*tmp';
end



% % [L,flag]=chol(I-Ky*Ky'-Ky'*Ky,'lower');
% % if flag
% %     warning(['KK non positive: min eigenvalue= ',num2str(min(eig(I-Kx*Kx'-Kx'*Kx)))])
% %     [u,s,~]=svd(I-Kx*Kx'-Kx'*Kx);
% %     tmp=u*s*u'; %svd symm
% %     L=chol(tmp,'lower');
% % end
% % LL{6}=[Ky',Ky,L];


%% 2 y back-forw with obliq.
%complete
if ~flag_analytic
    A=[I  ,  Kx , O; ...
        Kx',  I ,  Kx; ...
        O  ,  Kx', I   ];
    B1=[Ku',Ky',Kd] ;%system 1 (y-back)
    B2=[Kd',Ky,Ku] ;%system2 (y-up)
    
    T=corr_transform(A,[B1;B2],2);
    T(abs(T)<tol)=0;
    LL{2}=T;
    
    %1st x band
    T=corr_transform(A(nkl+1:end,nkl+1:end),[B1(:,nkl+1:end);B2(:,nkl+1:end)],2);
    T(abs(T)<tol)=0;
    LL{3}=T;
    
    %Last x band
    T=corr_transform(A(1:end-nkl,1:end-nkl),[B1(:,1:end-nkl);B2(:,1:end-nkl)],2);
    T(abs(T)<tol)=0;
    LL{4}=T;
else  %analytics...regarder les notes
    X2=(Ky'-Kd*Kx'-Ku'*Kx)*Z;
    X1=Ku'-X2*Kx';
    X3=Kd-X2*Kx;
    X5=(Ky-Ku*Kx'-Kd'*Kx)*Z;
    X4=Kd'-X5*Kx';
    X6=Ku-X5*Kx;
    X=[X1,X2,X3,X4,X5,X6];
    B=[Ku',Ky',Kd,Kd',Ky,Ku];
    tmp=eye(size(B,1)) - X*B';
    tmp=tril(tmp,-1)+tril(tmp)'; %%%impose symmetricity...1e-10 problems...
    [L,flag_chol]=chol(tmp,'lower' );
    if flag_chol;L=chol_correct(tmp);end
    T=[X,L];T(abs(T)<tol)=0;
    LL{2}=T;
    
    %1st x band
    X2=(Ky'-Kd*Kx')/(I-Kx*Kx');
    X3=Kd-X2*Kx;
    X5=(Ky-Ku*Kx')/(I-Kx*Kx');
    X6=Ku-X5*Kx;
    X=[X2,X3,X5,X6];
    B=[Ky',Kd,Ky,Ku];
    tmp=eye(size(B,1)) - X*B';
    tmp=tril(tmp,-1)+tril(tmp)'; %%%impose symmetricity...1e-10 problems...
    [L,flag_chol]=chol(tmp,'lower' );
    if flag_chol;L=chol_correct(tmp);end
    T=[X,L];T(abs(T)<tol)=0;
    LL{3}=T;
    
    %Last x band
    X2=(Ky'-Ku'*Kx)/(I-Kx'*Kx);
    X1=Ku'-X2*Kx';
    X5=(Ky-Kd'*Kx)/(I-Kx'*Kx);
    X4=Kd'-X5*Kx';
    X=[X1,X2,X4,X5];
    B=[Ku',Ky',Kd',Ky];
    tmp=eye(size(B,1)) - X*B';
    tmp=tril(tmp,-1)+tril(tmp)'; %%%impose symmetricity...1e-10 problems...
    [L,flag_chol]=chol(tmp,'lower' );
    if flag_chol;L=chol_correct(tmp);end
    T=[X,L];T(abs(T)<tol)=0;
    LL{4}=T;
    
    
end

%% Completely surrunded

if ~flag_analytic
    
    A=[I,  Kx,    O,    O,  O,    O,  Ky,    O       ; ...
        Kx',I,     Kx,   O,  O,    O,  Kd',   Ku    ; ...
        O,  Kx',   I,    O,  O,    O,  O,     Ky      ; ...
        O,  O,     O,    I,  Kx,   O,  Ky',   O       ; ...
        O,  O,     O,    Kx',I,    Kx, Ku',   Kd    ; ...
        O,  O,     O,    O,  Kx',  I,   O,     Ky'     ; ...
        Ky',Kd,  O,    Ky,   Ku,    O,  I,     O       ; ...
        O,  Ku', Ky',  O,    Kd',   Ky, O,     I      ];
    
    B=[Ku',  Ky',  Kd,  Kd'  ,Ky  ,Ku  ,Kx'  ,Kx];
    
    T=corr_transform(A,B,0);
    T(abs(T)<tol)=0;
    LL{5}=T;
else %%%voir les notes
    A2=Ky'-Ku'*Kx-Kd*Kx';
    A5=Ky -Kd'*Kx-Ku*Kx';
    A3=Ky'*Kx -Kd;
    A4=Ky'*Kx'  -Ku';
    A6=Ky*Kx-Ku;
    A7=Ky*Kx'-Kd';
    
    B1=Kx-Kd*Ky -Ku*Ky' + A2*Z*(Kx*Ky-Ku) + A5*Z*(Kx*Ky'-Kd);
    B2=A3*Z*(Kx*Ky-Ku) + A6*Z*(Kx*Ky'-Kd);
    
    W1=inv_chol(I-A4*Z*(Kx*Ky-Ku)-A7*Z*(Kx*Ky'-Kd)-Ky*Ky'-Ky'*Ky);
    W2=inv_chol(I-A3*Z*(Kx'*Ky-Kd')-A6*Z*(Kx'*Ky'-Ku')-Ky*Ky'-Ky'*Ky  - B2*W1*(A4*Z*(Kx'*Ky-Kd') + A7*Z*(Kx'*Ky'-Ku')) );
    
    X7=(Kx'-Ku'*Ky-Kd'*Ky' + A2*Z*(Kx'*Ky-Kd') + A5*Z*(Kx'*Ky'-Ku') +B1*W1*(A4*Z*(Kx'*Ky-Kd')+A7*Z*(Kx'*Ky'-Ku'))  ) *W2;
        
    X8=(B1+X7*B2)*W1;
    
    X2=(A2+ X7*A3+ X8*A4)*Z;
    X5=(A5+ X7*A6+ X8*A7)*Z;
    
    X3=Kd -X8*Ky'-X2*Kx;
    X1=Ku'-X7*Ky'-X2*Kx';
    X6=Ku -X8*Ky -X5*Kx;
    X4=Kd'-X7*Ky -X5*Kx';
    
    
    X=[X1,X2,X3,X4,X5,X6,X7,X8];
    B=[Ku',  Ky',  Kd,  Kd'  ,Ky  ,Ku  ,Kx'  ,Kx];
    tmp=eye(size(B,1)) - X*B';
    tmp=tril(tmp,-1)+tril(tmp)'; %%%impose symmetricity...1e-10 problems...
    [L,flag_chol]=chol(tmp,'lower' );
    if flag_chol;L=chol_correct(tmp);end
    T=[X,L];T(abs(T)<tol)=0;
    LL{5}=T;
    
    
end





return






