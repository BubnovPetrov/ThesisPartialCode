function [LL,kk,err]=kl_prolongation(VAL,VEC,Cov,Cov2,flag_k)

if ~exist('flag_k','var');flag_k=0;end

if ~exist('Cov2','var')
    Cov2=extend_cov(Cov,2);
    Cov2=toeplitz(spectral_condition(Cov2(1,:)));
end

CC=Cov2(1:length(Cov2)/2,length(Cov2)/2+1:end); %target cross-corr function

kk=diag(1./sqrt(VAL))*VEC'*CC*VEC*diag(1./sqrt(VAL));%sa projection sur base base de KL

%%%matrices de passage....voir note pour le calcul
I=eye(numel(VAL));
%backward -forward
if flag_k; %%%just want kk
    err=0;
    LL=0;
else
    [LL,flag]=chol(I-kk'*kk-kk*kk','lower'); %from backward correlation
    if flag
        warning(['KK non positive: min eigenvalue= ',num2str(min(eig(I-kk'*kk-kk*kk')))])
        [u,s,~]=svd(I-kk'*kk-kk*kk');
        tmp=u*s*u'; %svd symm
        LL=chol(tmp,'lower');
    end
    
    % % RR=chol([I,zeros(length(I)),kk;zeros(length(I)),I,kk';kk',kk,I],'lower');
    % % R1=RR(end-numel(VAL)+1:end,1:numel(VAL));
    % % R2=RR(end-numel(VAL)+1:end,1+numel(VAL):2*numel(VAL));
    % % R3=RR(end-numel(VAL)+1:end,1+numel(VAL)*2:end);
    % % LL=R3;
    
    v1=VEC*VEC(1,:)';
    v2=VEC*VEC(end,:)';
    err=1 - v2'*CC*v1 / CC(end,1) ;
    
end

%evolution of error
phi_L0=VEC(end,:)'*VEC(1,:);
KK=kk.*(sqrt(VAL*VAL'));
AA=KK.*phi_L0;
for ii=1:length(AA);
    er(ii)=1-sum(sum(AA(1:ii,1:ii)));
end

%% evolution approx
E=exp(1i*(pi*(0:length(Cov)-1)'*(1:length(VAL))/length(Cov)) )  /sqrt(length(Cov))*sqrt(2);
E=E.*(ones(length(Cov),1)*sign(VEC(1,:)));

V=(diag(1./sqrt(VAL))*VEC')';


return
L=100;
nv=size(VEC,2);
ds=L/size(Cov,1);
t=0:ds:L-ds;
ff=linspace(1/L,1/ds/2,1000);
vcos=(cos(2*pi*t'*ff));
vcos=vcos./repmat(sqrt(sum(vcos.^2)),[size(vcos,1),1]);
v1c=vcos(:,1:nv)*vcos(1,1:nv)';
v2c=vcos(:,1:nv)*vcos(end,1:nv)';





end





