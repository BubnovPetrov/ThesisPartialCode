clc;clear;close all;
set(0,'DefaultFigureColor',[1 1 1]);set(0,'DefaultTextInterpreter','latex');set(0,'defaultlinelinewidth',2);set(0,'DefaultAxesFontSize',20)
set(0,'DefaultFigureColormap',jet());close

rng(33)

%% Data
L=100; %length
ds=1;  %stepsize
lc=10; %correlation length
N=1000;   %number of generated sample
np=2;    %number of partitions
corr_model='exp'; %correlation model (see below)
err_kl=1e-3;    %where stop the KLE

%% For plotting
plot_on=1;
cmap=hot;cmap=flipud([cmap(:,1),cmap(:,2),cmap(:,3)]);close;
savefig_on=0;
%% Correlation
t=0:ds:L-ds;
t2=0:ds:2*L-ds;
switch corr_model
    case 'exp'
        corrfun=@(t,lc) exp(-t/lc);
    case 'gauss'
        corrfun=@(t,lc) exp(-(t/lc).^2);
    case 'sin'
        corrfun=@(t,lc) (t~=0).*sin(t/(lc/10))./(t/(lc/10)) ;
    case 'cos'
        corrfun=@(t,lc) exp(-t/lc).*cos(t/lc);
    case 'linear'
        corrfun=@(t,lc) ( abs(t)<=lc ).*(lc-abs(t))/lc;
    case 'bess'
        corrfun=@(t,lc) besselj ( 0, t/(lc/10) );
    case 'white'
        corrfun=@(t,lc) zeros(size(t))+1*(t==0);
    otherwise
       corrfun=@(t,lc) exp(-t/lc);
end
C=feval(corrfun,t,lc);  C(1)=1;
C2=feval(corrfun,t2,lc); C2(1)=1;
%% KL decomposition
[phi,lam]=eig(toeplitz(C)); [lam,ord]=sort(diag(lam),'descend');phi=phi(:,ord);
lam(lam<0)=0; %%just in case lam=-eps
nkl=find(1-cumsum(lam)/sum(lam)<err_kl,1);
lam=lam(1:nkl);phi=phi(:,1:nkl);

PSD=real(fft([C,fliplr(C(2:end))])); PSD(PSD<=0)=1e-300;


if plot_on
    figure(1);box on;hold on;grid on;grid minor
    plot(t,C,'b','linewidth',4)
    xlabel('Lag [-]');ylabel('Correlation function [-]')
    set(gca,'xtick',0:0.25*L/ds:L/ds);
    set(gca,'xticklabel',{'0','0.25L','0.5L','0.75L','L'})
    
    figure(2);box on;hold on;grid on;grid minor
    plot(linspace(0,1/ds/2,numel(C)),10*log(PSD(1:numel(C))),'r','linewidth',4)
    xlabel('Frequency [-]'),ylabel('Power spectral density [-], in dB')
    set(gca,'xtick',linspace(0,1/ds/2,5));xlim([0,1/ds/2])
    set(gca,'TickLabelInterpreter', 'latex');
    set(gca,'xticklabel',{'0','$\frac{n_s}{8L}$','$\frac{n_s}{4L}$','$\frac{3n_s}{8L}$','$\frac{n_s}{2L}$'})
    %ylim([max(-10,min(10*log10(PSD))),inf])

    figure(3);box on; hold on;grid on;grid minor
    plot(t,phi(:,1),'--','linewidth',1.5)
    plot(t,phi(:,2),'--','linewidth',1.5)
    plot(t,phi(:,3),'--','linewidth',1.5)
    plot(t,phi(:,4),'--','linewidth',1.5)
    xlabel('t [-]');ylabel('KL Mode [-]')
    
% % %     if savefig_on
% % %        figure(1)
% % %        savefig(['./fig2/',corr_model,'_corr'])
% % %        export_fig(gcf,['./fig2/',corr_model,'_corr'],'-pdf','-transparent')
% % %        figure(2)
% % %        savefig(['./fig2/',corr_model,'_psd'])
% % %        export_fig(gcf,['./fig2/',corr_model,'_psd'],'-pdf','-transparent')  
% % %     end
% % %     
    
end

%% Prolongate generation
[LL,kk,err]=kl_prolongation(lam,phi,toeplitz(C),toeplitz(C2)); %%matrice prolong

%initialisation
F=zeros(L/ds*np,N);
F_stick=zeros(L/ds*np,N);

%germs
etag=randn(nkl,N,np); %indipendent
etac=etag; %conditionated

%conditioning the germs

for ii=2:2:np
    if ii==np
        etac (:,:,ii)= kk'*etag(:,:,ii-1) + chol(eye(nkl)-kk'*kk,'lower')*etag(:,:,ii);
    else        
        etac (:,:,ii)= kk'*etag(:,:,ii-1) + kk*etag(:,:,ii+1) + LL*etag(:,:,ii); %apply transform forw-back...ii<np for backwarding
    end
end

for ii=1:np
    F((ii-1)*L/ds +1 :ii*L/ds,:)      =phi*sqrt(diag(lam))*etac(:,:,ii);
    F_stick((ii-1)*L/ds +1 :ii*L/ds,:)=phi*sqrt(diag(lam))*etag(:,:,ii);
end

K2=cov_est(L*np,ds,F',0); %covariance
K_stick=cov_est(L*np,ds,F_stick',0);
t2=0:ds:L*np-ds;

KK_err=zeros(length(C));
for ii=1:L/ds
   KK_err(ii,:)=(C-K2(ii,ii:L/ds+ii-1)); 
end

Kj=K2(L/ds,L/ds:L/ds*2-1); %covariance a la jonction

if plot_on
    
%     figure(1);box on;hold on;grid on;grid minor
%     plot(t,Kj,'r-.','linewidth',2)
%     legend('Target','Break','location','best')
    
    figure(5);box on;hold on;grid on;title('correlation after conditioning')
    imagesc(t2,t2,K2);colormap(cmap);colorbar;caxis([min(C),max(C)]);
    xlabel('t [-]');ylabel('t [-]');
    xlim([-inf inf]);ylim([-inf inf]);set(gca,'xtick',0:L*np/5:L*np);set(gca,'ytick',0:L*np/5:L*np)
    for ii=1:np %white subdomain boxes
        plot([t2((ii-1)*L/ds+1),t2((ii-1)*L/ds+1),t2(ii*L/ds),t2(ii*L/ds),t2((ii-1)*L/ds+1)],[t2((ii-1)*L/ds+1),t2(ii*L/ds),t2(ii*L/ds),t2((ii-1)*L/ds+1),t2((ii-1)*L/ds+1)],'g--');
    end
    figure(6);box on;hold on;grid on;title('correlation before conditioning')
    imagesc(t2,t2,K_stick);colormap(cmap);colorbar;caxis([min(C),max(C)]);
    xlabel('t [-]');ylabel('t [-]');
    xlim([-inf inf]);ylim([-inf inf]);set(gca,'xtick',0:L*np/5:L*np);set(gca,'ytick',0:L*np/5:L*np)
    for ii=1:np %white subdomain boxes
        plot([t2((ii-1)*L/ds+1),t2((ii-1)*L/ds+1),t2(ii*L/ds),t2(ii*L/ds),t2((ii-1)*L/ds+1)],[t2((ii-1)*L/ds+1),t2(ii*L/ds),t2(ii*L/ds),t2((ii-1)*L/ds+1),t2((ii-1)*L/ds+1)],'g--');
    end
    
    figure(7);hold on;box on;grid on;grid minor
    plot(t,K2(1,1:L/ds)-Kj);ylabel('Error at the jonction [-]')
   
    figure(8);box on;grid on;hold on
    F_stick(L/ds*(1:np-1),:)=nan;
    ii=1;%randi(N);
    plot(F_stick(:,ii),'k');plot(F(:,ii),'r--');
    set(gca,'xtick',0:25/ds:L/ds*np)
    aa=get(gca,'xticklabel');
    for ii=1:numel(aa);aa{ii}='';end
    aa{1}='0';
    for ii=1:np
       aa{4*ii+1}= [num2str(ii),'L'] ;
    end
    set(gca,'xticklabel',aa);        
    set(gca,'yticklabel',{})
    xlabel('s [-]'),ylabel('f [-]')
    
    figure(9);box on; title('conditioning matrix')
    imagesc_grid(kk);colorbar;caxis([-2.5,2.5]/10);grid off
    set(gca,'ydir','normal');

% % %     if savefig_on
% % %        figure(8)
% % %        savefig(['./fig2/',corr_model,'_ex'])
% % %        export_fig(gcf,['./fig2/',corr_model,'_ex'],'-pdf','-transparent')
% % %        figure(9)
% % %        savefig(['./fig2/',corr_model,'_kk'])
% % %        export_fig(gcf,['./fig2/',corr_model,'_kk'],'-png','-transparent')  
% % %     end
end





