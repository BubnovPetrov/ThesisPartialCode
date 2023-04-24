clc;clear;close all;
set(0,'DefaultFigureColor',[1 1 1]);set(0,'DefaultTextInterpreter','latex');set(0,'defaultlinelinewidth',4);set(0,'DefaultAxesFontSize',20)
set(0,'DefaultFigureColormap',feval('jet'));close
rng(10)

%% Data
L=10; % subdomain side length
ds=1; % KLE discretization spacing
lc1=L*0.5;lc2=L*0.5; % correlation lengths (with respect to subdomain)
N=1000; % samples of random field 
np=5; % number of domain partition (number of subdomains each dimension)
tens=1; % tensorizable corr fcn?
corr_model{1}='exp';
corr_model{2}='exp';
dist_norm='L2'; 
err_kl=1e-1; % KLE accuracy of representation
elaps=(zeros(1,3)); % elapsed times [KLE decomp, prologation/conditioning, random field generation]

%% For plotting
plot_on=1;
video_on=0; video_intermediate=0;filename_gif = 'aaa.gif';
cmap=jet;cmap=flipud([cmap(:,1),cmap(:,2),cmap(:,3)]);close;
%% Correlation
t=0:ds:L-ds;
t2=0:ds:2*L-ds;
for cc=1:2
    switch corr_model{cc}
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
    cfun{cc}=corrfun; %#ok<*SAGROW>
end
if tens;
    Cx=feval(cfun{1},t,lc1);  Cx(1)=1;  
    Cy=feval(cfun{2},t,lc2);  Cy(1)=1;
    C2x=feval(cfun{1},t2,lc1); C2x(1)=1;
    C2y=feval(cfun{2},t2,lc2); C2y(1)=1;
    C2=C2x'*C2y; % tensor product of the 1D correlation generated above
    C=Cx'*Cy; % tensor product of the 1D correlation generated above
else
    switch dist_norm
        case 'L2'
            C=feval(cfun{1},sqrt( (ones(L/ds,1)*t/lc1).^2 +(t'/lc2*ones(1,L/ds)).^2 ),1);C(1)=1;
            C2=feval(cfun{1},sqrt( (ones(L/ds*2,1)*t2/lc1).^2 +(t2'/lc2*ones(1,L/ds*2)).^2 ),1);C2(1)=1;
        case 'L1'
            C=feval(cfun{1},abs( (ones(L/ds,1)*t/lc1) +(t'/lc2*ones(1,L/ds)) ),1);C(1)=1;
            C2=feval(cfun{1},abs( (ones(L/ds*2,1)*t2/lc1) +(t2'/lc2*ones(1,L/ds*2)) ),1);C2(1)=1;
        otherwise
            C=feval(cfun{1},sqrt( (ones(L/ds,1)*t/lc1).^2 +(t'/lc2*ones(1,L/ds)).^2 ),1);C(1)=1;
            C2=feval(cfun{1},sqrt( (ones(L/ds*2,1)*t2/lc1).^2 +(t2'/lc2*ones(1,L/ds*2)).^2 ),1);C2(1)=1;
    end
end

%% KL decomposition
tic;
if ~tens
    CC=corr2D(C);
    [phi,lam]=eig(CC); [lam,ord]=sort(diag(lam),'descend');phi=phi(:,ord);
    lam(lam<0)=0; %%just in case lam=-eps
    nkl=find(1-cumsum(lam)/sum(lam)<err_kl,1);
    lam=lam(1:nkl);phi=phi(:,1:nkl);
else
    %Tensorisation; if corr1=corr2 one eigenval can be associated with 2 eigvect
    [p1,l1]=eig(toeplitz(Cx));[l1,ord]=sort(diag(l1),'descend');p1=p1(:,ord); % discrete solution for KLE eigenpairs in x direction
    [p2,l2]=eig(toeplitz(Cy));[l2,ord]=sort(diag(l2),'descend');p2=p2(:,ord); % discrete solution for KLE eigenpairs in y direction
    
    ll=l1*l2'; [ll,ind]=sort(ll(:),'descend');
    nkl0=find(1-cumsum(ll)/sum(ll)<err_kl,1);
    ll=ll(1:nkl0);phi=zeros((L/ds)^2,nkl0); % finding tensor product of eigenvalues upto certain KLE accuracy
    % finding the corresponding row and col related to the final
    % eigenvalues, as row = x mode and col = y mode (tensor product wise).
    [I,J] = ind2sub([L/ds L/ds],ind); IJ=[I(1:nkl0),J(1:nkl0)]; %indices des modes
    % IJ are the mode pairs of tensor product eigenvalues
    for ii=1:numel(ll);
        tmp=p1(:,I(ii))*p2(:,J(ii))'; % tensor product of the 2 mode eigenvectors
        phi(:,ii)=tmp(:);
    end
    lam=ll;nkl=nkl0;
end

elaps(1)=toc;
disp(elaps)

if plot_on;
    figure(1);set(gcf,'units','normalized','position',[0,0,1,1]);box on; hold on;grid on;grid minor
    for mm=1:25
        if mm<=size(phi,2)
            subplot(5,5,mm);imagesc(reshape(phi(:,mm),[L/ds,L/ds]));set(gca,'xtick',[],'ytick',[]);colormap(cmap);pbaspect([1 1 1]);
        end
    end
end

%% Prolongation - germs
tic;
% producing the X and L coupling matrices required for random field
% correlation extension. All are in the LL cell vector
if tens
    LL=kl_prolongation2D_tens(nkl,C2x,C2y,p1,p2,l1,l2,IJ); %%matrice prolong
else
    LL=kl_prolongation2D(C2,lam,phi); %%KL 
end
elaps(2)=toc;
disp(elaps)

tic;


F2=nan(L/ds*np,L/ds*np,N); % random field samples 
etag=randn(numel(lam),np,np,N); % random field coefficients for each sample
%%%jj = direction x; ii = direction y





% 0 - They Don't change (no need for condition)
for ii=1:2:np; % y position coeff
    for jj=1:2:np % x position coeff
        tmp=phi*sqrt(diag(lam))*squeeze(etag(:,jj,ii,:)); %generations
        F2((jj-1)*L/ds+1:jj*L/ds,(ii-1)*L/ds+1:ii*L/ds,:)=reshape(tmp,[L/ds,L/ds,N]);
    end
end




% 1 - Along x
for ii=1:2:np
    for jj=2:2:np
        if video_on && video_intermediate
            F3((jj-1)*L/ds+1:jj*L/ds,(ii-1)*L/ds+1:ii*L/ds,:)=reshape(phi*sqrt(diag(lam))*squeeze(etag(:,jj,ii,:)),[L/ds,L/ds,N]);
        end
        etag(:,jj,ii,:)=LL{1}*squeeze([etag(:,jj-1,ii,:);etag(:,jj+1,ii,:);etag(:,jj,ii,:)]);
        tmp=phi*sqrt(diag(lam))*squeeze(etag(:,jj,ii,:)); %generations
        F2((jj-1)*L/ds+1:jj*L/ds,(ii-1)*L/ds+1:ii*L/ds,:)=reshape(tmp,[L/ds,L/ds,N]);
    end
    
end
%%


%%
% 2 - Along y

for ii=2:2:np
    for jj=1:2:np
        if video_on && video_intermediate; F3((jj-1)*L/ds+1:jj*L/ds,(ii-1)*L/ds+1:ii*L/ds,:)=reshape(phi*sqrt(diag(lam))*squeeze(etag(:,jj,ii,:)),[L/ds,L/ds,N]);end
        if jj==1 % subdomain near x = 0, case of 4 neighbors
            tmp=[etag(:,jj,ii-1,:);etag(:,jj+1,ii-1,:);...
                etag(:,jj,ii+1,:); etag(:,jj+1,ii+1,:);...
                etag(:,jj,ii,:)];
            etag(:,jj,ii,:)=LL{3}*squeeze(tmp);
        elseif jj==np % subdomain near x = end, case of 4 neighbors
            tmp=[etag(:,jj-1,ii-1,:);etag(:,jj,ii-1,:);...
                etag(:,jj-1,ii+1,:); etag(:,jj,ii+1,:);...
                etag(:,jj,ii,:)];
            etag(:,jj,ii,:)=LL{4}*squeeze(tmp);
        else % subdomain in middle, case of 6 neighbors
            tmp=[etag(:,jj-1,ii-1,:);etag(:,jj,ii-1,:);etag(:,jj+1,ii-1,:);...
                etag(:,jj-1,ii+1,:); etag(:,jj,ii+1,:);etag(:,jj+1,ii+1,:);...
                etag(:,jj,ii,:)];
            etag(:,jj,ii,:)=LL{2}*squeeze(tmp);
            %            etag(:,jj,ii,:)=LL{6}*squeeze([etag(:,jj,ii-1,:);etag(:,jj,ii+1,:);etag(:,jj,ii,:)]);
        end
        
        tmp=phi*sqrt(diag(lam))*squeeze(etag(:,jj,ii,:)); %generations
        F2((jj-1)*L/ds+1:jj*L/ds,(ii-1)*L/ds+1:ii*L/ds,:)=reshape(tmp,[L/ds,L/ds,N]);
        
    end
end
%%

%%
% 3- Filling (Surrounded cells) ---> case of 8 neighbors
if video_intermediate;F3=F2;end
for ii=2:2:np
    for jj=2:2:np
        if video_on && video_intermediate; F3((jj-1)*L/ds+1:jj*L/ds,(ii-1)*L/ds+1:ii*L/ds,:)=reshape(phi*sqrt(diag(lam))*squeeze(etag(:,jj,ii,:)),[L/ds,L/ds,N]);end
        tmp=[etag(:,jj-1,ii-1,:);etag(:,jj,ii-1,:);etag(:,jj+1,ii-1,:);...
            etag(:,jj-1,ii+1,:);etag(:,jj,ii+1,:);etag(:,jj+1,ii+1,:);...
            etag(:,jj-1,ii,:)  ;etag(:,jj+1,ii,:);                    ...
            etag(:,jj,ii,:)];
        etag(:,jj,ii,:)=LL{5}*squeeze(tmp);
        
        tmp=phi*sqrt(diag(lam))*squeeze(etag(:,jj,ii,:)); %generations
        F2((jj-1)*L/ds+1:jj*L/ds,(ii-1)*L/ds+1:ii*L/ds,:)=reshape(tmp,[L/ds,L/ds,N]);
    end
    
end

elaps(3)=toc;

disp(elaps)

%% Writing video and plot

% plotting the first min(9, samples) random field generated

if plot_on
    figure(3);set(gcf,'units','normalized','position',[0,0,1,1]);
    nn=min(9,N);
    for mm=1:nn
        subplot(round(sqrt(nn)),ceil(sqrt(nn)),mm);box on;hold on;
        imagesc(F2(:,:,mm));set(gca,'xtick',[],'ytick',[]);xlim([-inf,inf]);ylim([-inf,inf]);colormap(cmap);
        caxis([-4.5,4.5])
        for ll=1:np-1
            plot([L/ds*ll,L/ds*ll],[1,L/ds*np],'w--','linewidth',2)
            plot([1,L/ds*np],[L/ds*ll,L/ds*ll],'w--','linewidth',2)
        end
    end
end





%% check mean and variance

for i = 1:N
    mu(i) = mean(F2(:,:,i),'all');
    varsamp(i) = var(F2(:,:,i),0,'all');
end

size(F2)
size(etag)
mean(mu)
mean(sqrt(varsamp))

