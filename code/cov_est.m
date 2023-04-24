function COVs = cov_est(L,ds,Obs,flag)

nu=size(Obs,1);
%Obs=Obs-repmat(mean(Obs,2),[1 size(Obs,2)]);
%%%% STATIONARY ESTIMATOR
COVs=Obs'*Obs/(nu-1);
if flag %moyenne les diagonal( on est stationnarie)
    tmp=COVs;
    for n=1:L/ds
        tmp(length(tmp)*(n-1)+1:length(tmp)+1:numel(tmp))=sum(tmp(length(tmp)*(n-1)+1:length(tmp)+1:numel(tmp)))/(L/ds-n+1);
    end
    for n=2:L/ds
        tmp(n:length(tmp)+1:numel(tmp)-length(tmp)*(n-1))=sum(tmp(n:length(tmp)+1:numel(tmp)-length(tmp)*(n-1)))/(L/ds-n+1);
    end
    COVs=tmp; % this is the covariance function
end
% ind=find([0 (COVs(1,2:end).*COVs(1,1:end-1))<0]);
% for i=1:length(ind)
%     e=sum(COVs(1,ind(i)+1:end).^2)/sum(COVs(1,:).^2);
%     if normcdf(e*sqrt(nu))<0.975
%         k=ind(i);
%         break;
%     end
% end
% end
% k=max(k,100*4);
% tmp=COVs;
% for n=k:L/ds
%     tmp(length(tmp)*(n-1)+1:length(tmp)+1:numel(tmp))=0;
% end
% for n=k:L/ds
%     tmp(n:length(tmp)+1:numel(tmp)-length(tmp)*(n-1))=0;
% end
% COVs=tmp;
