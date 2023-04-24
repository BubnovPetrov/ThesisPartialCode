function [finelepair,finE2sim,finE4str,kmean,finvecvalmean] = storagecycle_v1_Eandnu(imat,jmat,kele,Esampled,samples,freedofs)        
% the simulation is done using input samples

% TheMean = 1; %mean of random field
% TheCOV = 0.35; %Guest cantilever literature cantilever (large var)

% RFinput.LNMean = TheMean  ;
% RFinput.LNStdv = TheCOV * RFinput.LNMean ; 
% FGStdv = sqrt(log((RFinput.LNStdv / RFinput.LNMean)^2 +1));
% FGMean =  log(RFinput.LNMean) - 0.5 * (FGStdv)^2;

% lc = [100,30]; %2D long correlation length for cantilever beam
% E0 = 1; % wheel

% generate E statistics 
parpool(4);
nelem = size(Esampled,1);
kmean = zeros(size(kele(:,:,:,1)));
parfor i = 1:samples
    kcur = kele(:,:,:,i);
    temp = zeros(size(kcur));
    for j = 1:nelem
        temp(:,:,j) = Esampled(j,i).*kcur(:,:,j);
    end
    kmean = kmean + temp;
end
kmean = kmean./samples;

%% trace(E[K^2]) storage
elemnum = ones(size(kele(:,:,:,1)));
for i = 1:size(kele,3)
    elemnum(:,:,i) = elemnum(:,:,i)*i;
end

% initial mapping
colj = jmat(:); rowi = imat(:); elemnum = elemnum(:); 
consts = reshape(kele,[numel(colj),size(kele,4)]);
consts = consts.*Esampled(elemnum,:);

clear Esampled 

%removing fixed D.O.F entries
freeidx = ismember(colj,freedofs)&ismember(rowi,freedofs);
colj = colj(freeidx);
rowi = rowi(freeidx);
consts = consts(freeidx(:),:);
elemnum = elemnum(freeidx);

% column sorting
[colj,tempindx] = sort(colj);
rowi = rowi(tempindx);
consts = consts(tempindx(:),:);
constm = mean(consts,2); %%%%
elemnum = elemnum(tempindx);

finvecval = zeros(1e6,samples);
finelepair = zeros(1e6,2);
finvecvalmean = zeros(1e6,1);

% traverse through sorted column, populate row entries
coljindx = unique(colj);

indxctr = 0; 

for j = 1:length(coljindx)
    % determine the non-zero rows in column coljindx(j) 
    tempi = rowi(colj == coljindx(j));
    tempcon = consts(colj == coljindx(j),:);
    tempconm = constm(colj == coljindx(j)); %%%%
    tempnum = elemnum(colj == coljindx(j));

    % sort the row index
    [tempi,tempindx] = sort(tempi);
    tempcon = tempcon(tempindx(:),:);
    tempconm = tempconm(tempindx(:)); %%%%
    tempnum = tempnum(tempindx);

    %indxset = zeros(18,2);
    eleminfoa = zeros(18,4,samples);
    eleminfob = zeros(18,4);
    ctr = 1; % num of unique row index corresponding to current column
    ctr2 = zeros(18,1); % num of same entry terms in entry (i,j)
    for k = 1:length(tempi) % traverse through all row entries for current column
        ctr2(ctr) = ctr2(ctr) + 1; % number of elements interacting with current d.o.f (i,j)
        %indxset(ctr,:) = [coljindx(j),tempi(k)]; % (col#, row#) tupple
        eleminfoa(ctr,ctr2(ctr),:) = tempcon(k,:); % constant
        eleminfom(ctr,ctr2(ctr)) = tempconm(k); %%%%
        eleminfob(ctr,ctr2(ctr)) = tempnum(k); % ele num
        % (maximum of 4 elements can interact at each d.o.f for 2D case)
        % (hence eleminfo's 3rd dimension is restricted to 4)

        if k < length(tempi) && tempi(k) ~= tempi(k+1) 
            % reset counter when you move on to the next entry
            ctr = ctr+1;
        end
    end

    for k = 1:ctr % traverse through the unique (i,j) entries
        curindxcon = eleminfoa(k,1:ctr2(k),:); 
        if size(curindxcon,1) ~= 1
            print('flag')
        end
        t = squeeze(curindxcon);%%%%
        curindxcon = reshape(t,[size(curindxcon,2),size(curindxcon,3)]);%%%%
        curindxconm = eleminfom(k,1:ctr2(k)); %%%%
        curindxnum = eleminfob(k,1:ctr2(k));

        for i = 1:ctr2(k) % traverse through contributing elements to current (i,j) entry
            for m = i:ctr2(k) 
                if i == m % identical element entry
                    %ccmean = mean(curindxcon(i,:))*mean(curindxcon(m,:));
                    ccmean = curindxconm(i)*curindxconm(m); %%%%
                    cc = curindxcon(i,:).*curindxcon(m,:);
                else % cross element entry, so multiplied by 2
                    %ccmean = mean(curindxcon(i,:))*mean(curindxcon(m,:))*2;
                    ccmean = curindxconm(i)*curindxconm(m)*2; %%%%
                    cc = curindxcon(i,:).*curindxcon(m,:).*2;
                end
            
                % only store the element pairs and values
                indxctr = indxctr + 1;
                finvecval(indxctr,:) = cc;
                finvecvalmean(indxctr) = ccmean;
                finelepair(indxctr,:) = [curindxnum(i),curindxnum(m)];
            end
        end
    end
end

clear consts constm elemnum colj rowi freeidx coljindx

finelepair = finelepair(1:indxctr,:);
finvecval = finvecval(1:indxctr,:);
finvecvalmean = finvecvalmean(1:indxctr);



%parpool(4);
finelepair = sort(finelepair,2);
pairele = unique(finelepair,'rows');
valele = zeros(length(pairele),samples);
valmeanele = zeros(length(pairele),1);
parfor i = 1:size(pairele,1)
    mask = (pairele(i,1) == finelepair(:,1))&(pairele(i,2) == finelepair(:,2));
    valele(i,:) = sum(finvecval(mask,:),1);
    valmeanele(i) = sum(finvecvalmean(mask));
end

finelepair = pairele;
finvecval = valele;
finvecvalmean = valmeanele;
clear pairele valele valmeanele



% generate E^2 statistics
finE4str = zeros(length(finvecvalmean),samples);
finE2sim = zeros(size(finvecvalmean));
parfor i = 1:size(finvecval,1)
    finE2sim(i) = sum(finvecval(i,:))/samples;
    finE4str(i,:) = finvecval(i,:);
end

delete(gcp('nocreate'));

end


%finvecval is the vector of constants
%finelepair is the vector element pair