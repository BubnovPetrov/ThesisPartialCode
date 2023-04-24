function CC=corr2D(C,flag)

if ~exist('flag','var')
    flag=1;
end

c=cell(length(C),1);
[c{:}] = deal(zeros(length(C)));

for ii=1:length(C)
    c{ii}=toeplitz(C(:,ii)');
end

if flag;
    CC=cell2mat(c(toeplitz(1:numel(c))));
else 
    CC=c;
end
return

