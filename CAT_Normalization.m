function CAT = CAT_Normalization( CAT )
%CAT_NORMALIZATION Summary of this function goes here
%   Detailed explanation goes here
[C,IC]=unique(CAT.gene);
[x,y]=size(CAT.data);
cdslenth=zeros(numel(C),1);
counts=zeros(numel(C),y);
for i=1:numel(C);
qwert=find(strcmp(CAT.gene,C{i}));
if numel(qwert)>1;
counts(i,:)=max(CAT.data(qwert,:));
cdslenth(i,:)=max(CAT.cdslenth(qwert,:));
else
counts(i,:)=(CAT.data(qwert,:));
cdslenth(i,:)=(CAT.cdslenth(qwert,:));
end;
end
CAT.gene=C;
CAT.data=counts;
CAT.cdslenth=cdslenth;
maskc=(sum(CAT.data)<CAT.par.Normalize.maxcounts).*(sum(CAT.data)>CAT.par.Normalize.mincounts).*(sum(CAT.data>0)<CAT.par.Normalize.maxgenes).*(sum(CAT.data>0)>CAT.par.Normalize.mingenes).*(CAT.mitoPer<CAT.par.Normalize.mitoPer)>0;
CAT.Cellid=CAT.Cellid(maskc);
CAT.counts=CAT.data;
if ~isempty(CAT.label)
CAT.label=CAT.label(maskc);
end
CAT.data=CAT.data(:,maskc);
switch CAT.par.Normalize.Type
    case {'CPM'}
        CAT.data=CAT.data./(ones(size(CAT.data,1),1)*sum(CAT.data))*1000000;
    case {'TPM'}
        CAT.data=CAT.data./(CAT.cdslenth*sum(CAT.data))*100000000;
    otherwise
        CAT.data=CAT.data;
end
switch CAT.par.Normalize.transform
    case 'log2'
        CAT.data=log2(CAT.data+1);
    case 'log10'
        CAT.data=log10(CAT.data+1);   
    otherwise
        CAT.data=CAT.data;
end
switch numel(CAT.par.Normalize.maskr)
    case 0
        maskr=sum(CAT.data>2,2)>2;
    case 1
        maskr=sum(CAT.data>CAT.par.Normalize.maskr,2)>CAT.par.Normalize.maskr;
    case 2
        maskr=sum(CAT.data>CAT.par.Normalize.maskr(1),2)>CAT.par.Normalize.maskr(2);
end
CAT.data=CAT.data(maskr,:);
CAT.gene=CAT.gene(maskr);
CAT.counts=CAT.counts(maskr,maskc);

end

