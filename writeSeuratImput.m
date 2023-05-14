function writeSeuratImput( CAT,filename,fieldname )
%WRITESEURATIMPUT Summary of this function goes here
%   Detailed explanation goes here
if nargin < 3
    fieldname = 'data'
end
counts=CAT.(fieldname);
gene=CAT.gene;
Cellid=CAT.Cellid;
fid=fopen(filename,'w+');
CC='gene';
[x,y]=size(counts);
for i=1:y
    CC=[CC,',',Cellid{i}];
end
CC=[CC,char(10)];
CCC=cell(x,1);
parfor i=1:x;
    qwert=gene{i};
    for j=1:y
        qwert=[qwert,',',num2str(counts(i,j),10)];
    end
    qwert=[qwert,char(10)];
    CCC{i}=qwert;
end
fwrite(fid,CC);
for i=1:x
    fwrite(fid,CCC{i});
end

end

