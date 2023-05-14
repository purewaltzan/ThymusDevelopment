function CAT = CAT_ImportExpression(CAT)
%CAT_ImportExpression Import expression matrix from files
%   CAT = CAT_ImportExpression(CAT) returns an object CAT with adding
%   CAT.data, CAT.gene, CAT.Cellid, CAT.label, CAT.labelname according to
%   input files. The feature name will be detected and output to 
%   CAT.par.annotation.genetype
%
%   CAT must have CAT.filename attribute
%
%Example:
%   CAT = CAT_ImportExpression(CAT)

%Reading files by row
D = importdata(CAT.filename,'\v');
ad=regexp(CAT.filename,'\.','split');
for i = 1:numel(D)
    qwert=D{i};
    D{i}=qwert(qwert~='"');
end

%Determine the seperon according to file type
if isempty(CAT.splittype)
    switch ad{end}
        case 'csv'
            CAT.splittype = ',';
        case {'xlsx','xls','tsv','txt'}
            CAT.splittype = '\t' ;
        otherwise
            error('CAT only support .csv, Excel, .txt formation or any files which was given vertical and horizantal split charactors');
    end
end

%Parse data to matrix and write to CAT object
splittype=CAT.splittype;
switch iscell(D)*10+isstruct(D)
    case 10
        datastart=findstart(D)+1;        
        CAT.data=dlmread(CAT.filename,CAT.splittype,datastart-1,1);
        [x,y]=size(CAT.data) ;
        gene=cell(x,1);        
        for i=datastart:numel(D)
            qwe=regexp(D{i},splittype,'split');
            gene(i-datastart+1)=qwe(1);
        end
        CAT.gene=gene;      
        if datastart>2
            for i=1:datastart-2
                l(i,:)=regexp(D{i},CAT.splittype,'split');
            end
            for j=1:y
                wes=['CAT.label(',num2str(j),')=struct('];
                for k=1:datastart-2
                    wes=[wes,'l{',num2str(k),',1},l{',num2str(k),',',num2str(j+1),'},'];
                end
                wes=[wes(1:end-1),');'];
                eval(wes);
            end
        else
            CAT.label=[];
            CAT.labelname=[];
        end
        qwe=regexp(D{datastart-1},CAT.splittype,'split');
        CAT.Cellid=qwe(2:end);
        
    case 1
        if isstruct(D.data)
            fe=fieldnames(D.data);
            D.data=D.data.(fe{1});
            D.textdata=D.textdata.(fe{1});
        end
        datastart=find(strcmp(D.textdata(:,1),'gene'))+1;
        CAT.data=D.data;
        CAT.gene=D.textdata(datastart:end,1);
        CAT.Cellid=D.textdata(datastart-1,2:end);
        if datastart>2
           for i=1:datastart-2
               l=D.textdata(i,:);
               for j=1:y
                   CAT.label(j)=struct(l{1},l{j+1});
               end
           end
        else
            CAT.label=[];
            CAT.labelname=[];
        end       
            
        
end

%Determine feature type in imput files
if sum(CAT.gene{1}=='.')~=0;
    for i=1:numel(CAT.gene)
        qwert=CAT.gene{i};
        qwert=regexp(qwert,'\.','split');
        CAT.gene{i}=qwert{1};
    end
end
ee=CAT.gene{1};
if numel(ee)>13
    switch (ee(numel(ee)-11))
        case 'T'
            CAT.par.annotation.genetype='TranscriptID';
        case 'G'
            CAT.par.annotation.genetype='GeneID';
        otherwise
            CAT.par.annotation.genetype='Symbol';
    end
else
    CAT.par.annotation.genetype='Symbol';
end

end

function n=findstart(D);
%findstart find where is 'gene' in the table of imput files
n=1;
for i=1:numel(D)
    qwe=D{i};
    if ~strcmp(upper(qwe(1:4)),'GENE')
        n=n+1;
    else
    return
    end
end
end
