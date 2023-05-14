%%
%%Parse data from files downloaded from ThymusE1012. Write a csv file
%%following the required formats of futher procedure
clear all
filename='rawcounts.csv'
D = importdata(filename);
gene = D.textdata(2:end,1);
data = D.data;
[x,y]=size(data);
CC = 'gene'
for i=1:y
    CC = [CC,',',D.textdata{1,i+1}];
end
CC=[CC,char(10)];
parfor i=1:x
    qwert = gene{i};
    for j=1:y;
        qwert = [qwert,',',num2str(data(i,j),10)];
    end
    qwert = [qwert,char(10)];
    CCC{i}=qwert;
end

fid = fopen('originalcounts.csv','w+');
fwrite(fid,CC);
for i=1:x
    fwrite(fid,CCC{i});
end
fclose(fid)
clear all

%%
% Preprocessing of ThymusE1012
para=initialPara;
para.species='Mus';
para.normalizeType='CPM';
para.maxcounts=400000;
para.mincounts=40000;
para.maxgenes=5000;
para.mingenes=500;
para.maskr=[2,2];
para.perplex=30;
para.withRibo=0;
para.mito=50;
[CATn,CAT] = preprocessingData( 'ThymusE1012', para )
save('ThymusE10.mat')