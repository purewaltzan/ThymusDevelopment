function para = initialPara
%INITIALPARA Initialize parameters for preprocessing 

% default paramaters of dataset
para.species='Mus';
para.normalizeType='TPM';

% default quality control paramater
para.maxcounts=inf;
para.mincounts=0;
para.maxgenes=inf;
para.mingenes=0;
para.mitoPer=0.5;
para.maskr=[2,2];

% default perpelx for t-SNE transformation
para.perplex=30;

end

