function [CATn,CAT] = preprocessingData( project, para, outputname )
%PREPROCESSING this function was used to precess data include quality
%control, feature selection, t-SNE transformation, and export files
%   Detailed explanation goes here

% import counts data and annotation
if nargin<3
    outputname = 'Seurat.csv'
end
filePath=[''];
CAT=CAT_Initial([filePath,'originalcounts.csv'],para.species);
CAT=CAT_ImportExpression(CAT);
load(['Mouse_anno_filtered'],'anno_withRibo','anno_withoutRibo','anno_mito');
CAT.par.Normalize.Type=para.normalizeType;
try 
    if (para.mitoPer~=1)&(~isnan(para.mitoPer))
        CATm = CAT_UniformName( CAT,anno_mito );
        CAT.mitoPer=sum(CATm.data)./sum(CAT.data);
    end
catch
    CAT.mitoPer = zeros(1,numel(CAT.Cellid));
end
if para.withRibo==0
    CAT = CAT_UniformName( CAT,anno_withoutRibo );
else
    CAT = CAT_UniformName( CAT,anno_withRibo );
end
% quality control and calculate expression matrix
CAT.par.Normalize.maxcounts=para.maxcounts;
CAT.par.Normalize.mincounts=para.mincounts;
CAT.par.Normalize.maxgenes=para.maxgenes;
CAT.par.Normalize.mingenes=para.mingenes;
CAT.par.Normalize.mitoPer=para.mitoPer;
CAT.par.Normalize.maskr=para.maskr;
CATn = CAT;
CATn.mitoPer
CATn = CAT_Normalization( CAT );

writeSeuratImput( CATn, [filePath,outputname]); % export preprocessed matrix to csv file
end
