function CAT = CAT_Initial(filename,species)
%CAT_Initial Initial CAT object and set default parameters
%   CAT = CAT_Initial(filename,species) returns an object CAT with
%   parameters required in further manipulation. filename should be a table
%   file of gene expression matrix. This file could be txt, xls, xlsx, csv
%   and any other format but should give the seperons. There must be a unit
%   in the first column having a string 'gene' to show the gene list will
%   be start in the next row. species is the organism used in the
%   exprenments.
%
%   attributes of CAT:
%   CAT.filename:       Origin file of gene expression matrix
%   CAT.splittype:      seperon in the file, defaults are according to file
%                       type
%   CAT.par:            parameter required in further analysis
%       .annotation:    parameter for filtering genes by annotation
%           .species:   organism species used in the experiments
%           .selectTranscriptType:  
%                       whether filtered records according to transcript
%                       type. Default: 1
%           .selectSource: 
%                       whether filtered records according to their source.
%                       Default: 1
%           .excludeGene:
%                       whether remove records according to their special
%                       description. Default: 1
%           .Unitype:   According to what for combining redundant records.
%                       Default: 'Gene'
%           .genetype:  Will be detected automatically according to feature
%                       name
%           .typeselect:
%               .selectTranscriptType:
%                       aim transcript type
%               .selectSource:
%                       aim record source
%               .excludedType:
%                       remove genes with given discription
%       .Normalize:     parameter for quality control and normalization
%               .maxcounts/.mincounts:
%                       upper/lower limit of counts number
%               .maxgenes/.mingenes:
%                       upper/lower limit of feature number
%               .Type:  normalization methods.('CPM'/'TPM')
%               .maskr: an vector with 2 elements [a,b]. It means more than
%                       a cells with CPM/TPM larger than b
%               .transform:
%                       whether to transformation. Default:'Null'
%   CAT.data:           Count matrix(Before normalization)/Expression
%                       matrix(After normalization)
%   CAT.counts:         Count matrix (After normalization)
%   CAT.gene:           gene symbol
%   CAT.label:          metadata of data
%   CAT.Cellid:         Cell identification in expression files
%   CAT.orinalID:       feature name before symbol transformation
%   CAT.cdslenth:       Length of cds region
%
%Excample:
%   CAT = CAT_Initial(filename,species)

CAT.filename=filename;
CAT.splittype=[];
CAT.par.annotation.species=species;
CAT.par.annotation.selectTranscriptType=1;
CAT.par.annotation.selectSource=1;
CAT.par.annotation.excludeGene=1;
CAT.par.annotation.Uniontype='Gene'
CAT.par.Normalize.maxcounts=inf;
CAT.par.Normalize.mincounts=-inf;
CAT.par.Normalize.maxgenes=inf;
CAT.par.Normalize.mingenes=-inf;
CAT.par.Normalize.Type='CPM'
CAT.par.Normalize.maskr=[];
CAT.par.Normalize.transform='Null';
end

