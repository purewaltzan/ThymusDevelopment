function CAT = CAT_UniformName( CAT,anno )
%CAT_UNIFORMNAME Summary of this function goes here
%   Detailed explanation goes here
fieldsname = fields(anno.content);
CAT.par.annotation.genetype
switch CAT.par.annotation.genetype
    case 'TranscriptID'        
        if sum(strcmp(fieldsname,'Transcript_stable_ID'))==0;
            error('There did not exists Transcript stable ID column in ensembl annotation file')
        else
            KK='Transcript_stable_ID';

        end
    case 'GeneID'
        if sum(strcmp(fieldsname,'Gene_stable_ID'))==0;
            error('There did not exists Transcript stable ID column in ensembl annotation file')
        else
            KK='Gene_stable_ID';
        end
    case 'Symbol'
        if sum(strcmp(fieldsname,'Gene_name'))==0;
            error('There did not exists Gene name column in ensembl annotation file')
        else
            KK='Gene_name';
        end
end
annotationgene = getfield(anno.content,KK);
annotationgene=upper(annotationgene);
[C,IC,IA]=intersect(upper(CAT.gene),annotationgene);
for i=1:numel(C)
    L=0;
    qwert=find(strcmp(C{i},upper(getfield(anno.content,KK))));
    if numel(qwert)> 0
        L = anno.content.Transcript_length_including_UTRs_and_CDS(qwert);
    else
        L = 0;
    end
    cdsLength(i)=max(L);
end
CAT.orinalID=CAT.gene(IC);
qwert=getfield(anno.content,'Gene_name');
CAT.gene=qwert(IA,:);
CAT.data=CAT.data(IC,:);
CAT.cdslenth=cdsLength';
[aaa,bbb]=sort(CAT.gene);
CAT.gene=CAT.gene(bbb,:);
CAT.data=CAT.data(bbb,:);
CAT.orinalID=CAT.orinalID(bbb,:);
CAT.cdslenth=CAT.cdslenth(bbb,:);
CAT.par.annotation.typeselect=rmfield(anno,'content');
end

