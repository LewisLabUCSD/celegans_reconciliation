function genes = getallgenes_KEGG(data,org_code,remove_noKO)
% genes = getallgenes_KEGG(org_code,remove_noKO,web_options)
% finds all the genes which belong to the organism in KEGG
%
% INPUT:
% data: unrefined data obtained from KEGG API using webread
% remove_noKO: true or false, removes any gene with no KO assigned.
% web_options: options for webread
%
% OUTPUT:
% genes: column1-gene entry, column2-KO, column3-gene name, and
% column4-definition
%
ix = strfind(data,strcat(org_code,':'));
gene_list = cell(length(ix),1);
genes = cell(length(ix),1);
for i=1:length(ix)
    if i==length(ix)
        gene_list{i,1} = data(ix(i)+4:end);
    else
        gene_list{i,1} = data(ix(i)+4:ix(i+1)-1);
    end
    j = regexp(gene_list{i,1},'\t');
    for k = 1:length(j)+1
        if k==1
            genes{i,k} = strtrim(gene_list{i,1}(1:j(k)));
        elseif k==length(j)+1
            genes{i,k} = strtrim(gene_list{i,1}(j(k-1)+1:end));
        else
            genes{i,k} = strtrim(gene_list{i,1}(j(k-1)+1:j(k)));
        end
    end
    genes(i,1) = strrep(genes(i,1),'CELE_','');
    i1 = regexp(genes{i,2},' \| ');
    genes{i,3} = strtrim(genes{i,2}(i1+3:end));
    genes{i,2} = strtrim(genes{i,2}(1:i1));
    i1 = regexp(genes{i,3},';');
    genes{i,4} = strtrim(genes{i,3}(i1+1:end));
    genes{i,3} = strtrim(genes{i,3}(1:i1-1));
    genes{i,3} = strtrim(strrep(genes{i,3},'(RefSeq)',''));
end
if remove_noKO
    j = strcmp(genes(:,2),'no KO assigned');
    genes(j,:)=[];
end


