function [model,reconcile] = reconcile_gene_associations(model,rxn)

rid1 = ismember(model.rxns,rxn(1));
rid2 = ismember(model.rxns,rxn(2));

genes1 = find(model.rxnGeneMat(rid1,:)~=0);
genes2 = find(model.rxnGeneMat(rid2,:)~=0);

gene_list = setdiff(genes2,genes1);

if isempty(strfind(model.grRules{rid1,1},'and')) && isempty(strfind(model.grRules{rid2,1},'and'))
    for i=1:length(gene_list)
        grRules = strcat(model.grRules(rid1),{' or '},model.genes(gene_list(i)));
        rules = strcat(model.rules(rid1),{' or x('},num2str(gene_list(i)),{')'});
        model.grRules{rid1,1} = grRules{1,1};
        model.rxnGeneMat(rid1,gene_list(i)) = 1;
        model.rules{rid1,1} = rules{1,1};
    end
    reconcile = true;
elseif isempty(strfind(model.grRules{rid1,1},'or')) && isempty(strfind(model.grRules{rid2,1},'or'))
    for i=1:length(gene_list)
        grRules = strcat(model.grRules(rid1),{' and '},model.genes(gene_list(i)));
        rules = strcat(model.rules(rid1),{' and x('},num2str(gene_list(i)),{')'});
        model.grRules{rid1,1} = grRules{1,1};
        model.rxnGeneMat(rid1,gene_list(i)) = 1;
        model.rules{rid1,1} = rules{1,1};
    end
    reconcile = true;
else
    fprintf('\nGene reaction was not resolved.\n');
    reconcile = false;
end