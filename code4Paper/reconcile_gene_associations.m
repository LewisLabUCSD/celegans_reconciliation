function [model,reconcile] = reconcile_gene_associations(model,rxn,disp_flag)
% reconciles and changes gene properties within the model of two reactions
% written by chintan joshi on 2/8/2017 at ucsd

rid1 = ismember(model.rxns,rxn(1));
rid2 = ismember(model.rxns,rxn(2));

genes1 = find(model.rxnGeneMat(rid1,:)~=0);
genes2 = find(model.rxnGeneMat(rid2,:)~=0);

gene_list = setdiff(genes2,genes1);

% fprintf('Old Rule for %s:\t%s\n',model.rxns{rid1},model.grRules{rid1});

if isempty(strfind(model.grRules{rid1,1},'and')) && isempty(strfind(model.grRules{rid2,1},'and')) % if there is only "or" relationship in the entire rule of the reaction
    for i=1:length(gene_list)
        grRules = strcat(model.grRules(rid1),{' or ('},model.genes(gene_list(i)),')');
        rules = strcat(model.rules(rid1),{' | (x('},num2str(gene_list(i)),{'))'});
        model.grRules{rid1,1} = grRules{1,1};
        model.rxnGeneMat(rid1,gene_list(i)) = 1;
        model.rules{rid1,1} = rules{1,1};
    end
%     fprintf('New Rule for %s:\t%s\n',model.rxns{rid1},grRules{1,1});
    reconcile = true;
elseif isempty(strfind(model.grRules{rid1,1},'or')) && isempty(strfind(model.grRules{rid2,1},'or')) % if there is only "and" relationship in the entire rule of the reaction
    for i=1:length(gene_list)
        grRules = strcat(model.grRules(rid1),{' and ('},model.genes(gene_list(i)),')');
        rules = strcat(model.rules(rid1),{' & (x('},num2str(gene_list(i)),{'))'});
        model.grRules{rid1,1} = grRules{1,1};
        model.rxnGeneMat(rid1,gene_list(i)) = 1;
        model.rules{rid1,1} = rules{1,1};
    end
%     fprintf('New Rule for %s:\t%s\n',model.rxns{rid1},grRules{1,1});
    reconcile = true;
else
    if disp_flag
        fprintf('\nGene reaction was not resolved.\n');
    end
%     fprintf('New Rule for %s:\t%s\n',model.rxns{rid1},'Same as old rule');
    reconcile = false;
end