function model = update_gene_properties(model,gene_rep,rem_flag)
% model = update_gene_properties(model,gindex)
% Changes only gene properties.

% INPUT:
% model: model to be changed
% index: index, max(index) <= size(model.genes,1)
%
% OUTPUT:
% model: changed model

if rem_flag
    % % change gene properties from gene 1 to gene 2 and remove gene 1
    x_old = regexprep(strcat({'x'},num2str([1:1:length(model.genes)]')),' ','');
    x_old(ismember(model.genes,gene_rep(:,1))) = [];
    for i=1:size(gene_rep,1)
%         fprintf('Changing %s to %s.\n',gene_rep{i,1},gene_rep{i,2});
        model.rxnGeneMat(:,strcmp(model.genes,gene_rep{i,2})) = any(model.rxnGeneMat(:,strcmp(model.genes,gene_rep{i,1})) + ...
            model.rxnGeneMat(:,strcmp(model.genes,gene_rep{i,2})),3);
        model.grRules = regexprep(model.grRules,gene_rep{i,1},gene_rep{i,2});
        model.rxnGeneMat(:,strcmp(model.genes,gene_rep{i,1})) = [];
        model.genes(strcmp(model.genes,gene_rep{i,1})) = [];
        
    end
    x_new = regexprep(strcat({'x'},num2str([1:1:length(model.genes)]')),' ','');
    for i=1:length(x_new)
        model.rules = regexprep(model.rules,x_old{i,1},x_new{i,1});
    end
else
    % % change gene name from gene 1 to gene 2
    for i=1:size(gene_rep,1)
        model.grRules = regexprep(model.grRules,model.genes(strcmp(model.genes,gene_rep(i,1))),gene_rep(i,2));
        model.genes(strcmp(model.genes,gene_rep(i,1))) = gene_rep(i,2);
%         fprintf('%d. Total  = %d, unique = %d\n',i,length(model.genes),unique(length(model.genes)));
    end
end
