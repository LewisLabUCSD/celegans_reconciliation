function model = update_gene_properties(model,gene_rep,rem_flag)
% model = update_gene_properties(model,gindex)
% Changes only gene properties.

% INPUT:
% model: model to be changed
% index: index, max(index) <= size(model.genes,1)
%
% OUTPUT:
% model: changed model
% model.modelVersion
% find(ismember(model.genes,gene_rep(:,1)))
% find(ismember(model.genes,gene_rep(:,2)))
if rem_flag
    % % change gene properties from gene 1 to gene 2 and remove gene 1
    for i=1:size(gene_rep,1)
        x_old = regexprep(strcat({'x\('},num2str([1:1:length(model.genes)]'),'\)'),' ','');
        n_new = [1:1:length(model.genes)]';
        n_new(ismember(model.genes,gene_rep(i,1))) = n_new(ismember(model.genes,gene_rep(i,2))); % reassign genes
        irem = find(ismember(model.genes,gene_rep(i,1))); ikeep = find(ismember(model.genes,gene_rep(i,2)));
        fprintf('Changing %s to %s.\n',gene_rep{i,1},gene_rep{i,2});
        % % fixing rxnGeneMat
        model.rxnGeneMat(:,strcmp(model.genes,gene_rep{i,2})) = any(model.rxnGeneMat(:,strcmp(model.genes,gene_rep{i,1})) + ...
            model.rxnGeneMat(:,strcmp(model.genes,gene_rep{i,2})),3);
        model.rxnGeneMat(:,strcmp(model.genes,gene_rep{i,1})) = [];
        n_new(n_new > irem) = n_new(n_new > irem) - 1;
        x_new = regexprep(strcat({'g\('},num2str(n_new),'\)'),' ','');
        model.rules = regexprep(model.rules,x_old,x_new); % update rules field
        model.rules = regexprep(model.rules,'g','x');
        model.genes(strcmp(model.genes,gene_rep{i,1})) = []; % remove genes
    end
else
    % % change gene name from gene 1 to gene 2
    for i=1:size(gene_rep,1)
        model.genes(strcmp(model.genes,gene_rep(i,1))) = gene_rep(i,2);
    end
end
% update grRules field
model = grRulesModel(model);