function [impr,rem_rxns,gene_list,gra_index] = reconcile_gene_reaction(impr)
% [impr,rem_rxns] = reconcile_gene_reaction(impr)
% updates the gra to the version that has the larger coverage and removes
% the duplicate

% INPUT:
% impr: improvement strcuture
% OUTPUT:
% impr: new improvement structure
% rem_rxns: reactions removed as a consequence of updating GRA
%
% Comment: Assumption is that bigger gene converage of the reaction is more
% accurate. This can't always be true, but works for now.

diff_gra = impr.diff_gra;
model = impr.model;
gene_list = cell(size(diff_gra,1),1);
rem_rxns = cell(1,1);
% find the number of genes in original and duplicate
c = 0;
general_gene = {'NA';'Unknown';'TBD';'ND'};
for i=1:size(diff_gra,1)
    gene_list{i,1} = model.genes(model.rxnGeneMat(strcmp(model.rxns,diff_gra{i,1}),:)~=0);
    gene_list{i,2} = model.genes(model.rxnGeneMat(strcmp(model.rxns,diff_gra{i,2}),:)~=0);
    if isequal(union(gene_list{i,1},gene_list{i,2}),gene_list{i,1}) 
        c = c+1;
        rem_rxns{c,1} = diff_gra{i,1}; rem_rxns{c,2} = diff_gra{i,2}; 
        gra_index(c,1) = i;
    elseif isequal(union(gene_list{i,1},gene_list{i,2}),gene_list{i,2}) ...
            || (~isempty(find(ismember(gene_list{i,1},general_gene))) && ~isempty(gene_list{i,2}))
        c = c+1;
        rem_rxns{c,1} = diff_gra{i,2}; rem_rxns{c,2} = diff_gra{i,1};
        gra_index(c,1) = i;
    end
end
rem_index = find(ismember(model.rxns,rem_rxns(:,2))); % finding the index for removal
gene_list(gra_index,:) = [];
model1 = update_reaction_properties(model,rem_index); % removing from the model
impr.rem_rxns = [impr.rem_rxns;rem_rxns]; % adding to removed reactions
impr.model = model1;
impr.diff_gra(gra_index,:) = []; % removing from diff_gra
% removing these from dup_rxns
tot = [impr.diff_s;impr.diff_gra];
k1 = ismember(impr.dup_rxns(:,1),tot(:,1)); k2 = ismember(impr.dup_rxns(:,2),tot(:,2));
impr.dup_rxns(~(k1&k2),:) = []; 
impr.belong_dups(~(k1&k2),:) = [];