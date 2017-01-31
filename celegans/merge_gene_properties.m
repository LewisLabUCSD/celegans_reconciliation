function  [model] = merge_gene_properties(model)

% % change gene properties from gene 1 to gene 2 and remove gene 1
x_old = regexprep(strcat({'x'},num2str([1:1:length(model.genes)]')),' ','');

genes = unique(model.genes);
rem_genes = {'NA';'ND';'Unknown';'TBD'};
genes = setdiff(genes,rem_genes);

for i=1:size(genes,1)
    gindex = find(strcmp(model.genes,genes{i,1}));
    if gindex > 1
        gmain = gindex(1);
        model.rxnGeneMat(:,gmain) = any(sum(model.rxnGeneMat(:,gindex),2),3);
        gindex(gindex==gmain) = [];
        model.rxnGeneMat(:,gindex) = [];
        model.genes(gindex) = [];
        x_old(gindex) = [];
    end
end
x_new = regexprep(strcat({'x'},num2str([1:1:length(model.genes)]')),' ','');
for i=1:length(x_new)
    model.rules = regexprep(model.rules,x_old{i,1},x_new{i,1});
end