function geneFields = listGeneFields(model)

% % get fieldnames % %
fdNames = fieldnames(model);
fdNames(strcmp(fdNames,'S')) = [];

geneFields = []; r = 0;
for i=1:length(fdNames)
    if eval(strcat('size(model.',fdNames{i,1},',1) == size(model.genes,1)'))
        r = r + 1;
        geneFields{r,1} = fdNames{i,1};
    end
end
