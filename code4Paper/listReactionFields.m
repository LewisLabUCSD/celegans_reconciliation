function rxnFields = listReactionFields(model)

% % get fieldnames % %
fdNames = fieldnames(model);
fdNames(strcmp(fdNames,'S')) = [];

rxnFields = []; r = 0;
for i=1:length(fdNames)
    if eval(strcat('size(model.',fdNames{i,1},',1) == size(model.rxns,1)'))
        r = r + 1;
        rxnFields{r,1} = fdNames{i,1};
    end
end
