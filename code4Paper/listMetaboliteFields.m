function metFields = listMetaboliteFields(model)

% % get fieldnames % %
fdNames = fieldnames(model);
fdNames(strcmp(fdNames,'S')) = [];

metFields = []; m = 0;
for i=1:length(fdNames)
    if eval(strcat('size(model.',fdNames{i,1},',1) == size(model.mets,1)'))
        m = m + 1;
        metFields{m,1} = fdNames{i,1};
    end
end
