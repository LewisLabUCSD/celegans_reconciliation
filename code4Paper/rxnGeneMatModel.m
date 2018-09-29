function [model,rxnGeneMat] = rxnGeneMatModel(model)
% creates a rxnGeneMat field for model

% written by chinta joshi on 9/12/2017 at ucsd

rxnGeneMat = zeros(length(model.rxns),length(model.genes));
g = regexprep(strcat({'x('},num2str([1:1:length(model.genes)]'),')'),' ','');
for i=1:length(g)
    indices = find(~cellfun(@isempty,strfind(model.rules,g{i})));
    if ~isempty(indices)
        rxnGeneMat(indices,i) = 1;
    end
end
model.rxnGeneMat = rxnGeneMat;