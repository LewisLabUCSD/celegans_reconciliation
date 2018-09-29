function [model,ruleReplace] = updateModelGenes(model,geneList)
% updates the list of model genes as geneList. Make sure model genes are at
% least a subset of the geneList

% written by chintan joshi on 9/14/2017 at ucsd

newX = cell(length(model.genes),1);
oldGenes = model.genes;
model.oldGenes = oldGenes;
newRxnGeneMat = zeros(size(model.rxnGeneMat));
model.genes = geneList;
oldX = regexprep(strcat({'x('},num2str([1:1:length(oldGenes)]'),')'),' ','');
for i=1:length(oldGenes)
    index = find(strcmp(geneList,oldGenes{i,1}));
    newX{i,1} = strcat('g(',num2str(index),')');
    model.rules = strrep(model.rules,oldX{i,1},newX{i,1});
    newRxnGeneMat(:,index) = model.rxnGeneMat(:,i);
end
model.rules = strrep(model.rules,'g','x');
model.rxnGeneMat = newRxnGeneMat;
ruleReplace = [strrep(oldX,'x','old') strrep(newX,'g','new')];
