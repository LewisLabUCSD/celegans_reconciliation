function [model] = rulesModel(model)

% L = length(model.genes) - length(genes);
model.genes = unique(model.genes);
x = regexprep(strcat({'x\('},num2str([1:1:length(model.genes)]'),'\)'),' ','');
model.rules = regexprep(model.grRules,strcat(model.genes,'(?!\w)'),x);
model.rules = regexprep(strrep(model.rules,'or','|'),'and','&');
model = rxnGeneMatModel(model);
% [~,ix] = unique(model.genes);
% ix = setdiff([1:1:length(model.genes)],ix);
% for i=1:length(ix)
%     dupInd = find(strcmp(model.genes,model.genes(ix(i))));
%     keepInd = dupInd(dupInd==ix(i));
%     
%     model.rxnGeneMat(:,dupInd) = [];
%     model.genes(dupInd) = [];
%     rules = 
% end
    
    