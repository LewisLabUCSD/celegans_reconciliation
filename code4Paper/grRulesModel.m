function [model,grRules] = grRulesModel(model)
% creates a grRules field for model

% written by chinta joshi on 9/12/2017 at ucsd

grRules = model.rules;
g = regexprep(strcat({'x\('},num2str([1:1:length(model.genes)]'),'\)'),' ','');
grRules = regexprep(grRules,g,model.genes);
grRules = regexprep(regexprep(grRules,'\|','or'),'\&','and');
model.grRules = grRules;