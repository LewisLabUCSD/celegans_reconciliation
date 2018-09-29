function [met_list] = convert_reaction_formula2list(formula_list)

met_list = [];
for i=1:length(formula_list)
    if ~isempty(strfind(formula_list{i,1},' -> '))
        X = regexp(formula_list{i,1},' -> ','split');
    else
        X = regexp(formula_list{i,1},' <=> ','split');
    end
    X1 = regexprep(regexprep(regexp(X{1,1},'+','split'),'[0-9] ',''),' ','');
    X2 = regexprep(regexprep(regexp(X{1,2},'+','split'),'[0-9] ',''),' ','');
    X = unique([X1,X2]');
    met_list = union(met_list,X);
end
