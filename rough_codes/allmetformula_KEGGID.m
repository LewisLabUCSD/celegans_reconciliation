function formula = allmetformula_KEGGID(all_keggid)
% formula = allmetformula_KEGGID(all_keggid)
% INPUT: 
% all_keggid: list of KEGG compound ids
% OUTPUT:
% formula: list of compound formulas retrieved from KEGG

% The code needs further improvement by
% 1. convert all_keggid to a unique list
% 2. use a while loop

formula = cell(length(all_keggid),1);

for i=1:length(all_keggid)
    if ~isempty(all_keggid{i,1}) && ~strncmp(all_keggid{i,1},'CC',2)
        fprintf('%d. Looking for Formula for %s: ',i,all_keggid{i,1});
        all_keggid{i,1} = regexprep(all_keggid(i,1),'[(a-z)]','');
        all_keggid{i,1} = strrep(all_keggid{i,1},'[]','');
        formula{i,1} = getmetformula_KEGGID(all_keggid{i,1}{1,1});
        fprintf('%s\n',formula{i,1});
    end
end