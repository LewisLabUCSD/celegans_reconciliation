function met_formula = getmetformula_KEGGID(KEGGid)
% met_formula = getmetformula_KEGGID(KEGGid)
% INPUT: 
% KEGGid: KEGG compound id
% OUTPUT:
% met_formula: compound formula retrieved from KEGG

% It is used by "allmetformula_KEGGID.m" but still needs to be improved to
% shorten retrieval and fails to retrieve sometimes

data = webread(strcat('http://www.genome.jp/dbget-bin/www_bget?',KEGGid));
% data = webread(strcat('http://www.kegg.jp/entry/',KEGGid));
if isempty(strfind(data,'No such data'))
    y = strfind(data,'<td class="td20" style="border-color:#000; border-width: 1px 1px 0px 1px; border-style: solid"><div style="width:555px;overflow-x:auto;overflow-y:hidden">');
    data = data(y(1)+154:end);
    y = strfind(data,'<br>');
    data = data(1:y(1)-1);
    met_formula = strtrim(data);
else
    met_formula = [];
    fprintf('No such KEGG data.\n');
end