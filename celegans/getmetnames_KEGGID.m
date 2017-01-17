function met_names = getmetnames_KEGGID(KEGGid)
% met_names = getmetnames_KEGGID(KEGGid)
% INPUT: 
% KEGGid: KEGG compound id
% OUTPUT:
% met_names: list of all compound names retrieved from KEGG

% It still needs to be improved to
% shorten retrieval time and fails to retrieve sometimes.

data = webread(strcat('http://www.genome.jp/dbget-bin/www_bget?',KEGGid));
x = strfind(data,'</div></div></td></tr>');
y = strfind(data,'<div style="width:555px;overflow-x:auto;overflow-y:hidden"><div style="width:555px;overflow-x:auto;overflow-y:hidden">');
data = data(y+118:x(1)-1);
data = strrep(data,'<br>','');
data = regexp(data,';','split');
met_names = strtrim(data');