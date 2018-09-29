function met_reactions = getmetreactions_KEGGID(KEGGid)
% met_reactions = getmetreactions_KEGGID(KEGGid)
% INPUT: 
% KEGGid: KEGG compound id
% OUTPUT:
% met_reactions: list of all KEGG reaction ids retrieved from KEGG
% compound, may use the same enzyme or not.

% It still needs to be improved to
% shorten retrieval time and fails to retrieve sometimes.

data = webread(strcat('http://www.genome.jp/dbget-bin/www_bget?',KEGGid));
y = strfind(data,'<tr><th class="th21" align="left" valign="top" style="border-color:#000; border-width: 1px 0px 0px 1px; border-style: solid"><nobr>Reaction</nobr></th>');
data = data(y(1):end);
y = strfind(data,'<div style="width:555px;overflow-x:auto;overflow-y:hidden">');
data = data(y(1):end);
data = strrep(data,'<div style="width:555px;overflow-x:auto;overflow-y:hidden">','');
y = strfind(data,'</div></td></tr>');
data = strtrim(data(1:y(1)-1));
data = strtrim(regexp(data,'</a> ','split'))';

for i=1:length(data)
    k = data{i,1};
    met_reactions{i,1} = k(strfind(k,'rn:')+3:strfind(k,'">')-1);
end