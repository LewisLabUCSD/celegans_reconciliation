function met_enzymes = getmetenzymes_KEGGID(KEGGid)
% met_enzymes = getmetenzymes_KEGGID(KEGGid)
% INPUT: 
% KEGGid: KEGG compound id
% OUTPUT:
% met_enzymes: list of all EC numbers which use this retrieved from KEGG

% It still needs to be improved to
% shorten retrieval time and fails to retrieve sometimes.

data = webread(strcat('http://www.genome.jp/dbget-bin/www_bget?',KEGGid));
y = strfind(data,'<nobr>Enzyme</nobr></th>');
data = data(y(1):end);
y = strfind(data,'<div style="width:555px;overflow-x:auto;overflow-y:hidden">');
data = data(y(1):end);
data = strrep(data,'<div style="width:555px;overflow-x:auto;overflow-y:hidden">','');
y = strfind(data,'</div></td></tr>');
data = strtrim(data(1:y(1)-1));
data = strtrim(regexp(data,'&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;','split'))';
for i=1:length(data)
    k = data{i,1};
    if ~isempty(strfind(k,'ec:'))
        met_enzymes{i,1} = k(strfind(k,'ec:')+3:strfind(k,'">')-1);
    else
        met_enzymes{i,1} = k;
    end
end