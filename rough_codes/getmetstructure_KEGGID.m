function getmetstructure_KEGGID(KEGGid)
% getmetnames_KEGGID(KEGGid)
% INPUT: 
% KEGGid: KEGG compound id
% OUTPUT:
% generates a MATLAB figure of the structure of the compound.

% It still needs to be improved to
% shorten retrieval time and fails to retrieve sometimes.

data = webread(strcat('http://www.genome.jp/Fig/compound/',KEGGid,'.gif'));
data(data~=mode(data))=0;
figure
imshow(data)