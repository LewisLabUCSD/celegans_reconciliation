function [org_table] = kegg_organisms
% [org_table] = kegg_organisms
% generates a table with rows in order of T-code, Letter code, Name &
% strain, and Type.

% OUTPUT:
% org_table: actual table as a cell array of strings

data = webread('http://rest.kegg.jp/list/organism');
ix = regexp(data,'T[0-9][0-9][0-9][0-9][0-9]')';
org = cell(length(ix),1);
org_table = cell(length(ix),4);
for i=1:length(ix)
    if i==length(ix)
        org{i,1} = data(ix(i):end);
    else
        org{i,1} = data(ix(i):ix(i+1)-1);
    end
    j = regexp(org{i,1},'\t');
    for k = 1:length(j)+1
        if k==1
            org_table{i,k} = strtrim(org{i,1}(1:j(k)));
        elseif k==length(j)+1
            org_table{i,k} = strtrim(org{i,1}(j(k-1)+1:end));
        else
            org_table{i,k} = strtrim(org{i,1}(j(k-1)+1:j(k)));
        end
    end
end
