function codes = findorgcode_KEGG(org_fullname)
% codes = findorgcode_KEGG(org_fullname)
% finds and returns the all the basic info for organism and strain of
% org_fullname
%
% INPUT:
% org_fullname: full name of the organism (make sure spelling is corect)

% OUTPUT:
% codes: list of 3-4 letter codes for the organism and strains

org_table = kegg_organisms;
org_name_list = org_table(:,3);
c = 0;
for i=1:length(org_name_list)
    if ~isempty(org_name_list{i,1})
        ix = regexp(org_name_list{i,1},org_fullname,'ONCE');
        if ~isempty(ix)
            c = c+1;
            I(c,1) = i;
        end
    end
end
fprintf('%d organism(s) found.\n',c);
org_table = org_table(I,:);
codes = org_table(:,2);
fprintf('No.\tT-code\t\t\tLetter code\t\t\tName & strain\t\t\tType\n');
fprintf('------------------------------------------------------------------------------\n')
for i=1:size(org_table,1)
    fprintf('%d.\t%s\t\t\t%s\t\t\t%s\t\t\t%s\n',i,org_table{i,1},org_table{i,2},org_table{i,3},org_table{i,4});
end

