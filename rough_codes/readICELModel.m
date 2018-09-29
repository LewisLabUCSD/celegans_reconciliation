function [icel,genes_found] = readICELModel(filename,all_hrid,all_mrid,all_mets,all_bigg,all_kegg,all_names,all_formulas)
% [icel] = readICELModel(filename,all_hrid,all_mrid,all_mets,all_bigg,all_kegg,all_names)
% reads the icel model and converts it in the desired format

% INPUT:
% filename: iCEL1273_1.xml, provided the file is in working directory
% all_hrid: import coloumn A of "ICEL & Common" tab of
% "E:\Dropbox\Sean-Chintan\chintan\Metabolite List.xlsx"
% all_mrid: import coloumn B of "ICEL & Common" tab of
% "E:\Dropbox\Sean-Chintan\chintan\Metabolite List.xlsx"
% all_mets: import coloumn A of "All Metabolites (repeats too)" tab of
% "E:\Dropbox\Sean-Chintan\chintan\Metabolite List.xlsx"
% all_bigg: import coloumn D of "All Metabolites (repeats too)" tab in
% "E:\Dropbox\Sean-Chintan\chintan\Metabolite List.xlsx"
% all_kegg: import coloumn B of "All Metabolites (repeats too)" tab in
% "E:\Dropbox\Sean-Chintan\chintan\Metabolite List.xlsx"
% all_names: import coloumn C of "All Metabolites (repeats too)" tab in
% "E:\Dropbox\Sean-Chintan\chintan\Metabolite List.xlsx"
% all_formulas: import coloumn E of "All Metabolites (repeats too)" tab in
% "E:\Dropbox\Sean-Chintan\chintan\Metabolite List.xlsx"

% OUTPUT:
% icel: iCEL1273 model in Cobra format

icel = readCbModel(filename);
icel.mets = strrep(icel.mets,'[Cytosol]','[c]');
icel.mets = strrep(icel.mets,'[Extracellular]','[e]');
icel.mets = strrep(icel.mets,'[Mitochondria]','[m]');
for i=1:length(icel.mets)
    if ~isempty(regexp(icel.mets{i,1},'E')) && ~isempty(regexp(icel.mets{i,1},'[c]'))
        icel.mets{i,1} = strrep(icel.mets{i,1},'[c]','[e]');
    end
    if ~isempty(regexp(icel.mets{i,1},'E')) && ~isempty(regexp(icel.mets{i,1},'[m]'))
        icel.mets{i,1} = strrep(icel.mets{i,1},'[m]','[e]');
    end
    if ~isempty(regexp(icel.mets{i,1},'M')) && ~isempty(regexp(icel.mets{i,1},'[c]'))
        icel.mets{i,1} = strrep(icel.mets{i,1},'[c]','[m]');
    end
    if ~isempty(regexp(icel.mets{i,1},'M')) && ~isempty(regexp(icel.mets{i,1},'[e]'))
        icel.mets{i,1} = strrep(icel.mets{i,1},'[e]','[m]');
    end
    if ~isempty(regexp(icel.mets{i,1},'C')) && ~isempty(regexp(icel.mets{i,1},'[m]'))
        icel.mets{i,1} = strrep(icel.mets{i,1},'[m]','[c]');
    end
    if ~isempty(regexp(icel.mets{i,1},'C')) && ~isempty(regexp(icel.mets{i,1},'[e]'))
        icel.mets{i,1} = strrep(icel.mets{i,1},'[e]','[c]');
    end
    if ~isempty(regexp(icel.mets{i,1},'CC')) && ~isempty(regexp(icel.mets{i,1},'[m]'))
        icel.mets{i,1} = strrep(icel.mets{i,1},'[m]','[c]');
    end
    if ~isempty(regexp(icel.mets{i,1},'MC')) && ~isempty(regexp(icel.mets{i,1},'[c]'))
        icel.mets{i,1} = strrep(icel.mets{i,1},'[c]','[m]');
    end
    if ~isempty(regexp(icel.mets{i,1},'EC')) && ~isempty(regexp(icel.mets{i,1},'[c]'))
        icel.mets{i,1} = strrep(icel.mets{i,1},'[c]','[e]');
    end
end
icel.mets = strrep(icel.mets,'E','C');
icel.mets = strrep(icel.mets,'M','C');
mets = icel.mets;

all_hrid = strrep(all_hrid,'[c]','');
all_hrid = strrep(all_hrid,'[e]','');
all_hrid = strrep(all_hrid,'[m]','');
all_hrid = strrep(all_hrid,'[n]','');

all_mrid = strrep(all_mrid,'[c]','');
all_mrid = strrep(all_mrid,'[e]','');
all_mrid = strrep(all_mrid,'[m]','');
all_mrid = strrep(all_mrid,'[n]','');

all_mets = strrep(all_mets,'[c]','');
all_mets = strrep(all_mets,'[e]','');
all_mets = strrep(all_mets,'[m]','');
all_mets = strrep(all_mets,'[n]','');

all_bigg = strrep(all_bigg,'[c]','');
all_bigg = strrep(all_bigg,'[e]','');
all_bigg = strrep(all_bigg,'[m]','');
all_bigg = strrep(all_bigg,'[n]','');

for i=1:length(all_kegg)
    if isempty(all_kegg{i,1})
        all_kegg{i,1} = ' ';
    else
        all_kegg{i,1} = all_kegg{i,1};
    end
end
all_kegg = strrep(all_kegg,'[c]','');
all_kegg = strrep(all_kegg,'[e]','');
all_kegg = strrep(all_kegg,'[m]','');
all_kegg = strrep(all_kegg,'[n]','');

for i=1:length(all_formulas)
    if isempty(all_formulas{i,1})
        all_formulas{i,1} = ' ';
    else
        all_formulas{i,1} = all_formulas{i,1};
    end
end

compartments = {'c';'e';'m'};
% compartments = {'m'};
acc = 0;
oldmets = [];
for i=1:length(compartments)
    all_kegg_ = strcat(all_mrid,'[',compartments{i,1},']');
    fprintf('Finding all [%s] metabolites..',compartments{i,1});
    [A,ia,ib] = intersect(all_kegg_,mets);
    fprintf('\t%d\n',length(ia));
    acc = acc + length(ia);
    for j=1:length(ia)
        comp = regexp(mets(ib(j)),'[','split');
        comp = strrep(comp{1,1}{1,2},']','');
        oldmets{ib(j),1} = strcat(all_hrid{ia(j),1},'[',comp,']');
    end
end
if ~isempty(oldmets)
    icel.metKEGGID = icel.mets; icel.mets = oldmets;
    fprintf('%d metabolites out of %d accounted for.\n',acc,length(icel.mets));
end

fprintf('Assigning BiGG IDs as metabolite identifiers.\n');
% making sure they are assigned BiGG ids
mets = icel.mets;
acc = 0;
newmets = [];
for i=1:length(compartments)
    all_mets_ = strcat(all_mets,'[',compartments{i,1},']');
    fprintf('Finding all [%s] metabolites..',compartments{i,1});
    [A,ia,ib] = intersect(all_mets_,mets);
    fprintf('\t%d\n',length(ia));
    acc = acc + length(ia);
    for j=1:length(ia)
        comp = regexp(mets(ib(j)),'[','split');
        comp = strrep(comp{1,1}{1,2},']','');
        newmets{ib(j),1} = strcat(all_bigg{ia(j),1},'[',comp,']');
        names{ib(j),1} = all_names{ia(j),1};
        keggid{ib(j),1} = all_kegg{ia(j),1};
        formula{ib(j),1} = all_formulas{ia(j),1};
    end
end
if ~isempty(newmets)
    icel.oldmets = icel.mets; icel.mets = newmets;
    fprintf('%d metabolites out of %d accounted for.\n',acc,length(icel.oldmets));
    icel.metNames = names;
    icel.metKEGGID = keggid;
    icel.metFormulas = formula;
end
model1 = icel;
% % find duplicate genes, uses WormBase and KEGG subsequently to look up information
filename_local = 'E:\Downloads\c_elegans.canonical_bioproject.current\functional_descriptions.txt';
fprintf('Looking for genes in WormBase and fixing:');
[genes_found1,genes_not_found_wb,twice_present] = getgeneinfo_WormBase(filename_local,model1.genes);
genes_found1(:,[1 4 5]) = [];
genes_found1(:,2) = strrep(genes_found1(:,2),'not known','');
empty_ind = find(cellfun(@isempty,genes_found1(:,2)));
genes_found1(empty_ind,2)=genes_found1(empty_ind,1);
genes_found1(:,2) = regexprep(genes_found1(:,2),'[a-z]','');
genes_not_found_wb(ismember(genes_not_found_wb,genes_found1(:,1))==1) = [];
genes_not_found_wb(ismember(genes_not_found_wb,genes_found1(:,2))==1) = [];
gf1 = ismember(genes_found1(:,1),twice_present(:,1)); gf2 = ismember(genes_found1(:,2),twice_present(:,2));
genes_found1(gf1 & gf2,:) = [];

fprintf('Finished\n');
fprintf('%d genes found in WormBase.\n',length(genes_found1));
fprintf('Looking for %d genes, not found in WormBase, in KEGG and fixing:',length(genes_not_found_wb));
[genes_found2,genes_not_found] = getgeneinfo_KEGG('cel',genes_not_found_wb);
genes_found2(:,[2 4]) = []; genes_found2(:,[2 1]) = genes_found2(:,[1 2]);
genes_found2(cellfun(@isempty,genes_found2(:,1))==1,:)=[];
genes_not_found(ismember(genes_not_found,genes_found2(:,1))==1) = [];
genes_not_found(ismember(genes_not_found,genes_found2(:,2))==1) = [];
model1 = update_gene_properties(model1,twice_present,1); % substitution, merge, and remove
model1 = merge_gene_properties(model1); % merge and remove
fprintf('Finished\n');
genes_found = [genes_found1;genes_found2];
model1 = update_gene_properties(model1,genes_found(ismember(genes_found(:,1),model1.genes),:),0); % substitution

% % check if duplicate genes have been updated correctly
[~,~,del_match] = check_gene_account(icel,model1,twice_present);
if sum(del_match==del_match)~=length(del_match)
    fprintf('Genes have not been updated correctly and are accounted for.\n');
end
if ~isempty(genes_not_found)
    fprintf('%d genes were not found in either database, and were left as it is.\n',size(genes_not_found,1));
end
if ~isempty(twice_present)
    fprintf('%d genes were duplicates.\n',size(twice_present,1));
end
fprintf('%d genes were finally reduced to %d genes using KEGG and WormBase.\n',length(icel.genes),length(model1.genes));

icel = model1;