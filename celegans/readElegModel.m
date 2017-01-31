function [eleg,genes_found,twice_present,genes_not_found] = readElegModel(filename,all_mets,all_bigg,all_kegg,all_names,all_formulas)
% [eleg] = readElegModel(filename,all_mets,all_bigg,all_kegg,all_names)
% reads the icel model and converts it in the necessary format

% INPUT:
% filename: iCEL1273_1.xml, provided the file is in working directory
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
% eleg: Kaleta model in COBRA format

eleg = readCbModel(filename);
eleg.mets = strrep(eleg.mets,'_c[Cytosol]','[c]');
eleg.mets = strrep(eleg.mets,'_e[Cytosol]','[e]');
eleg.mets = strrep(eleg.mets,'_m[Cytosol]','[m]');
eleg.mets = strrep(eleg.mets,'_n[Cytosol]','[n]');
eleg.mets = strrep(eleg.mets,'_c[Extraorganism]','[c]');
eleg.mets = strrep(eleg.mets,'_e[Extraorganism]','[e]');
eleg.mets = strrep(eleg.mets,'_m[Extraorganism]','[m]');
eleg.mets = strrep(eleg.mets,'_n[Extraorganism]','[n]');
eleg.mets = strrep(eleg.mets,'_c[Mitochondrion]','[c]');
eleg.mets = strrep(eleg.mets,'_e[Mitochondrion]','[e]');
eleg.mets = strrep(eleg.mets,'_m[Mitochondrion]','[m]');
eleg.mets = strrep(eleg.mets,'_n[Mitochondrion]','[n]');
eleg.mets = strrep(eleg.mets,'_c[Nucleus]','[c]');
eleg.mets = strrep(eleg.mets,'_e[Nucleus]','[e]');
eleg.mets = strrep(eleg.mets,'_m[Nucleus]','[m]');
eleg.mets = strrep(eleg.mets,'_n[Nucleus]','[n]');
eleg.mets = strrep(eleg.mets,'_i[Cytosol]','_i[c]');

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

mets = eleg.mets;
compartments = {'c';'e';'m';'n'};
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
    eleg.oldmets = eleg.mets; eleg.mets = newmets;
    fprintf('%d metabolites out of %d accounted for.\n',acc,length(eleg.oldmets));
    eleg.metNames = names;
    eleg.metKEGGID = keggid;
    eleg.metFormulas = formula;
end
printRxnFormula(eleg,{'MDH_m';'MDH_c'},true,true,false,1,true,true);
% find(strcmp(eleg.genes,'F46E10.10'))
% find(strcmp(eleg.genes,'Y49E10.11'))
model1 = eleg;
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
model1 = update_gene_properties(model1,genes_found(ismember(genes_found(:,1),model1.genes),:),0); % final substitution

% % check if duplicate genes have been updated correctly
[~,~,del_match] = check_gene_account(eleg,model1,twice_present);
if sum(del_match==del_match)~=length(del_match)
    fprintf('Genes have not been updated correctly and are accounted for.\n');
end
if ~isempty(genes_not_found)
    fprintf('%d genes were not found in either database, and were left as it is.\n',size(genes_not_found,1));
else
    fprintf('All genes were found in at least one of the databases.\n');
end
if ~isempty(twice_present)
    fprintf('%d genes were duplicates.\n',size(twice_present,1));
end
fprintf('%d genes were finally reduced to %d genes using KEGG and WormBase.\n',length(eleg.genes),length(model1.genes));

eleg = model1;