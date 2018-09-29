function [eleg,genes_found,genes_not_found] = readElegModel(filename,all_mets,all_bigg,all_kegg,all_names,all_formulas)
% [eleg] = readElegModel(filename,all_mets,all_bigg,all_kegg,all_names)
% reads the Celegans (axenic or ecoli) model and converts it in the necessary format

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
eleg = grRulesModel(eleg); eleg = rxnGeneMatModel(eleg);
% find(strcmp(eleg.genes,'mthf-1'))
% find(strcmp(eleg.genes,'C06A8.1'))
% find(strcmp(eleg.genes,'C06A8.1a'))

% remove unused metabolites
unusedMets = find(~any(eleg.S,2));
metFields = listMetaboliteFields(eleg);
for i=1:length(metFields)
    eval(strcat('eleg.',metFields{i,1},'(unusedMets) = [];'));
end
eleg.S(unusedMets,:) = [];
orig_eleg = eleg;

eleg.mets = strrep(eleg.mets,'[Cytosol]','[c]');
eleg.mets = strrep(eleg.mets,'[Extraorganism]','[e]');
eleg.mets = strrep(eleg.mets,'[Mitochondrion]','[m]');
eleg.mets = strrep(eleg.mets,'[Nucleus]','[n]');

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
    [~,ia,ib] = intersect(all_mets_,mets);
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
model_unref = eleg;
model1 = eleg;

% % find duplicate genes, uses WormBase and KEGG subsequently to look up information
currentFolder = pwd;
filename_local = strcat(currentFolder,'\c_elegans.canonical_bioproject.current\functional_descriptions.txt');
fprintf('Looking for genes in WormBase and fixing:');
[genes_found1,genes_not_found_wb] = getgeneinfo_WormBase(filename_local,model1.genes);
genes_found1(:,[4 5]) = [];
% find those for which generic name was not found
genes_found1(:,2) = strrep(genes_found1(:,2),'not known','');
empty_ind = find(cellfun(@isempty,genes_found1(:,2)));
genes_found1(empty_ind,2)=genes_found1(empty_ind,1);
% genes_found1(:,2) = regexprep(genes_found1(:,2),'[a-z]','');
genes_not_found_wb(ismember(genes_not_found_wb,genes_found1(:,1))==1) = [];
genes_not_found_wb(ismember(genes_not_found_wb,genes_found1(:,2))==1) = [];
genes_found1 = [genes_found1(:,2) genes_found1(:,1);genes_found1(:,3) genes_found1(:,1)];
% gf1 = ismember(genes_found1(:,1),twice_present(:,1)); 
% gf2 = ismember(genes_found1(:,2),twice_present(:,2));
% genes_found1(gf1 & gf2,:) = [];

fprintf('Finished\n');
fprintf('%d genes found in WormBase.\n',length(genes_found1));
fprintf('Looking for %d genes, not found in WormBase, in KEGG and fixing:\n',length(genes_not_found_wb));
[genes_found2,genes_not_found] = getgeneinfo_KEGG('cel',genes_not_found_wb);
genes_found2(:,[2 4]) = []; genes_found2(:,[2 1]) = genes_found2(:,[1 2]);
genes_found2(cellfun(@isempty,genes_found2(:,1))==1,:)=[];
genes_not_found(ismember(genes_not_found,genes_found2(:,1))==1) = [];
genes_not_found(ismember(genes_not_found,genes_found2(:,2))==1) = [];
% model1 = update_gene_properties(model1,twice_present,1); % substitution, merge, and remove
fprintf('Finished\n');
genes_found = [genes_found1;genes_found2];
model1 = update_gene_properties(model1,genes_found(ismember(genes_found(:,1),model1.genes),:),0); % final substitution
if length(model1.genes) ~= length(unique(model1.genes))
    model1 = rulesModel(model1); % merge and remove
end

% % check if duplicate genes have been updated correctly
% [~,~,del_match] = check_gene_account(eleg,model1,twice_present);
% if sum(del_match==del_match)~=length(del_match)
%     fprintf('Genes have not been updated correctly and are accounted for.\n');
% end
if ~isempty(genes_not_found)
    fprintf('%d genes were not found in either database, and were left as it is.\n',size(genes_not_found,1));
else
    fprintf('All genes were found in at least one of the databases.\n');
end
% if ~isempty(twice_present)
%     fprintf('%d genes were duplicates.\n',size(twice_present,1)+L);
% end
fprintf('%d genes were finally reduced to %d genes using KEGG and WormBase.\n',length(eleg.genes),length(model1.genes));

eleg = model1;

% resolve non-unique metabolites
indexm = find(cellfun(@isempty,eleg.mets));
indextm = zeros(length(indexm),1);
for i=1:length(indexm)
    indextm(i) = setdiff(find(strcmp(orig_eleg.mets,orig_eleg.mets{indexm(i),1})),indexm);
    eleg.S(indextm(i),:) = eleg.S(indextm(i),:) + eleg.S(indexm(i),:);
end
fprintf('List of non-unique metabolites that were resolved:\n');
fprintf('%s\n',strjoin(orig_eleg.mets(indexm),', '));
eleg.mets(indexm) = []; eleg.metNames(indexm) = []; eleg.metNotes(indexm) = [];
eleg.metFormulas(indexm) = []; eleg.metCharges(indexm) = []; eleg.b(indexm) = [];
eleg.metKEGGID(indexm) = []; eleg.oldmets(indexm) = [];
eleg.csense(indexm) = []; eleg.S(indexm,:) = [];
fprintf('Non-unique metabolites resolved.\n');
% find(strcmp(eleg.genes,'mthf-1'))
% find(strcmp(eleg.genes,'C06A8.1'))
% find(strcmp(eleg.genes,'C06A8.1a'))
end

function [model,grRules] = grRulesModel(model)
% generate grRules field
grRules = model.rules;
g = regexprep(strcat({'x\('},num2str([1:1:length(model.genes)]'),'\)'),' ','');
grRules = regexprep(grRules,g,model.genes);
grRules = regexprep(regexprep(grRules,'\|','or'),'\&','and');
model.grRules = grRules;
end

function [model,rxnGeneMat] = rxnGeneMatModel(model)
% generate rxnGeneMat field
rxnGeneMat = zeros(length(model.rxns),length(model.genes));
gnum = regexprep(strcat({'x('},num2str([1:1:length(model.genes)]'),')'),' ','');
for i=1:length(gnum)
    indices = find(~cellfun(@isempty,strfind(model.rules,gnum{i})));
    if ~isempty(indices)
        rxnGeneMat(indices,i) = 1;
    end
end
model.rxnGeneMat = rxnGeneMat;
end