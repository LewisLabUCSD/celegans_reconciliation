function [eleg] = readElegModel(filename,all_mets,all_bigg,all_kegg,all_names,all_formulas)
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