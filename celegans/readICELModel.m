function [icel] = readICELModel(filename,all_hrid,all_mrid,all_mets,all_bigg,all_kegg,all_names,all_formulas)
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