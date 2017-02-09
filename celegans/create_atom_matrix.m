function atoms_mat = create_atom_matrix(model,mets,disp_flag)
% atoms_mat = create_atom_matrix(model,mets)
% creates a cell matrix which can be evaluated using eval to determine atom balancing
% User-defined functions used:
% 1. create_atom_numbers
%
% INPUT:
% model: model to be used which contains the info
% mets: list of metabolites
% disp_flag: true, print summary; false, don't print summary
%
% OUTPUT:
% atoms_mat: a cell matrix with atom numbers listed as string
% Columns = mets, Rows = elements
%
if nargin < 2
    mets = model.mets;
    metFormulas = model.metFormulas;
else
    for i=1:length(mets)
        metFormulas{i,1} = model.metFormulas{ismember(model.mets,mets{i,1}),1};
    end
end
if nargin < 3
    disp_flag = false;
end

elements = {'C';'Ca';'Cl';'Co';'Cu';'Fe';'H';'I';'K';'Mg';'Mo';'N';'Na';'O';'P';'R';'S';'Se';'Zn'};
atoms_mat = cell(size(mets,1),size(elements,1));
for i=1:length(mets)
    if isempty(regexp(metFormulas{i,1},' |Generic'))
        x = regexp(metFormulas{i,1},'\)[a-z]|\)[1-9]');
        if isempty(x)
            for j=1:length(elements)
                atoms_mat{i,j} = create_atom_numbers(metFormulas{i,1},elements{j,1},false);
            end
        else
            multiplier = metFormulas{i,1}(x+1);
            metFormulas_ = regexp(metFormulas{i,1},strcat('\(|\)',multiplier),'split')';
            metFormulas_(cellfun(@isempty,metFormulas_)) = [];
            for j=1:length(elements)
                for k=1:length(metFormulas_)
                    atoms_mat_{k,j} = create_atom_numbers(metFormulas_{k,1},elements{j,1},false);
                    if ~isempty(regexp(metFormulas{i,1},strcat('\(',metFormulas_{k,1},'\)',multiplier)))
                        if ~isempty(atoms_mat_{k,j})
                            atoms_mat_{k,j} = strcat(atoms_mat_{k,j},'*',multiplier);
                        end
                    end
                    if k==1
                        atoms_mat{i,j} = strcat(atoms_mat{i,j},atoms_mat_{k,j});
                    else
                        if ~isempty(atoms_mat_{k,j})
                            atoms_mat{i,j} = strcat(atoms_mat{i,j},'+',atoms_mat_{k,j});
                        end
                    end
                end
            end
            atoms_mat(i,:) = strrep(atoms_mat(i,:),strcat('+*',multiplier),'');
        end
    end
    if sum(cellfun(@isempty,atoms_mat(i,:))) == length(elements)
        met_contain(i,1) = 0;
    else
        met_contain(i,1) = 1;
    end
end
if disp_flag
    fprintf('%d metabolites do not have formulas.\n',length(find(met_contain==0)));
    fprintf('%d metabolites have formulas.\n',length(find(met_contain==1)));
end