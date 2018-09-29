function [stat,bal] = atom_balance(model,rxns,disp_flag)
% [stat,t] = atom_balance(model,rxns)
% gives the elemental imbalance in the rxns
% User-defined function used:
% 1. create_atom_numbers
% 2. create_atom_matrix
%
% INPUT:
% model: model to be checked
% rxns: list of reactions, belonging to model, to be checked
% disp_flag: print output on command window (true).
%
% OUTPUT:
% stat: balance stats (imbalance = 0, balance = 1, not checked = -1)
% bal: elemental imbalance in reactions
% all_elements = {'C';'Ca';'Cl';'Co';'Cu';'Fe';'H';'I';'K';'Mg';'Mo';'N';'Na';'O';'P';'R';'S';'Se';'Zn'}';
% rows = rxns, columns = all_elements
%
if nargin < 2
    rxns = model.rxns;
end
stat = repmat(-2,size(rxns,1),1);
if disp_flag
    fprintf('Rxn no. \tReaction name\t\t\tComment\n------------------------------------------------\n');
end
sym('n');
for i=1:length(rxns)
    rxn_index = find(strcmp(model.rxns,rxns{i,1}));
    if ~isempty(rxn_index)
        met_index = find(model.S(:,rxn_index)~=0); metp_index = find(model.S(:,rxn_index) > 0); metn_index = find(model.S(:,rxn_index) < 0);
        mets = model.mets(met_index);
        atoms_mat = create_atom_matrix(model,mets,false);
        %     printRxnFormula(model,model.rxns(i));
        if sum(cellfun(@isempty,atoms_mat),2)==size(atoms_mat,2)
            stat(i,1) = -1;
            if disp_flag
                fprintf('%d. %s: One or more metabolites missing formula.\n',i,rxns{i,1});
            end
        else
            Sj = zeros(1,size(atoms_mat,2));
            n=50;
            for j=1:length(met_index)
                atoms_mat(j,cellfun(@isempty,atoms_mat(j,:))) = {'0'};
                atoms_mat(j,strcmp(atoms_mat(j,:),'*n')) = {'0'};
                Sj = Sj + model.S(met_index(j),rxn_index)*cellfun(@eval,atoms_mat(j,:));
            end
            bal(i,:) = Sj;
            if sum(bal(i,:)) == 0 || (~isempty(metp_index) && isempty(metn_index)) || (isempty(metp_index) && ~isempty(metn_index))
                stat(i,1) = 1;
                if disp_flag
                    fprintf('%d. %s: Atoms balanced.\n',i,rxns{i,1});
                end
            else
                stat(i,1) = 0;
                if disp_flag
                    fprintf('%d. %s: Atoms not balanced.\n',i,rxns{i,1});
                end
            end
        end
    else
        stat(i,1) = -2;
        if disp_flag
            fprintf('%d. %s: Reaction not present in the model.\n',i,rxns{i,1});
        end
    end
end
stat = stat';