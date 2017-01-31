function [model1,model_unref,dup_rxns,rem_rxns,res_rxns,gra_rxns,prep_rem,account,ch] = refine_merged_model_v2(modelA,modelB,obj_rxn)
% [model_unref,impr,t] = refine_merged_model(model)
% refines the merged model
% user-defined functions used:
% 1. mergeTwoModels
% 2. update_metabolite_properties
% 3. getgeneinfo_WormBase
% 4. getgeneinfo_KEGG
% 5. check_gene_account
% 6. charge_balance
% 7. update_reaction_properties
% 8. ......more to come.......

% INPUT:
% modelA: 1st model
% modelB: 2nd column
% obj_rxn: objective reaction, may belong to either models
% org_code: 3-4 letter organism code

% OUTPUT:
% model1: new refined model
% impr: improvements made
% t: time taken to create the improved model (model1)

% COMMENTS: Does not take genes and stoichiometry into account! Needs to be included

tic
% % merge two models to create a preliminary model
model_unref = mergeTwoModels(modelA,modelB,obj_rxn,true);
% % initialize new_model
model1 = model_unref;

% % remove any unused metabolites and remove metabolite properties that exist
fprintf('Looking for unused metabolites:');
unused_mets = model_unref.mets(sum(any(model_unref.S,3),2)==0);
[~,iunmets] = intersect(model1.mets,unused_mets);
model1 = update_metabolite_properties(model1,iunmets);
fprintf('Finished.\n');

% % % find duplicate genes, uses WormBase and KEGG subsequently to look up information
% filename_local = 'E:\Downloads\c_elegans.canonical_bioproject.current\functional_descriptions.txt';
% fprintf('Looking for genes in WormBase and fixing...');
% [genes_found1,genes_not_found_wb,twice_present] = getgeneinfo_WormBase(filename_local,model1.genes);
% genes_found1(:,[1 4 5]) = [];
% genes_found1(:,2) = strrep(genes_found1(:,2),'not known','');
% empty_ind = find(cellfun(@isempty,genes_found1(:,2)));
% genes_found1(empty_ind,2)=genes_found1(empty_ind,1);
% genes_found1(:,2) = regexprep(genes_found1(:,2),'[a-z]','');
% genes_not_found_wb(ismember(genes_not_found_wb,genes_found1(:,1))==1) = [];
% genes_not_found_wb(ismember(genes_not_found_wb,genes_found1(:,2))==1) = [];
% gf1 = ismember(genes_found1(:,1),twice_present(:,1)); gf2 = ismember(genes_found1(:,2),twice_present(:,2));
% genes_found1(gf1 & gf2,:) = [];
% % model1 = update_gene_properties(model1,genes_found1(ismember(genes_found1(:,1),model1.genes),:),0); % substitution
% 
% fprintf('Finished\nLooking for %d genes, not found in WormBase, in KEGG and fixing...',length(genes_not_found_wb));
% [genes_found2,genes_not_found] = getgeneinfo_KEGG('cel',genes_not_found_wb);
% genes_found2(:,[2 4]) = []; genes_found2(:,[2 1]) = genes_found2(:,[1 2]);
% genes_found2(cellfun(@isempty,genes_found2(:,1))==1,:)=[];
% genes_not_found(ismember(genes_not_found,genes_found2(:,1))==1) = [];
% genes_not_found(ismember(genes_not_found,genes_found2(:,2))==1) = [];
% % model1 = update_gene_properties(model1,genes_found2(ismember(genes_found2(:,1),model1.genes),:),0); % substitution
% model1 = update_gene_properties(model1,twice_present,1); % substitution, merge, and remove
% model1 = merge_gene_properties(model1); % merge and remove
% fprintf('Finished\n');
% 
% % % check if duplicate genes have been updated correctly
% [~,~,del_match] = check_gene_account(model_unref,model1,twice_present);
% if sum(del_match==del_match)~=length(del_match)
%     fprintf('Genes have not been updated correctly and are accounted for.\n');
% end
% % impr.modelg = model1;

S = full(any(model1.S,3));
h_mets = {'h[c]';'h[m]';'h[n]';'h[e]'};
S_woh = S;
S_woh(ismember(model1.mets,h_mets),:) = [];
h = waitbar(0,'Checking reactions for duplicity...');
fprintf('Looking for duplicate reactions:');
b = 0; k = 0; r = 0; ch = 0; p = 0; account = zeros(1,1);
dup_rxns = cell(1,2); rem_rxns = cell(1,1); res_rxns = cell(1,2);
for i=1:size(S,2)
    [~,indices1] = ismember(S(:,i)',S','rows');
    if sum(S_woh(:,i))~=0
        [~,indices2] = ismember(S_woh(:,i)',S_woh','rows');
    else
        indices2 = [];
    end
    indices = unique([indices1;indices2]);
    indices(indices==i) = [];
    if ~isempty(indices)
        if ~isempty(find(strcmp(modelA.rxns,model1.rxns{i,1})))
            belong1 = modelA.description;
        else
            belong1 = modelB.description;
        end
        charge1 = charge_balance(model1,model1.rxns(i),false);
        atom1 = atom_balance(model1,model1.rxns(i),false);
        gra1 = model1.rxnGeneMat(i,:);
        for j=1:length(indices)
            if ~isempty(find(strcmp(modelA.rxns,model1.rxns{indices(j),1})))
                belong2 = modelA.description;
            else
                belong2 = modelB.description;
            end
            if ~strcmp(belong1,belong2)
                b = b+1;
                dup_rxns{b,1} = model1.rxns{i,1};
                dup_rxns{b,2} = model1.rxns{indices(j),1};
                account(b,1) = 0;
                charge2 = charge_balance(model1,model1.rxns(indices(j)),false);
                atom2 = atom_balance(model1,model1.rxns(indices(j)),false);
                gra2 = model1.rxnGeneMat(indices(j),:);
                if sum(gra1==gra2)==length(model1.genes) % pick either to keep, but we choose one
                    gra_rxns{b,1} = 'pick any'; % both gra is same
                elseif (length(union(find(gra1),find(gra2)))==length(find(gra1))) && (length(union(find(gra1),find(gra2)))~=length(find(gra2)))
                    gra_rxns{b,1} = 'pick one'; % gra2 is a subset of gra1
                elseif (length(union(find(gra1),find(gra2)))==length(find(gra2))) && (length(union(find(gra1),find(gra2)))~=length(find(gra1)))
                    gra_rxns{b,1} = 'pick two'; % gra1 is a subset of gra2
                else
                    gra_rxns{b,1} = 'pick none';
                end
                if ((charge1 && ~charge2) && (atom1 && ~atom2)) || ((charge1 && charge2) && (atom1 && atom2)) || ...
                        ((charge1 && ~charge2) && (atom1 && atom2)) || ((charge1 && charge2) && (atom1 && ~atom2))
                    account(b,1) = 1;
                    if strcmp(gra_rxns{b,1},'pick two')
                        k = k+1; ch = ch+1;
                        rem_rxns{k,1} = model1.rxns{indices(j),1}; % choose 1st
                        model1.rxnGeneMat(i,:) = gra2;
                        gra_rxns{b,1} = 'do nothing';
                    elseif strcmp(gra_rxns{b,1},'pick one') || strcmp(gra_rxns{b,1},'pick any')
                        k = k+1;
                        rem_rxns{k,1} = model1.rxns{indices(j),1}; % choose 1st
                        gra_rxns{b,1} = 'do nothing';
                    else
                        p = p+1;
                        prep_rem{p,1} = model1.rxns{indices(j),1}; % choose 1st
                    end
                elseif ((~charge1 && charge2) && (~atom1 && atom2)) || ((charge1 && charge2) && (~atom1 && atom2)) || ...
                        ((~charge1 && charge2) && (atom1 && atom2))
                    account(b,1) = 2;
                    if strcmp(gra_rxns{b,1},'pick one') || strcmp(gra_rxns{b,1},'pick any')
                        k = k+1; ch = ch+1;
                        rem_rxns{k,1} = model1.rxns{i,1}; % choose 2nd
                        model1.rxnGeneMat(indices(j),:) = gra1;
                        gra_rxns{b,1} = 'do nothing';
                    elseif strcmp(gra_rxns{b,1},'pick two')
                        k = k+1;
                        rem_rxns{k,1} = model1.rxns{i,1}; % choose 2nd
                        gra_rxns{b,1} = 'do nothing';
                    else
                        p = p+1;
                        prep_rem{p,1} = model1.rxns{i,1}; % choose 2nd
                    end
                elseif ((~charge1 && charge2) && (atom1 && ~atom2)) || ((charge1 && ~charge2) && (~atom1 && atom2)) || ...
                        (~charge1 && ~charge2) || (~atom1 && ~atom2)
                    r = r+1;
                    res_rxns{r,1} = model1.rxns{i,1};
                    res_rxns{r,2} = model1.rxns{indices(j),1};
                    account(b,1) = 3;
                end
            end
        end
    end
    waitbar(i/size(S,2));
end

close(h);
fprintf('Finished.\n');
fprintf('\nSUMMARY:\n---------------------------------------------------------\n');
fprintf('%d unused metabolties were found and removed.\n',length(iunmets));
fprintf('%d duplicate reactions could be found.\n',length(dup_rxns));
fprintf('%d duplicate reaction could be resolved by atom/charge balancing and GRA, and can be removed.\n',length(rem_rxns));
fprintf('%d duplicate reactions could be resolved by atom/charge balancing but not GRA.\n',length(prep_rem));
fprintf('%d duplicate reactions could not be resolved using atom/charge balancing.\n',length(res_rxns));
fprintf('%d duplicate reactions for which GRA was changed.\n',ch);
fprintf('%d duplicate reactions for which GRA has been resolved.\n',length(find(strcmp(gra_rxns,'do nothing'))));
fprintf('%d duplicate reactions for which GRA has not been resolved.\n',length(find(strcmp(gra_rxns,'pick none'))));
fprintf('%d duplicate reactions that still need to be removed.\n',length(prep_rem)+length(res_rxns));

t = toc;