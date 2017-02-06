function [model1,model_unref,dup_rxns,rem_rxns,res_rxns,gra_rxns,prep_rem,account,ch] = refine_merged_model_v2(modelA,modelB,obj_rxn,remFlag)
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

% % check duplicate reactions
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
                    gra_rxns{b,1} = 'not resolved';
                end
                if ((charge1 && ~charge2) && (atom1 && ~atom2)) || ((charge1 && charge2) && (atom1 && atom2)) || ...
                        ((charge1 && ~charge2) && (atom1 && atom2)) || ((charge1 && charge2) && (atom1 && ~atom2))
                    account(b,1) = 1;
                    if strcmp(gra_rxns{b,1},'pick two')
                        k = k+1; ch = ch+1;
                        rem_rxns{k,1} = model1.rxns{indices(j),1}; % choose 2nd
                        model1.rxnGeneMat(i,:) = gra2;
                        gra_rxns{b,1} = 'resolved';
                    elseif strcmp(gra_rxns{b,1},'pick one')
                        k = k+1;
                        rem_rxns{k,1} = model1.rxns{indices(j),1}; % choose 1st
                        gra_rxns{b,1} = 'resolved';
                    elseif strcmp(gra_rxns{b,1},'pick any')
                        k = k+1;
                        rem_rxns{k,1} = model1.rxns{indices(j),1}; % choose 1st
                        gra_rxns{b,1} = 'resolving not needed';
                    else
                        [model1,recon] = reconcile_gene_associations(model1,{model1.rxns{i,1} model1.rxns{indices(j),1}});
                        if ~recon
                            p = p+1;
                            prep_rem{p,1} = model1.rxns{i,1};
                            prep_rem{p,2} = model1.rxns{indices(j),1};
                        else
                            k = k+1;
                            rem_rxns{k,1} = model1.rxns{indices(j),1};
                            gra_rxns{b,1} = 'resolved';
                        end
                    end
                elseif ((~charge1 && charge2) && (~atom1 && atom2)) || ((charge1 && charge2) && (~atom1 && atom2)) || ...
                        ((~charge1 && charge2) && (atom1 && atom2))
                    account(b,1) = 2;
                    if strcmp(gra_rxns{b,1},'pick one')
                        k = k+1; ch = ch+1;
                        rem_rxns{k,1} = model1.rxns{i,1}; % choose 1st
                        model1.rxnGeneMat(indices(j),:) = gra1;
                        gra_rxns{b,1} = 'resolved';
                    elseif strcmp(gra_rxns{b,1},'pick any')
                        k = k+1;
                        rem_rxns{k,1} = model1.rxns{i,1}; % choose 1st
                        model1.rxnGeneMat(indices(j),:) = gra1;
                        gra_rxns{b,1} = 'resolving not needed';
                    elseif strcmp(gra_rxns{b,1},'pick two')
                        k = k+1;
                        rem_rxns{k,1} = model1.rxns{i,1}; % choose 2nd
                        gra_rxns{b,1} = 'resolved';
                    else
                        [model1,recon] = reconcile_gene_associations(model1,{model1.rxns{indices(j),1} model1.rxns{i,1}});
                        if ~recon
                            p = p+1;
                            prep_rem{p,1} = model1.rxns{i,1};
                            prep_rem{p,2} = model1.rxns{indices(j),1};
                        else
                            k = k+1;
                            rem_rxns{k,1} = model1.rxns{i,1};
                            gra_rxns{b,1} = 'resolved';
                        end
                    end
                elseif (~charge1 && ~charge2) || (~atom1 && ~atom2)
                    account(b,1) = 3;
                    if strcmp(gra_rxns{b,1},'pick two')
                        k = k+1; ch = ch+1;
                        rem_rxns{k,1} = model1.rxns{indices(j),1}; % choose 2nd
                        model1.rxnGeneMat(i,:) = gra2;
                        gra_rxns{b,1} = 'resolved';
                    elseif strcmp(gra_rxns{b,1},'pick one')
                        k = k+1;
                        rem_rxns{k,1} = model1.rxns{indices(j),1}; % choose 1st
                        gra_rxns{b,1} = 'resolved';
                    elseif strcmp(gra_rxns{b,1},'pick any')
                        k = k+1;
                        rem_rxns{k,1} = model1.rxns{indices(j),1}; % choose 1st
                        gra_rxns{b,1} = 'resolving not needed';
                    else
                        [model1,recon] = reconcile_gene_associations(model1,{model1.rxns{i,1} model1.rxns{indices(j),1}});
                        if ~recon
                            p = p+1;
                            prep_rem{p,1} = model1.rxns{i,1};
                            prep_rem{p,2} = model1.rxns{indices(j),1};
                        else
                            k = k+1;
                            rem_rxns{k,1} = model1.rxns{indices(j),1};
                            gra_rxns{b,1} = 'resolved';
                        end
                    end
                elseif ((~charge1 && charge2) && (atom1 && ~atom2)) || ((charge1 && ~charge2) && (~atom1 && atom2))
                    r = r+1;
                    res_rxns{r,1} = model1.rxns{i,1};
                    res_rxns{r,2} = model1.rxns{indices(j),1};
                    account(b,1) = 3;
                    if strcmp(gra_rxns{b,1},'pick two')
                        ch = ch+1;
                        model1.rxnGeneMat(i,:) = gra2;
                        gra_rxns{b,1} = 'resolved';
                    elseif strcmp(gra_rxns{b,1},'pick one')
                        gra_rxns{b,1} = 'resolved';
                    elseif strcmp(gra_rxns{b,1},'pick any')
                        gra_rxns{b,1} = 'resolving not needed';
                    else
                        gra_rxns{b,1} = 'not resolved';
                    end
                end
            end
        end
    end
    waitbar(i/size(S,2));
end
res_rxns(ismember(res_rxns(:,1),rem_rxns) | ismember(res_rxns(:,2),rem_rxns),:) = [];
close(h);
fprintf('Looking for duplicate reactions: Finished.\n');
fprintf('\nSUMMARY:\n---------------------------------------------------------\n');
fprintf('%d unused metabolties were found and removed.\n',length(iunmets));
fprintf('%d duplicate reactions could be found.\n',length(dup_rxns));
fprintf('%d duplicate reaction could be resolved by atom/charge balancing and GRA, and can be removed.\n',length(find(account==1|account==2)));
fprintf('%d duplicate reactions resolved solely by GRA because stoichoimery of reactions was wrong.\n',length(find(account==3|strcmp(gra_rxns,'resolved'))));
fprintf('%d duplicate reactions could be resolved by atom/charge balancing but not GRA.\n',size(prep_rem,1));
fprintf('%d duplicate reactions could not be resolved using atom/charge balancing.\n',length(res_rxns));
fprintf('%d duplicate reactions for which GRA was changed.\n',ch);
fprintf('%d duplicate reactions for which GRA does not need changes.\n',length(find(strcmp(gra_rxns,'resolving not needed'))));
fprintf('%d duplicate reactions for which GRA has been resolved.\n',length(find(strcmp(gra_rxns,'resolved'))));
% fprintf('%d duplicate reactions for which neither GRA, nor stoichiometry been resolved.\n',length(find(strcmp(gra_rxns,'not resolved')))-size(prep_rem,1));
fprintf('%d duplicate reactions that still need to be resolved.\n',size(prep_rem,1)+length(res_rxns));

if remFlag
    rindex = find(ismember(model1.rxns,rem_rxns));
    model1 = update_reaction_properties(model1,rindex);
end
t = toc;