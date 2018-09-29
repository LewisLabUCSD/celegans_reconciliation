function [model1,model_unref,duplicateBelongs,allIndices,remainDuplicates] = refine_merged_model_v3(modelA,modelB,obj_rxn,remFlag)
% [model_unref,impr,t] = refine_merged_model(model)
% refines the merged model
% user-defined functions used:
% 1. mergeTwoModels
% 2. update_metabolite_properties
% 3. check_gene_account
% 4. charge_balance
% 5. update_reaction_properties
% 8. ......more to come.......

% INPUT:
% modelA: 1st model
% modelB: 2nd column
% obj_rxn: objective reaction, may belong to either models
% remFlag: remove reactions

% OUTPUT:
% model1: new refined model
% impr: improvements made
% t: time taken to create the improved model (model1)

tic
if nargin < 4
    remFlag = true;
end
% % merge two models to create a preliminary model
model_unref = mergeModelInitial(modelA,modelB,obj_rxn);
% % initialize new_model
model1 = model_unref;

% % remove any duplicate metabolites and remove metabolite properties that
% exist
fprintf('Looking for duplicate metabolites:');
[~,idupmets] = unique(model_unref.mets);
idupmets = setdiff([1:1:length(model_unref.mets)],idupmets);
if ~isempty(idupmets)
    for i=1:length(idupmets)
        m = setdiff(find(strcmp(model_unref.mets,model_unref.mets(idupmets(i)))),idupmets(i));
        model_unref.S(m,:) = model_unref.S(m,:) + model_unref.S(idupmets(i),:);
        model_unref.S(idupmets(i),:) = zeros(1,size(model_unref.S,2));
    end
end
fprintf('Finished.\n');

% % remove any unused metabolites and remove metabolite properties that exist
fprintf('Looking for unused metabolites:');
unused_mets = model_unref.mets(sum(any(model_unref.S,3),2)==0);
[~,iunmets] = intersect(model1.mets,unused_mets);
if ~isempty(iunmets)
    model1 = update_metabolite_properties(model1,iunmets);
end
fprintf('Finished.\n');

al = length(modelA.rxns);

% % check duplicate reactions
S = full(any(model1.S,3));
h_mets = {'h[c]';'h[m]';'h[n]';'h[e]'};
S_woh = S;
S_woh(ismember(model1.mets,h_mets),:) = []; % remove protons to account for duplicate differently balanced reactions
h = waitbar(0,'Checking reactions for duplicity...');
fprintf('Looking for duplicate reactions:\n');
b = 0; k = 0; r = 0; ch = 0; p = 0; account = zeros(1,1); %allIndices = zeros(al,1); % initialize arrays and counters
dup_rxns = []; rem_rxns = []; res_rxns = []; prep_rem = []; mult_dups = [];
for i=1:al
%     model1.rxns{i}
    [~,indices1] = ismember(S(:,i)',S(:,al+1:end)','rows'); % check if any reactions in model1 are similar to reactions in model 2
    if sum(any(S_woh(:,i)))~=0
        [~,indices2] = ismember(S_woh(:,i)',S_woh(:,al+1:end)','rows');
    else
        indices2 = [];
    end
    indices = unique([indices1;indices2]); % combine both indices
    if ~isempty(indices)
        allIndices{i,1} = i;
        allIndices{i,2} = al + indices;
    end
    if ~isempty(indices)
        indices = indices + al;
        if ~isempty(find(strcmp(modelA.rxns,model1.rxns{i,1})))
            belong1 = modelA.description;
        else
            belong1 = modelB.description;
        end
        charge1 = charge_balance(model1,model1.rxns(i),false);
        atom1 = atom_balance(model1,model1.rxns(i),false);
        gra1 = model1.rxnGeneMat(i,:); grRule1 = model1.grRules(i);
        for j=1:length(indices)
            if ~isempty(find(strcmp(modelA.rxns,model1.rxns{indices(j),1})))
                belong2 = modelA.description;
            else
                belong2 = modelB.description;
            end
            if ~strcmp(belong1,belong2) %|| isempty(find(strcmp(rem_rxns,model1.rxns{indices(j),1})))
                b = b+1;
                dup_rxns{b,1} = model1.rxns{i,1};
                dup_rxns{b,2} = model1.rxns{indices(j),1};
                account(b,1) = 0;
                charge2 = charge_balance(model1,model1.rxns(indices(j)),false);
                atom2 = atom_balance(model1,model1.rxns(indices(j)),false);
                gra2 = model1.rxnGeneMat(indices(j),:); grRule2 = model1.grRules(indices(j));
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
                        model1.rxnGeneMat(i,:) = gra2; model1.grRules(i) = grRule2;
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
%                         {model1.rxns{i,1} model1.rxns{indices(j),1}}
                        [model1,recon] = reconcile_gene_associations(model1,{model1.rxns{i,1} model1.rxns{indices(j),1}},false);
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
%                         {model1.rxns{indices(j),1} model1.rxns{i,1}}
                        [model1,recon] = reconcile_gene_associations(model1,{model1.rxns{indices(j),1} model1.rxns{i,1}},false);
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
                        model1.rxnGeneMat(i,:) = gra2; model1.grRules(i) = grRule2;
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
%                         {model1.rxns{i,1} model1.rxns{indices(j),1}}
                        [model1,recon] = reconcile_gene_associations(model1,{model1.rxns{i,1} model1.rxns{indices(j),1}},false);
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
                        model1.rxnGeneMat(i,:) = gra2; model1.grRules(i) = grRule2;
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
    waitbar(i/al);
end

fprintf('Finished.\n');
close(h);
fprintf('Checking if revised charge and formula information from BiGG resolves reactions:');
if ~isempty(res_rxns)
    [unatom1,uncharge1,model1] = check_charge_formula_improvement(model1,res_rxns(:,1));
    [unatom2,uncharge2,model1] = check_charge_formula_improvement(model1,res_rxns(:,2));
    rem_rxns = [rem_rxns; res_rxns(uncharge1 & unatom1',1)];
    rem_rxns = [rem_rxns; res_rxns(uncharge2 & unatom2',2)];
    res_rxns((uncharge1 & unatom1') | (uncharge2 & unatom2'),:) = [];
    fprintf('Finished.\n');
    mult_dups = res_rxns(ismember(res_rxns(:,1),rem_rxns) | ismember(res_rxns(:,2),rem_rxns),:);
    res_rxns(ismember(res_rxns(:,1),rem_rxns) | ismember(res_rxns(:,2),rem_rxns),:) = [];
end

fprintf('\nSUMMARY:\n---------------------------------------------------------\n');
fprintf('%d unused metabolties were found and removed.\n',length(iunmets));
if ~isempty(dup_rxns)
    fprintf('%d duplicate reactions could be found.\n',length(dup_rxns));
end
fprintf('%d duplicate reaction could be resolved by atom/charge balancing and GRA, and can be removed.\n',length(find(account==1|account==2)));
fprintf('%d duplicate reactions resolved solely by GRA because stoichoimery of reactions was wrong.\n',length(find(account==3|strcmp(gra_rxns,'resolved'))));
if ~isempty(prep_rem)
    fprintf('%d duplicate reactions could be resolved by atom/charge balancing but not GRA.\n',size(prep_rem,1));
end
if ~isempty(res_rxns)
    fprintf('%d duplicate reactions could not be resolved using atom/charge balancing.\n',size(res_rxns,1));
else
    fprintf('All duplicate reactions could be resolved using atom/charge balancing.\n');
end
fprintf('%d duplicate reactions for which GRA was changed.\n',ch);
fprintf('%d duplicate reactions for which GRA does not need changes.\n',length(find(strcmp(gra_rxns,'resolving not needed'))));
fprintf('%d duplicate reactions for which GRA has been resolved.\n',length(find(strcmp(gra_rxns,'resolved'))));

if ~isempty(res_rxns)
    fprintf('%d duplicate reactions that still need to be resolved.\n',size(prep_rem,1)+size(res_rxns,1));
else
    fprintf('%d duplicate reactions that still need to be resolved.\n',size(prep_rem,1));
end
if ~isempty(mult_dups)
    fprintf('%d reactions that have multiple duplicate pairs and were resolved.\n',size(mult_dups,1));
else
    fprintf('No multiple duplicates were found.\n');
end
model1.rules = regexprep(model1.rules,'and','\&');
model1.rules = regexprep(model1.rules,'or','\|');

if remFlag
    fprintf('Removing %d reactions\n',length(rem_rxns));
    rindex = unique(find(ismember(model1.rxns,rem_rxns)));
    model1 = update_reaction_properties(model1,rindex);
%     model1 = changeObjective(model1,obj_rxn); % turn this on if you want to set the objective function
    if sum(strcmp(modelA.rxns,obj_rxn))
        modelB.description
        fprintf('Removing the objective function of %s: %s\n',modelB.description,modelB.rxns{modelB.c==1});
        model1 = update_reaction_properties(model1,find(ismember(model1.rxns,modelB.rxns{modelB.c==1})));
    else
        fprintf('Removing the objective function of %s: %s\n',modelA.description,modelA.rxns{modelA.c==1});
        model1 = update_reaction_properties(model1,find(ismember(model1.rxns,modelA.rxns{modelA.c==1})));
    end
    fprintf('Finished.\n');
end

[~,otherDuplicates] = unique(full(model1.S'),'rows');
otherDuplicates = setdiff([1:1:length(model1.rxns)]',otherDuplicates);
duplicateBelongs = cell(length(otherDuplicates),1);
for i=1:length(otherDuplicates)
    fprintf('Original reaction: %s\n',model1.rxns{otherDuplicates(i)});
    printRxnFormula(model1,'rxnAbbrList',model1.rxns{otherDuplicates(i)},'gprFlag',true);
    foundDuplicates = find(ismember(full(model1.S)',full(model1.S(:,otherDuplicates(i))'),'rows'));
    foundDuplicates(foundDuplicates == otherDuplicates(i)) = [];
    originalDuplicates(i) = foundDuplicates;
    fprintf('Duplicate reactions: \n');
    printRxnFormula(model1,'rxnAbbrList',model1.rxns(foundDuplicates),'gprFlag',true);
    if otherDuplicates(i) > length(modelA.rxns)
        if ~isempty(find(foundDuplicates > length(modelA.rxns)))
            duplicateBelongs{i} = modelB.description;
        else
            duplicateBelongs{i} = 'different models';
        end
    else
        if ~isempty(find(foundDuplicates <= length(modelA.rxns)))
            duplicateBelongs{i} = modelA.description;
        else
            duplicateBelongs{i} = 'different model';
        end
    end
end
remainDuplicates = [originalDuplicates' otherDuplicates];
t = toc;