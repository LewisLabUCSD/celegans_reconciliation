function [model_unref,impr,t] = refine_merged_model(modelA,modelB,obj_rxn,org_code)
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
% merge two models to create a preliminary model
model_unref = mergeTwoModels(modelA,modelB,obj_rxn);
% initialize new_model
model1 = model_unref;

% % remove any unused metabolites and remove metabolite properties that exist
fprintf('Looking for unused metabolites...');
unused_mets = model_unref.mets(sum(any(model_unref.S,3),2)==0);
[~,iunmets] = intersect(model1.mets,unused_mets);
model1 = update_metabolite_properties(model1,iunmets);
fprintf('Finished.\n');

% % find duplicate genes, uses WormBase and KEGG subsequently to look up information
% filename_local = 'E:\Downloads\c_elegans.canonical_bioproject.current\functional_descriptions.txt';
fprintf('Looking for genes in WormBase...');
[~,genes_not_found_wb,twice_present] = getgeneinfo_WormBase('online',model1.genes);
fprintf('Finished\nLooking for genes, not found in WormBase, in KEGG...');
[genes_found,genes_not_found] = getgeneinfo_KEGG(org_code,genes_not_found_wb);
genes_found(:,[2 4]) = []; genes_found(:,[2 1]) = genes_found(:,[1 2]);
g1 = ismember(genes_found(:,1),model1.genes);
g2 = ismember(genes_found(:,2),model1.genes);
twice_present = [twice_present;genes_found(g1 & g2,:)];
fprintf('Finished\nCombining duplicate genes...');
h = waitbar(0,'Combining duplicate genes...');
steps = size(twice_present,1);
grRules = model1.grRules;
rules = model1.rules;
x = strcat({'x'},num2str([1:1:length(model1.genes)]'));
x = strrep(x,' ','');
gindex1 = zeros(size(twice_present,1),1);
gindex2 = zeros(size(twice_present,1),1);
% % first fixing refers to changing gene properties based on current rules
for i=1:size(twice_present,1)
    % % first fixing reaction-gene matrix
    gindex1(i,1) = find(strcmp(model1.genes,twice_present{i,1}));
    gindex2(i,1) = find(strcmp(model1.genes,twice_present{i,2}));
    gene_vec_ori = model1.rxnGeneMat(:,strcmp(model1.genes,twice_present{i,2}));
    gene_vec_ori = gene_vec_ori + model1.rxnGeneMat(:,strcmp(model1.genes,twice_present{i,1}));
    model1.rxnGeneMat(:,strcmp(model1.genes,twice_present{i,2})) = gene_vec_ori;
    % % first fixing grRules
    grRules = strrep(grRules,twice_present{i,1},twice_present{i,2});
    % % first fixing rules
    rules = strrep(rules,x(gindex1(i,1)),x(gindex2(i,1)));
    waitbar(i/steps);
end
% % remove duplicate gene indices from model gene list & reaction-gene matrix
model1.genes(gindex1) = []; x(gindex1) = [];
model1.rxnGeneMat(:,gindex1) = [];
close(h);
% % second fixing begins
fprintf('Finished\nUpdating rules...');
h = waitbar(0,'Updating rules...');
steps = size(twice_present,1);
rules = model1.rules;
y = strcat({'x'},num2str([1:1:length(model1.genes)]'));
y = strrep(x,' ','');
% % second fixing refers to changing only the gene rule property based on updated rules
for i=1:size(y,1)
    rules = strrep(rules,x{i,1},y{i,1});
    waitbar(i/steps);
end
close(h);
fprintf('Finished\n');

% % check if duplicate genes have been updated corectly
[~,~,del_match] = check_gene_account(model_unref,model1,twice_present);
if sum(del_match==del_match)==length(del_match)
    fprintf('Genes have updated correctly and are accounted for.\n');
else
    fprintf('Genes have not been updated correctly and are accounted for.\n');
end
fprintf('Looking for duplicate reactions...');
% % find duplicate reactions
r = 0; % % reversibility mismatch counter
kg = 0; % % net reaction removal involvement counter
k = 0; % % duplicate reaction but different gene reaction counter
s = 0; % % stoichiometry difference counter
d = 0; % % net duplicate reaction counter
c = 1; % % query reaction counter (WHILE loop)
h = waitbar(c,'Checking reactions for duplicity and reversibility....');
while c~=length(model1.rxns)
    steps = length(model1.rxns);
    S = model1.S; % % to account for h differences..e.g. A + B -> C (vs.) A + B -> C + h
    S(ismember(model1.mets,{'h[c]';'h[m]';'h[n]';'h[e]'}),:) = [];
    %         fprintf('Looking for a duplicate copy of %s....\n',model1.rxns{c,1});
    g = 0; % % mismatch removal involvement counter
    iremrxns = [];
    for i=c+1:length(model1.rxns)
        if (sum(any(S(:,i),2)==any(S(:,c),2))==size(S,1))
%             fprintf('Another copy of %s exists at %d: %s\n',model1.rxns{c,1},i,model1.rxns{i,1});
            X(1,1) = printRxnFormula(model1,model1.rxns{c,1},false);
%             fprintf('Original: %s: %s\n',model1.rxns{c,1},X{1,1});
            X(2,1) = printRxnFormula(model1,model1.rxns{i,1},false);
%             fprintf('Duplicate: %s: %s\n',model1.rxns{i,1},X{2,1});
            d = d+1;
            dup_rxns{d,1} = model1.rxns{c,1}; % % stores original
            dup_rxns{d,2} = model1.rxns{i,1}; % % stores duplicate
            if ismember(model1.rxns{c,1},modelA.rxns)
                belong_dups{d,1} = modelA.description;
            else
                belong_dups{d,1} = modelB.description;
            end
            if ismember(model1.rxns{i,1},modelA.rxns)
                belong_dups{d,2} = modelA.description;
            else
                belong_dups{d,2} = modelB.description;
            end
            % % check if stoichiometry is same
            if sum(roundn(abs(S(:,c))/max(abs(S(:,c))),-2)==roundn(abs(S(:,i))/max(abs(S(:,i))),-2))~=size(S,1)
%                 fprintf('The stoichiometry is different for duplicate.\n');
                s = s+1; % % track all differences in stoichiometry
                diff_s{s,1} = model1.rxns{c,1}; % % stores original
                diff_s{s,2} = model1.rxns{i,1}; % % stores duplicate
            else
                % % check if gene reaction involvement is same
                if (sum(model1.rxnGeneMat(c,:)==model1.rxnGeneMat(i,:))==size(model1.genes,1))
                    if (sum(model1.rxnGeneMat(i,:))~=0)
                        if isempty(strmatch(belong_dups{d,1},belong_dups{d,2}))
%                             fprintf('Different models, same stoichiometry, same gene reaction association and will be removed.\n');
                            g = g+1;
                            kg = kg+1;
                            % % decide which reaction to keep
                            if charge_balance(model1,model1.rxns(i))==1
                                rem_rxns{kg,1} = model1.rxns{i,1}; % % stores original
                                rem_rxns{kg,2} = model1.rxns{c,1}; % % stores duplicate
                                model1 = update_reaction_properties(model1,[i c]); % copy from reaction i to reaction c
                            else
                                rem_rxns{kg,1} = model1.rxns{c,1}; % % stores original
                                rem_rxns{kg,2} = model1.rxns{i,1}; % % stores duplicate
                            end
                            iremrxns(g,1) = i;
                            % % reactions will be removed even if reversibility is different
                            % % check if reversibilities are similar
                            if sum(cellfun(@isempty,strfind(X,'->')))==sum(cellfun(@isempty,strfind(X,'<=>')))
%                                 fprintf('Reversibility does not match, original reaction was made reversible.\n');
                                % % make the original reaction reversible
                                model1.rev(c,1) = 1;
                                model1.lb(c,1) = -1000;
                                model1.ub(c,1) = 1000;
                                r = r+1; % % track all the changed reactions
                                changed_rev{r,1} = model1.rxns{c,1}; % % stores original
                                changed_rev{r,2} = model1.rxns{i,1}; % % stores duplicate
%                             else
%                                 fprintf('Reversibility of the duplicate reaction matches that of original.\n');
                            end
%                         else
%                             fprintf('Duplicate reaction belongs to the same model, and has same gene reaction association.\n');
                        end
%                     else
%                         fprintf('No genes associated with these reactions.\n');
                    end
                else
                    k = k+1;
                    diff_gra{k,1} = model1.rxns{c,1}; % % stores original
                    diff_gra{k,2} = model1.rxns{i,1}; % % stores duplicate
%                     fprintf('The gene reaction association is different and reaction cannot be removed.\n');
                end
            end
        end
    end
    if g~=0
        % % remove respective reaction properties
        model1 = update_reaction_properties(model1,iremrxns);
%         fprintf('Round %d removed %d reactions.\n',c,g);
    end
    c = c+1;
    waitbar(c/steps);
end
close(h);
fprintf('Finished\n');
% % identify any duplicate metabolites
[~,i] = unique(model1.mets);
if ~isempty(setdiff([1:1:length(model1.mets)]',i))
    warning('Metabolite list is not unique. There may be errors.\n');
end
% % identify if all genes have a reaction associated
if ~isempty(model1.genes(sum(model1.rxnGeneMat,1)'==0))
    warning('%d genes have no reaction associated with them.\n',length(model1.genes(sum(model1.rxnGeneMat,1)'==0)));
end

impr.unused_mets = unused_mets; % % list of ununsed metabolite ids (removed)
impr.dup_genes = twice_present; % % genes with public_name converted to molecular_name
impr.genes_not_found = genes_not_found; % % genes which were not found in WormBase or KEGG
impr.dup_rxns = dup_rxns; % % duplicate reactions (all - same products, same reactants)
impr.belong_dups = belong_dups; % % which model do duplicate reactions belong to
impr.changed_rev = changed_rev; % % may need resolving (duplicate reactions and removed reactions subset)
impr.diff_s = diff_s; % % different stoichiometry, needs to be resolved (duplicate reactions subset, not removed)
impr.diff_gra = diff_gra; % % different stoichiometry, needs to be resolved (duplicate reactions subset, not removed)
impr.rem_rxns = rem_rxns; % % removed reactions (duplicate reactions subset)
impr.model = model1;

% removing duplicate reaction pairs that originate from the same model
% assumption: no duplicates exist within the models, onlt amongst the
% models
index = find(cellfun(@strcmp,impr.belong_dups(:,1),impr.belong_dups(:,2)));
impr.dup_rxns(index,:) = [];
impr.belong_dups(index,:) = [];
% removing diff_gra pairs that belong to same model
ic1 = ismember(impr.diff_gra(:,1),modelA.rxns); ic2 = ismember(impr.diff_gra(:,2),modelA.rxns);
el1 = ismember(impr.diff_gra(:,1),modelB.rxns); el2 = ismember(impr.diff_gra(:,2),modelB.rxns);
impr.diff_gra((el1 & el2) | (ic1 & ic2),:) = [];

% removing uptake and exchange reactions from various lists
impr = remove_duplicate_uptake(impr,1);

% removing diff_s pairs that belong to same model or those that have been
% removed from the model
ic1 = ismember(impr.diff_s(:,1),modelA.rxns); ic2 = ismember(impr.diff_s(:,2),modelA.rxns);
el1 = ismember(impr.diff_s(:,1),modelB.rxns); el2 = ismember(impr.diff_s(:,2),modelB.rxns);
i1 = ismember(impr.diff_s(:,1),impr.model.rxns); i2 = ismember(impr.diff_s(:,2),impr.model.rxns);
impr.diff_s((el1 & el2) | (ic1 & ic2) | ~i1 | ~i2,:) = [];

fprintf('\nSUMMARY:\n---------------------------------------------------------\n');
fprintf('%d unused metabolties were found and removed.\n',length(iunmets));
fprintf('%d duplicate genes have been found, combined and removed.\n',size(impr.dup_genes,1));
fprintf('%d genes not found in WormBase or KEGG.\n',size(impr.genes_not_found,1));
fprintf('%d total duplicate reactions were found.\n',size(impr.dup_rxns,1));
fprintf('%d reactions were made reversible and duplicates removed.\n',size(impr.changed_rev,1));

is1 = ismember(impr.dup_rxns(:,1),impr.rem_rxns(:,2));
is2 = ismember(impr.dup_rxns(:,2),impr.rem_rxns(:,1));
tot = [impr.diff_s;impr.diff_gra];
k1 = ismember(impr.dup_rxns(:,1),tot(:,1)); k2 = ismember(impr.dup_rxns(:,2),tot(:,2));
is = find((is1&is2)|(~k1&k2));
impr.dup_rxns(is,:) = [];
impr.belong_dups(is,:) = [];
impr = reconcile_stoichiometry(impr);
% analyzing and resolving gene-reaction
impr = reconcile_gene_reaction(impr);

fprintf('%d duplicate reactions with different gene-reaction still remain.\n',size(impr.diff_gra,1));
fprintf('%d duplicate reactions with different stoichiometry still remain.\n',size(impr.diff_s,1));
fprintf('%d duplicate reactions could be removed.\n',size(impr.rem_rxns,1));
fprintf('%d duplicate reaction pairs still remain within the model.\n',size(impr.dup_rxns,1));

t = toc;