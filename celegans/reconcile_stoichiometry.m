function impr2 = reconcile_stoichiometry(impr)
% impr2 = reconcile_stoichiometry(impr)
% reconciles the stoichiometry of the model after it has been
% simple-processed(!?)
%
% INPUT:
% impr: improvement structure
%
% OUTPUT:
% impr2: improvement structure after reconciling stoichiometry

rxn_bal(:,1) = atom_balance(impr.model,impr.diff_s(:,1),false);
rxn_bal(:,2) = atom_balance(impr.model,impr.diff_s(:,2),false);

r10 = rxn_bal(:,1)==0; r21 = rxn_bal(:,2)==1;
r20 = rxn_bal(:,2)==0; r11 = rxn_bal(:,1)==1;

% reactions to remove and update impr to impr1
impr1 = impr;
rxns = [impr.diff_s(r11 & r20,1);impr.diff_s(r10 & r21,2)];
rindex = find(ismember(impr.model.rxns,rxns));
impr1.model = update_reaction_properties(impr.model,rindex);
% add the reaction pair to impr.rem_rxns
impr1.rem_rxns = [impr1.rem_rxns;impr1.diff_s((r11&r20)|(r10&r21),:)];
% remove the reaction pair from impr.diff_s
impr1.diff_s((r11&r20)|(r10&r21),:) = [];
% removing the reaction pairs from initially identified pair
tot = [impr1.diff_s;impr1.diff_gra];
k1 = ismember(impr1.dup_rxns(:,1),tot(:,1)); k2 = ismember(impr1.dup_rxns(:,2),tot(:,2));
impr1.dup_rxns(~(k1&k2),:) = []; 
impr1.belong_dups(~(k1&k2),:) = [];
k1 = []; k2 = [];

% removing the pairs for which original reaction (column 1) has been
% removed
impr2 = impr1;
bal(:,1) = atom_balance(impr1.model,impr1.diff_s(:,1),false);
rem_index = bal==-2;
% remove bal reactions from diff_s
% add the reaction pair to impr.rem_rxns
impr2.rem_rxns = [impr2.rem_rxns;impr2.diff_s(rem_index,:)];
% remove the reaction pair from impr.diff_s
impr2.diff_s(rem_index,:) = [];
% removing the reaction pairs from initially identified pair
tot = [impr2.diff_s;impr2.diff_gra];
k1 = ismember(impr2.dup_rxns(:,1),tot(:,1)); k2 = ismember(impr2.dup_rxns(:,2),tot(:,2));
impr2.dup_rxns(~(k1&k2),:) = []; 
impr2.belong_dups(~(k1&k2),:) = [];