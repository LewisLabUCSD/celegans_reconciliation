function [impr1,rems] = remove_duplicate_uptake(impr,rem_col)
% removing uptake and exchange reactions from various lists

[exch_rxns,up_rxns] = findExcRxns(impr.model,true,false);
exch_rxns = impr.model.rxns(exch_rxns);
up_rxns = impr.model.rxns(up_rxns);
upr(:,1) = ismember(impr.dup_rxns(:,1),up_rxns);
upr(:,2) = ismember(impr.dup_rxns(:,2),up_rxns);
exchr(:,1) = ismember(impr.dup_rxns(:,1),exch_rxns);
exchr(:,2) = ismember(impr.dup_rxns(:,2),exch_rxns);

rem_4mdups = find((upr(:,1) & exchr(:,1)) | (exchr(:,2) & upr(:,2)));
rem_dups = find((upr(:,1) & exchr(:,2)) | (exchr(:,1) & upr(:,2)));
rems = intersect(rem_4mdups,rem_dups);
if ~isempty(rems)
    impr1 = impr;
    impr1.version0 = impr;
    impr1.dup_rxns(rems,:) = [];
    impr1.belong_dups(rems,:) = [];
    rems = impr.dup_rxns(rems,:);
    if ~isempty(rem_dups)
        rem_dups = impr.dup_rxns(rem_dups,:);
        r = ismember(impr.model.rxns,rem_dups(:,rem_col));
        impr1.model = update_reaction_properties(impr.model,find(r));
        impr1.rem_rxns = [impr1.rem_rxns; rem_dups];
    end
    igra1 = ismember(impr1.diff_gra(:,1),rems(:,1)); is1 = ismember(impr1.diff_s(:,1),rems(:,1));
    igra2 = ismember(impr1.diff_gra(:,2),rems(:,2)); is2 = ismember(impr1.diff_s(:,2),rems(:,2));
    igra = igra1 & igra2; is = is1 & is2;
    impr1.diff_gra(igra,:) = []; impr1.diff_s(is,:) = [];
end