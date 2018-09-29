function [constrain_old,constrain_new,del_match] = check_gene_account(model_old,model_new,twice_present)

for i=1:size(twice_present,1)
    constrain_old1 = model_old.rxns(find(model_old.rxnGeneMat(:,strcmp(model_old.genes,twice_present{i,1}))));
    constrain_old2 = model_old.rxns(find(model_old.rxnGeneMat(:,strcmp(model_old.genes,twice_present{i,2}))));
    constrain_new = model_new.rxns(find(model_new.rxnGeneMat(:,strcmp(model_new.genes,twice_present{i,2}))));
    constrain_old = union(constrain_old1,constrain_old2);
    s(i,1) = sum(ismember(constrain_new,constrain_old));
    s(i,2) = sum(ismember(constrain_old,constrain_new));
    if s(i,1) == s(i,2)
        del_match(i,1) = true;
    else
        del_match(i,1) = false;
    end
end
    