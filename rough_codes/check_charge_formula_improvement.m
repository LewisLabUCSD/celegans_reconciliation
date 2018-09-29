function [unatom,uncharge,model1] = check_charge_formula_improvement(model,rxns)

formula_list = printRxnFormula(model,rxns,false,true,false,1,true,true);
met_list = convert_reaction_formula2list(formula_list);
prop_mets = metprop_BiGG(met_list,'all');
model1 = model;
for i=1:length(met_list)
    old_formula{i,1} = model.metFormulas{strcmp(model.mets,met_list{i,1}),1};
    new_formula{i,1} = prop_mets.formula{i,1};
    if strcmp(prop_mets.formula{i,1},'not found')
        model1.metFormulas{strcmp(model.mets,met_list{i,1}),1} = ' ';
    elseif ~strcmp(old_formula{i,1},new_formula{i,1})
        model1.metFormulas{strcmp(model.mets,met_list{i,1}),1} = prop_mets.formula{i,1};
    end
end
for i=1:length(met_list)
    old_charge(i,1) = model.metCharge(strcmp(model.mets,met_list{i,1}));
    new_charge(i,1) = prop_mets.charge{i,1};
    model1.metCharge(strcmp(model.mets,met_list{i,1})) = int32(prop_mets.charge{i,1});
end

old_stat = atom_balance(model,rxns,false);
rxns1 = rxns(old_stat~=1);
if ~isempty(rxns1)
    new_stat = atom_balance(model1,rxns1,false);
    unatom = rxns1(new_stat~=1);
else
    unatom = [];
end
old_stat = charge_balance(model,rxns,false);
rxns1 = rxns(old_stat~=1);
if ~isempty(rxns1)
    new_stat = charge_balance(model1,rxns1,false);
    uncharge = rxns1(new_stat~=1);
else
    uncharge = [];
end
uncharge = charge_balance(model1,rxns,false);
unatom = atom_balance(model1,rxns,false);