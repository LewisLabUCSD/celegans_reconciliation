function [charge_imb,charge] = charge_balance(model,rxns,disp_flag)
% [charge_imb,charge] = charge_balance(model,rxns,disp_flag)
% finds the net charge remaining on the reactions, and also a boolean-like
% vector suggesting if the reaction was charge balanced

% INPUT:
% model: model of interest
% rxns: list of reaction ids belonging to the model
% disp_flag: print on command window (true)
%
% OUTPUT:
% charge_imb: a boolean-like vector suggesting if the reaction was charge balanced
% charg: the net charge remaining on the reaction

if nargin < 2
    rxns = model.rxns;
end
if nargin < 3
    disp_flag = false;
end
charge_imb = false;
charge = zeros(length(rxns),1);
for i=1:length(rxns)
    if ~isempty(strcmp(model.rxns,rxns{i,1}))
        charge(i) = sum(full(model.S(:,strcmp(model.rxns,rxns{i,1}))'*double(model.metCharge)));
        if charge(i)==0
            charge_imb(i,1) = 1;
        else
            charge_imb(i,1) = 0;
        end
    else
        charge(i) = -100000;
        charge_imb(i,1) = -1;
    end
end
if disp_flag
    if sum(charge_imb)==size(rxns,1)
        fprintf('All reaction(s) is(are) charge balanced.\n');
    else
        fprintf('%d reaction(s) not charge balanced.\nReaction id\t\tNet charge\n---------------------------\n',sum(charge_imb))
        for i=1:length(rxns)
            fprintf('%s\t\t%d\n',rxns{i,1},charge(i));
        end
    end
end
charge = charge';