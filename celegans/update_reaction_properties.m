function model = update_reaction_properties(model,index)
% model = update_reaction_properties(model,index)
% Changes only reaction properties. Must change model.S and model.rxns,
% before

% INPUT:
% model: model to be changed
% index: index, max(index) <= size(model.rxns,1)
%
% OUTPUT:
% model: changed model
%
if max(index) <= size(model.rxns,1)
    if isfield(model,'rxnNames') % % reaction names
        model.rxnNames(index) = [];
    end
    if isfield(model,'lb')
        model.lb(index) = []; % % reaction lower bounds
        model.ub(index) = []; % % reaction upper bounds
    end
    if isfield(model,'subSystems') % % reaction subsystem
        model.subSystems(index) = [];
    end
    if isfield(model,'rxnReferences') % % reaction references
        model.rxnReferences(index) = [];
    end
    if isfield(model,'rxnECNumbers') % % reaction EC numbers
        model.rxnECNumbers(index) = [];
    end
    if isfield(model,'rxnNotes') % % reaction notes
        model.rxnNotes(index) = [];
    end
    if isfield(model,'c') % % reaction objective coefficient
        model.c(index) = [];
    end
    if isfield(model,'rev') % % reaction reversibility
        model.rev(index) = [];
    end
    if isfield(model,'rxnGeneMat') % % reaction gene association
        model.rxnGeneMat(index,:) = [];
    end
    if isfield(model,'grRules') % % reaction gene rules
        model.grRules(index) = [];
    end
    if isfield(model,'rules') % % reaction rules
        model.rules(index) = [];
    end
    if isfield(model,'confidenceScores') % % reaction confidence scores
        model.confidenceScores(index) = [];
    end
else
    warning('Model not updated, please double check indices vector');
end