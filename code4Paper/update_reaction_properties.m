function model = update_reaction_properties(model,index)
% model = update_reaction_properties(model,index)
% Changes only reaction properties.

% INPUT:
% model: model to be changed
% index: index, max(index) <= size(model.rxns,1)
%
% OUTPUT:
% model: changed model
%
rxnFields = listReactionFields(model);
if max(index) <= size(model.rxns,1)
    if size(index,2) > 1
        % % change reaction properties from reaction 1 to reaction 2
        for i=1:length(rxnFields)
            eval(strcat('model.',rxnFields{i},'(index(2)) = model.',rxnFields{i},'(index(1));'));
            if strcmp(rxnFields{i},'rxnGeneMat')
                eval(strcat('model.',rxnFields{i},'(index(2),:) = model.',rxnFields{i},'(index(1),:);'));
            end
        end
        model.S(:,index(2)) = model.S(:,index(1)); % % stoichiometric matrix
        
%         model.rxns(index(2)) = model.rxns(index(1)); % % reaction ids
%         if isfield(model,'rxnNames') % % reaction names
%             model.rxnNames(index(2)) = model.rxnNames(index(1));
%         end
%         if isfield(model,'lb')
%             model.lb(index(2)) = model.lb(index(1)); % % reaction lower bounds
%             model.ub(index(1)) = model.ub(index(1)); % % reaction upper bounds
%         end
%         if isfield(model,'subSystems') % % reaction subsystem
%             model.subSystems(index(2)) = model.subSystems(index(1));
%         end
%         if isfield(model,'rxnReferences') % % reaction references
%             model.rxnReferences(index(2)) = model.rxnReferences(index(1));
%         end
%         if isfield(model,'rxnECNumbers') % % reaction EC numbers
%             model.rxnECNumbers(index(2)) = model.rxnECNumbers(index(1));
%         end
%         if isfield(model,'rxnNotes') % % reaction notes
%             model.rxnNotes(index(2)) = model.rxnNotes(index(1));
%         end
%         if isfield(model,'c') % % reaction objective coefficient
%             model.c(index(2)) = model.c(index(1));
%         end
%         if isfield(model,'rev') % % reaction reversibility
%             model.rev(index(2)) = model.rev(index(1));
%         end
%         if isfield(model,'rxnGeneMat') % % reaction gene association
%             model.rxnGeneMat(index(2),:) = model.rxnGeneMat(index(1),:);
%         end
%         if isfield(model,'grRules') % % reaction gene rules
%             model.grRules(index(2)) = model.grRules(index(1));
%         end
%         if isfield(model,'rules') % % reaction rules
%             model.rules(index(2)) = model.rules(index(2));
%         end
%         if isfield(model,'confidenceScores') % % reaction confidence scores
%             model.confidenceScores(index(2)) = model.confidenceScores(index(1));
%         end
    else
        % % remove respective reaction properties
        for i=1:length(rxnFields)
            if strcmp(rxnFields{i},'rxnGeneMat')
                model.rxnGeneMat(index,:) = [];
            else
                eval(strcat('model.',rxnFields{i},'(index) = [];'));
            end
        end
        model.S(:,index) = []; % % stoichiometric matrix
        
%         model.rxns(index) = []; % % reaction ids
%         if isfield(model,'rxnNames') % % reaction names
%             model.rxnNames(index) = [];
%         end
%         if isfield(model,'lb')
%             model.lb(index) = []; % % reaction lower bounds
%             model.ub(index) = []; % % reaction upper bounds
%         end
%         if isfield(model,'subSystems') % % reaction subsystem
%             model.subSystems(index) = [];
%         end
%         if isfield(model,'rxnReferences') % % reaction references
%             model.rxnReferences(index) = [];
%         end
%         if isfield(model,'rxnECNumbers') % % reaction EC numbers
%             model.rxnECNumbers(index) = [];
%         end
%         if isfield(model,'rxnNotes') % % reaction notes
%             model.rxnNotes(index) = [];
%         end
%         if isfield(model,'c') % % reaction objective coefficient
%             model.c(index) = [];
%         end
%         if isfield(model,'rev') % % reaction reversibility
%             model.rev(index) = [];
%         end
%         if isfield(model,'rxnGeneMat') % % reaction gene association
%             model.rxnGeneMat(index,:) = [];
%         end
%         if isfield(model,'grRules') % % reaction gene rules
%             model.grRules(index) = [];
%         end
%         if isfield(model,'rules') % % reaction rules
%             model.rules(index) = [];
%         end
%         if isfield(model,'confidenceScores') % % reaction confidence scores
%             model.confidenceScores(index) = [];
%         end
    end
else
    warning('Model not updated, please double check indices vector');
end