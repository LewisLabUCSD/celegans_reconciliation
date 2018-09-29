function model = update_metabolite_properties(model,index)
% model = update_metabolite_properties(model,index)
% Changes only metabolite properties. Must change model.S and model.mets,
% before

% INPUT:
% model: model to be changed
% index: index, max(index) < size(model.mets,1)
%
% OUTPUT:
% model: changed model
%
metFields = listMetaboliteFields(model);
metFields(ismember(metFields,{'mets';'S'})) = [];
if max(index) <= size(model.mets,1)
    % % remove respective metabolite properties
    model.mets(index) = []; % metabolite ids
    model.S(index,:) = []; % stoichiometric matrix
    if ~isempty(metFields)
        for i=1:length(metFields)
            eval(strcat('model.',metFields{i,1},'(index) = [];'));
        end
    end
    %     if isfield(model,'metNames') % metabolite name
    %         model.metNames(index) = [];
    %     end
    %     if isfield(model,'metFormulas') % metabolite formula
    %         model.metFormulas(index) = [];
    %     end
    %     if isfield(model,'metCharges') % metabolite charge
    %         model.metCharges(index) = [];
    %     end
    %     if isfield(model,'b') % metabolite accumulation
    %         model.b(index) = [];
    %     end
    %     if isfield(model,'csense')
    %         model.csense(index) = [];
    %     end
    %     if isfield(model,'metKEGGID') % metabolite KEGG
    %         model.metKEGGID(index) = [];
    %     end
    %     if isfield(model,'metInChIString') % metabolite InChIString
    %         model.metInChIString(index) = [];
    %     end
    %     if isfield(model,'metChEBIID') % metabolite ChEBI
    %         model.metChEBIID(index) = [];
    %     end
    %     if isfield(model,'metPubChemID') % metabolite PubChem
    %         model.metPubChemID(index) = [];
    %     end
else
    warning('Model not updated, please double check indices vector');
end