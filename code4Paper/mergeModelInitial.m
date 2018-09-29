function [modelNew, otherFields] = mergeModelInitial(model1,model2,objRxn)
% creates an initial merge of the model

% written by chintan joshi at ucsd on 09/17/2017

% % get fieldnames % %
fdNames = intersect(fieldnames(model1),fieldnames(model2));
fdNames(strcmp(fdNames,'S')) = [];
fdNames(strcmp(fdNames,'rxnGeneMat')) = [];
fdNames(strcmp(fdNames,'osense')) = []; modelNew.osense = -1;
fdNames(strcmp(fdNames,'description')) = []; 
modelNew.description = [regexprep(model1.description,'.xml',''),' + ',regexprep(model2.description,'.xml','')];
fdNames(strcmp(fdNames,'modelVersion')) = []; % currently not being assigned, can be assigned later.

% % identify fieldnames associated with reactions, metabolites, genes, &
% % compartments % %
metFields = []; rxnFields = []; compFields = []; geneFields = []; otherFields = [];
m = 0; r = 0; g = 0; c = 0; o = 0;
for i=1:length(fdNames)
    if eval(strcat('size(model1.',fdNames{i,1},',1) == size(model1.mets,1)'))
        m = m + 1;
        metFields{m,1} = fdNames{i,1};
    elseif eval(strcat('size(model1.',fdNames{i,1},',1) == size(model1.rxns,1)'))
        r = r + 1;
        rxnFields{r,1} = fdNames{i,1};
    elseif eval(strcat('size(model1.',fdNames{i,1},',1) == size(model1.comps,1)'))
        c = c + 1;
        compFields{c,1} = fdNames{i,1};
    elseif eval(strcat('size(model1.',fdNames{i,1},',1) == size(model1.genes,1)'))
        g = g + 1;
        geneFields{g,1} = fdNames{i,1};
    else
        o = o + 1;
        otherFields{o,1} = fdNames{i,1};
    end
end

% % get reaction properties from models % %
fprintf('Combining reaction properties...');
[additionalRxns,ir] = setdiff(model2.rxns,model1.rxns);
modelNew.rxns = vertcat(model1.rxns,additionalRxns);
rxnFields(strcmp(rxnFields,'rxns')) = [];
if ~isempty(rxnFields)
    for i=1:length(rxnFields)
        eval(strcat('modelNew.',rxnFields{i,1},' = vertcat(model1.',rxnFields{i,1},',model2.',rxnFields{i,1},'(ir));'));
    end
end
% modelNew.rxnNames = vertcat(model1.rxnNames,model2.rxnNames(ir));
% modelNew.lb = vertcat(model1.lb,model2.lb(ir));
% modelNew.ub = vertcat(model1.ub,model2.ub(ir));
% modelNew.c = vertcat(model1.c,model2.c(ir));
% modelNew.subSystems = vertcat(model1.subSystems,model2.subSystems(ir));
% modelNew.rxnNotes = vertcat(model1.rxnNotes,model2.rxnNotes(ir));
% modelNew.rxnConfidenceScores = vertcat(model1.rxnConfidenceScores,model2.rxnConfidenceScores(ir));
% modelNew.rxnECNumbers = vertcat(model1.rxnECNumbers,model2.rxnECNumbers(ir));
fprintf('Finished!\n');

% % get metabolite properties from models % %
fprintf('Combining metabolite properties...');
[additionalMets,im] = setdiff(model2.mets,model1.mets);
modelNew.mets = vertcat(model1.mets,additionalMets);
metFields(strcmp(metFields,'mets')) = [];
if ~isempty(metFields)
    for i=1:length(metFields)
        eval(strcat('modelNew.',metFields{i,1},' = vertcat(model1.',metFields{i,1},',model2.',metFields{i,1},'(im));'));
    end
end
% modelNew.metNames = vertcat(model1.mets,model2.metNames(im));
% modelNew.metFormulas = vertcat(model1.metFormulas,model2.metFormulas(im));
% modelNew.metCharges = vertcat(model1.metCharges,model2.metCharges(im));
% modelNew.b = vertcat(model1.b,model2.b(im));
% modelNew.csense = vertcat(model1.csense,model2.csense(im));
% modelNew.metNotes = vertcat(model1.metNotes,model2.metNotes(im));
fprintf('Finished!\n');

% % get compartments from models % %
fprintf('Combining compartment properties...');
[additionalComps,ic] = setdiff(model2.comps,model1.comps);
modelNew.comps = vertcat(model1.comps,additionalComps);
compFields(strcmp(compFields,'comps')) = [];
if ~isempty(compFields)
    for i=1:length(compFields)
        eval(strcat('modelNew.',compFields{i,1},' = vertcat(model1.',compFields{i,1},',model2.',compFields{i,1},'(ic));'));
    end
end
% modelNew.compNames = vertcat(model1.comps,model2.compNames(ic));
fprintf('Finished!\n');

% % make stoichiometric matrix % %
fprintf('Constructing stoichiometric matrix...');
modelNew.S = zeros(length(modelNew.mets),length(modelNew.rxns));
modelNew.S(1:length(model1.mets),1:length(model1.rxns)) = model1.S;
for i=1:length(additionalRxns)
    rxnMets = find(model2.S(:,ir(i))~=0);
    metabs = model2.mets(rxnMets);
    for j=1:length(rxnMets)
        modelNew.S(strcmp(modelNew.mets,metabs{j}),length(model1.rxns)+i) = model2.S(rxnMets(j),ir(i));
    end
end
fprintf('Finished!\n');

% % make gene properties from the model % %
fprintf('Combining genes...');
[additionalGenes,ig] = setdiff(model2.genes,model1.genes);
modelNew.genes = vertcat(model1.genes,additionalGenes);
model1 = updateModelGenes(model1,modelNew.genes);
model2 = updateModelGenes(model2,modelNew.genes);
geneFields(strcmp(geneFields,'genes')) = [];
if ~isempty(geneFields)
    for i=1:length(geneFields)
        eval(strcat('modelNew.',geneFields{i,1},' = vertcat(model1.',geneFields{i,1},',model2.',geneFields{i,1},'(ig));'));
    end
end
fprintf('Finished!\n');

% % make rules of the model % %
fprintf('Creating rules and grRules...');
modelNew.grRules = vertcat(model1.grRules,model2.grRules(ir));
modelNew.rules = vertcat(model1.rules,model2.rules(ir));
fprintf('Finished!\n');

% % creating rxnGeneMat of the model % %
modelNew = rxnGeneMatModel(modelNew);

% % change/set objective
if nargin > 2
    modelNew = changeObjective(modelNew,objRxn);
end
fprintf('Finished!\n');

% % printing other details of the model % %
fprintf('Combined model name is...');
fprintf('%s\n',modelNew.description);
fprintf('# of reactions = %d\n',length(modelNew.rxns));
fprintf('# of metabolites = %d\n',length(modelNew.mets));
fprintf('# of genes = %d\n',length(modelNew.genes));