This folder contains all the codes required for generating an automated reconciliation for C. elegans model used for further curation.

The following lines will generate the automated reconciled model:
%%%%%%%
load('yilmaz_1273_new.mat','icel');
load('kaleta_new.mat','eleg');
[model1,model_unref,duplicateBelongs,allIndices,remainDuplicates] = refine_merged_model_v3(eleg,icel,'R_Biomass2',true);
%%%%%%%

FILES INFO:
Celegans_axenic.xml
iCEL1273_2.xml
The above two files are XML files for the previously published models used for the automated reconciliation.

kaleta_new.mat
yilmaz_1273_2.mat
The above two files contain manual ID conversions of genes (WBGene ids), reactions (BiGG format) and metabolites (BiGG format) for previously published models.

mergedModel.mat
Contains output of the code above. The merged model is model1 in the mat file. To access the model use the code below
%%%%%%%
load('mergedModel.mat','model1');
%%%%%%%

readElegModel.m
readICELModel.m
The above two files take conversions from the mat files and reads the above two existing models in desired model identifiers.

changeObjective.m
findRxnIDs.m
printRxnFormula.m
The above files are from COBRA toolbox and present here to facilitate using reconciliation codes without having access to COBRA toolbox.

c_elegans.canonical_bioproject.current
The list of files used for generating a list of wormbase genes and other identifiers associated with that genes

Celegans_axenic.xml
iCEL1273_2.xml
readElegModel
readICELModel
changeObjective
findRxnIDs
printRxnFormula
c_elegans.canonical_bioproject.current
atom_balance: 						generates elemental balance for a list of reactions
create_atom_matrix:					creates a cell matrix which is used to determine atom balance
create_atom_numbers:				gives number of atoms for elements in a given metbolite formula
getallgenes_KEGG:					finds all genes which belong to an organism from KEGG
getgeneinfo_KEGG:					gets information for a given list of genes from KEGG
getgeneinfo_WormBase:				gets gene information for a given list of genes from WormBase
update_gene_properties:				changes gene properties of a list of genes in the model
updateModelGenes:					updates & changes the list of genes in the model
getchargeinfo_MNX:					gets charge information for a metabolite from MetaNetX using a web link
getformulainfo_MNX:					gets formula information for a metabolite from MetNetX using a web link
metprop_BiGG:						finds a given information type for a metabolite in BiGG
update_metabolite_properties:		changes metabolite properties of a list of metabolites within a model
charge_balance:						finds the net charge remaining on the reactions, and also a if the reaction was charge balanced
check_charge_formula_improvement:	checks the improvements for charge and formula
convert_reaction_formula2list:		estimates the list of metabolites based on the formula provided
update_reaction_properties:			changes reaction properties of a list of reactions within a model
grRulesModel:						creates geRules field for a given model structure
rulesModel:							create rules field for a given model structure
rxnGeneMatModel:					create rxnGeneMat field for a given model structure
listGeneFields:						gets the list of gene-related fields used in the model structure
listMetaboliteFields:				gets the list of metabolite-related fields in the model structure
listReactionFields:					gets the list of reaction-related fields in th emodel structure
mergeModelInitial:					performs an initial merging of the model without identification of improvements or duplicates
reconcile_gene_association:			reconciles and changes gene properties within the model of two reactions
refine_merged_model_v3:				merge and refine model including identification of improvements and duplicates