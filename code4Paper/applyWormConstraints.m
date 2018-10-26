function [model,constrainedRxns] = applyWormConstraints(model,c)

%% CONSTRAIN SOME REACTIONS:
bio = model.rxns(~cellfun(@isempty,strfind(model.rxns,'BIO')));
model = changeRxnBounds(model,bio,0,'b');
model = changeRxnBounds(model,'EXC0050',0,'b');

%% GLOBAL CONSTRAINTS:
model.lb(model.lb==-999999) = -1000;
model.ub(model.ub==999999) = 1000;
modelOrig = model;

%% IMPORTS:
% THF
model = changeRxnBounds(model,'R_thf__exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_thf__exchange_LPAREN_e_RPAREN_',0,'u');
% L-alanine
model = changeRxnBounds(model,'EX00041',-c,'l');
model = changeRxnBounds(model,'EX00041',0,'u');
% S-adenosyl-L-methionine
model = changeRxnBounds(model,'EX00019',-c,'l');
model = changeRxnBounds(model,'EX00019',0,'u');
% arginine
model = changeRxnBounds(model,'EX00062',-c,'l');
model = changeRxnBounds(model,'EX00062',0,'u');
% asparagine
model = changeRxnBounds(model,'EX00152',-c,'l');
model = changeRxnBounds(model,'EX00152',0,'u');
% aspartate
model = changeRxnBounds(model,'EX00049',-c,'l');
model = changeRxnBounds(model,'EX00049',0,'u');
% calcium
model = changeRxnBounds(model,'R_Ca2_LPAREN_e_RPAREN_',0,'l');
model = changeRxnBounds(model,'R_Ca2_LPAREN_e_RPAREN_',c,'u');
% chlorine
model = changeRxnBounds(model,'R_cl_LPAREN_e_RPAREN_',0,'l');
model = changeRxnBounds(model,'R_cl_LPAREN_e_RPAREN_',c,'u');
% copper
model = changeRxnBounds(model,'R_cu2_LPAREN_e_RPAREN_',0,'l');
model = changeRxnBounds(model,'R_cu2_LPAREN_e_RPAREN_',c,'u');
% cystiene
model = changeRxnBounds(model,'R_Ex_C_cys_L_LSQBKT_c_RSQBKT__b',0,'l');
model = changeRxnBounds(model,'R_Ex_C_cys_L_LSQBKT_c_RSQBKT__b',c,'u');
% adenosine
model = changeRxnBounds(model,'R_Adenosine_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_Ex_C_cys_L_LSQBKT_c_RSQBKT__b',0,'u');
% cytidine
model = changeRxnBounds(model,'R_cytidine_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_cytidine_exchange_LPAREN_e_RPAREN_',0,'u');
% guanosine
model = changeRxnBounds(model,'R_Guanosine_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_Guanosine_exchange_LPAREN_e_RPAREN_',0,'u');
% thymidine
model = changeRxnBounds(model,'EX00214',-c,'l');
model = changeRxnBounds(model,'EX00214',0,'u');
% ferrous (iron)
model = changeRxnBounds(model,'R_fe2_LPAREN_e_RPAREN_',0,'l');
model = changeRxnBounds(model,'R_fe2_LPAREN_e_RPAREN_',c,'u');
% glutamine
model = changeRxnBounds(model,'R_Ex_C_gln_L_LSQBKT_c_RSQBKT__b',-c,'l');
model = changeRxnBounds(model,'R_Ex_C_gln_L_LSQBKT_c_RSQBKT__b',0,'u');
% glutamate
model = changeRxnBounds(model,'R_Ex_C_glu_L_LSQBKT_c_RSQBKT__b',-c,'l');
model = changeRxnBounds(model,'R_Ex_C_glu_L_LSQBKT_c_RSQBKT__b',0,'u');
% glycine
model = changeRxnBounds(model,'R_Ex_C_gly_LSQBKT_c_RSQBKT__b',-c,'l');
model = changeRxnBounds(model,'R_Ex_C_gly_LSQBKT_c_RSQBKT__b',0,'u');
% h2o
model = changeRxnBounds(model,'R_Water_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_Water_exchange_LPAREN_e_RPAREN_',0,'u');
% histidine
model = changeRxnBounds(model,'R_Histidine_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_Histidine_exchange_LPAREN_e_RPAREN_',0,'u');
% isoleucine
model = changeRxnBounds(model,'R_Isoleucine_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_Isoleucine_exchange_LPAREN_e_RPAREN_',0,'u');
% potassium
model = changeRxnBounds(model,'R_k_LPAREN_e_RPAREN_',0,'l');
model = changeRxnBounds(model,'R_k_LPAREN_e_RPAREN_',c,'u');
% leucine
model = changeRxnBounds(model,'R_Leucine_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_Leucine_exchange_LPAREN_e_RPAREN_',0,'u');
% lysine
model = changeRxnBounds(model,'R_Lysine_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_Lysine_exchange_LPAREN_e_RPAREN_',0,'u');
% methionine
model = changeRxnBounds(model,'R_Methionine_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_Methionine_exchange_LPAREN_e_RPAREN_',0,'u');
% magnesium
model = changeRxnBounds(model,'R_mg2_LPAREN_e_RPAREN_',0,'l');
model = changeRxnBounds(model,'R_mg2_LPAREN_e_RPAREN_',c,'u');
% folate
model = changeRxnBounds(model,'R_folate_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_folate_exchange_LPAREN_e_RPAREN_',0,'u');
% manganese
model = changeRxnBounds(model,'R_mn2_LPAREN_e_RPAREN_',0,'l');
model = changeRxnBounds(model,'R_mn2_LPAREN_e_RPAREN_',c,'u');
% phenylalanine
model = changeRxnBounds(model,'R_Phenylalanine_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_Phenylalanine_exchange_LPAREN_e_RPAREN_',0,'u');
% pheme
model = changeRxnBounds(model,'R_Heme_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_Heme_exchange_LPAREN_e_RPAREN_',0,'u');
% proline
model = changeRxnBounds(model,'R_Ex_C_pro_L_LSQBKT_c_RSQBKT__b',0,'l');
model = changeRxnBounds(model,'R_Ex_C_pro_L_LSQBKT_c_RSQBKT__b',c,'u');
% pyridoxine
model = changeRxnBounds(model,'R_pydxn_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_pydxn_LPAREN_e_RPAREN_',0,'u');
% pyridoxamine
model = changeRxnBounds(model,'R_pydam_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_pydam_LPAREN_e_RPAREN_',0,'u');
% pyridoxal
model = changeRxnBounds(model,'R_pydx_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_pydx_LPAREN_e_RPAREN_',0,'u');
% riboflavin
model = changeRxnBounds(model,'R_Riboflavin_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_Riboflavin_exchange_LPAREN_e_RPAREN_',0,'u');
% serine
model = changeRxnBounds(model,'R_Ex_C_ser_L_LSQBKT_c_RSQBKT__b',-c,'l');
model = changeRxnBounds(model,'R_Ex_C_ser_L_LSQBKT_c_RSQBKT__b',0,'u');
% sulfate
model = changeRxnBounds(model,'EX00059',-c,'l');
model = changeRxnBounds(model,'EX00059',0,'u');
% thiamin
model = changeRxnBounds(model,'R_thiamin_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_thiamin_exchange_LPAREN_e_RPAREN_',0,'u');
% threonine
model = changeRxnBounds(model,'R_Threonine_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_Threonine_exchange_LPAREN_e_RPAREN_',0,'u');
% tryptophan
model = changeRxnBounds(model,'R_Tryptophan_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_Tryptophan_exchange_LPAREN_e_RPAREN_',0,'u');
% tyrosine
model = changeRxnBounds(model,'R_Ex_C_tyr_L_LSQBKT_c_RSQBKT__b',-c,'l');
model = changeRxnBounds(model,'R_Ex_C_tyr_L_LSQBKT_c_RSQBKT__b',0,'u');
% valine
model = changeRxnBounds(model,'R_Valine_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_Valine_exchange_LPAREN_e_RPAREN_',0,'u');
% 5-methyl-tetrahydrofolate
model = changeRxnBounds(model,'EX00440',-c,'l');
model = changeRxnBounds(model,'EX00440',0,'u');
% coenzyme B12
model = changeRxnBounds(model,'R_Vitamin_B12_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_Vitamin_B12_exchange_LPAREN_e_RPAREN_',0,'u');
% biotin
model = changeRxnBounds(model,'R_biotin_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_biotin_exchange_LPAREN_e_RPAREN_',0,'u');
% glutathione
model = changeRxnBounds(model,'R_Ex_gthrdT_',-c,'l');
model = changeRxnBounds(model,'R_Ex_gthrdT_',0,'u');
% lipoate
model = changeRxnBounds(model,'R_LIPOIC_ACID_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_LIPOIC_ACID_LPAREN_e_RPAREN_',0,'u');
% uridine
model = changeRxnBounds(model,'R_Uridine_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_Uridine_exchange_LPAREN_e_RPAREN_',0,'u');
% D-alanine
model = changeRxnBounds(model,'EX00133',-c,'l');
model = changeRxnBounds(model,'EX00133',0,'u');
% N-acetyl-D-glucosamine
model = changeRxnBounds(model,'R_N_Acetylglucosamine_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_N_Acetylglucosamine_exchange_LPAREN_e_RPAREN_',0,'u');
% p-aminobenzoate
model = changeRxnBounds(model,'R_4ABZ_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_4ABZ_LPAREN_e_RPAREN_',0,'u');
% oxygen
model = changeRxnBounds(model,'R_O2_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_O2_exchange_LPAREN_e_RPAREN_',0,'u');
% nicotinate
model = changeRxnBounds(model,'R_Niacin_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_Niacin_exchange_LPAREN_e_RPAREN_',0,'u');
% nicotinamide
model = changeRxnBounds(model,'R_Niacinamide_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_Niacinamide_exchange_LPAREN_e_RPAREN_',0,'u');
% phosphate
model = changeRxnBounds(model,'R_Pi_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_Pi_exchange_LPAREN_e_RPAREN_',0,'u');
% sodium
model = changeRxnBounds(model,'R_na1_LPAREN_e_RPAREN_',0,'l');
model = changeRxnBounds(model,'R_na1_LPAREN_e_RPAREN_',c,'u');
% hydroxid ion
model = changeRxnBounds(model,'R_OH_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_OH_exchange_LPAREN_e_RPAREN_',0,'u');
% proton
model = changeRxnBounds(model,'R_Proton_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_Proton_exchange_LPAREN_e_RPAREN_',100,'u');
% myo-inositol
model = changeRxnBounds(model,'R_myoinositol_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_myoinositol_exchange_LPAREN_e_RPAREN_',0,'u');
% glucose
model = changeRxnBounds(model,'R_EX_GLC_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_EX_GLC_LPAREN_e_RPAREN_',0,'u');
% choline
model = changeRxnBounds(model,'R_Choline_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_Choline_exchange_LPAREN_e_RPAREN_',0,'u');
% citrate
model = changeRxnBounds(model,'R_Citrate_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_Citrate_exchange_LPAREN_e_RPAREN_',0,'u');
% sitosterol
model = changeRxnBounds(model,'R_Sitosterol_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_Sitosterol_exchange_LPAREN_e_RPAREN_',0,'u');
% glycogenin
model = changeRxnBounds(model,'R_Glycogenin_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_Glycogenin_exchange_LPAREN_e_RPAREN_',0,'u');
% peptide
model = changeRxnBounds(model,'R_Peptides_exchange_LPAREN_e_RPAREN_',-c,'l');
model = changeRxnBounds(model,'R_Peptides_exchange_LPAREN_e_RPAREN_',0,'u');

%% EXPORTS:
% CO2
model = changeRxnBounds(model,'R_CO2_exchange_LPAREN_e_RPAREN_',0,'l');
model = changeRxnBounds(model,'R_CO2_exchange_LPAREN_e_RPAREN_',1000,'u');
% urea
model = changeRxnBounds(model,'R_Urea_exchange_LPAREN_e_RPAREN_',0,'l');
model = changeRxnBounds(model,'R_Urea_exchange_LPAREN_e_RPAREN_',1000,'u');

constrainedRxns = model.rxns(model.lb~=modelOrig.lb | model.ub~=modelOrig.ub);