function metBiGG = getbiggChargeFormula(all_biggid, disp)
% metBiGG = getbiggChargeFormula(all_bigg, 1)
% INPUT: 
% all_bigg: list of Bigg compound ids
% disp: flag - 1 to show outpput - 0 - not to show 
% OUTPUT:
% metBiGG: list of Chemical Formulas and Charges from
%          Bigg 
tic;
metBiGG = cell(length(all_biggid),9);
options = weboptions('Timeout',30);

notInbigg = 0;
noCharge=0;
noFormula=0;
inchiFound = 0;
chResolved = 0;
nochResolved = 0;
formResolved = 0;
noformResolved = 0;

%check all Bigg IDs
for i=1:length(all_biggid)
    metBiGG{i,1}= all_biggid{i,1};
                
    data= webread(strcat('http://bigg.ucsd.edu/api/v2/universal/metabolites/', all_biggid{i,1}(1:end-3)),options);

    if isfield(data.database_links,'BioCyc')
        data.database_links.BioCyc.id
        metaCyc = webread( strcat('http://websvc.biocyc.org/getxml?id=', data.database_links.BioCyc.id), options);
        regex = regexp(metaCyc, '<gibbs-0 datatype="float" units="kcal\/mol"\>(-?\d+\.\d+)', 'tokens');
        if isempty(regex)~=1
            gibbs = regex{1};
            display(gibbs);
        end
    end
    try
        %if not same id as previous one

           
            
            data= webread(strcat('http://bigg.ucsd.edu/api/v2/universal/metabolites/', all_biggid{i,1}(1:end-3)),options);
            
            if isfield(data.database_links,'BioCyc')
                metaCyc = webread( data.database_links.BioCyc.link, options);
                regex = regexp(metaCyc, '<gibbs-0 datatype="float" units="kcal\/mol"\>(-?\d+\.\d+)', 'tokens');
                if isempty(regex)~=1
                    gibbs = regex{1};
                    display(gibbs);
                    metBiGG{i,9} = gibbs;
                end
            end
            
            %if there is a charge in BIGG
            if isempty(data.charges)~=1
                
                %getting all charges
                for k=1:length(data.charges)
                    if k == length(data.charges)
                        metBiGG{i,2} = strcat(metBiGG{i,2}, num2str(data.charges(k,1)));
                    else
                        metBiGG{i,2} = strcat(metBiGG{i,2}, num2str(data.charges(k,1)), ', ');
                        metBiGG{i,4} = length(data.charges);
                       
                        if isfield(data.database_links,'MetaNetX_MNX_Chemical')
                            metaNetX = webread( data.database_links.MetaNetX_MNX_Chemical.link);
                            regex = regexp(metaNetX, 'InChI=(.+)</td></tr>   <tr><td', 'tokens');
                            if isempty(regex)~=1
                                inchi = regex{1};
                            end
                        end
                        if isempty(inchi)~=1
                             metBiGG{i,6} = inchi;
                             inchiFound = inchiFound +1;
                        end
                        
                        if  isempty(metBiGG{i,6})~=1
                            metBiGG{i,7} = getChargeFromInChI(metBiGG{i,6});
                            chResolved = chResolved+1;
                        end
                    end
                end 
            else
                if isfield(data.database_links,'MetaNetX_MNX_Chemical')
                    metaNetX = webread( data.database_links.MetaNetX_MNX_Chemical.link);
                    regex = regexp(metaNetX, 'InChI=(.+)</td></tr>   <tr><td', 'tokens');
                    if isempty(regex)~=1
                        inchi = regex{1};
                    end
                end
                if isempty(inchi)~=1
                     metBiGG{i,6} = inchi;
                     inchiFound = inchiFound +1;
                end
                
                noCharge = noCharge+1;
                metBiGG{i,2} = 'No Charge In BiGG';
                if  isempty(metBiGG{i,6})~=1
                    metBiGG{i,7} = getChargeFromInChI(metBiGG{i,6});
                    nochResolved = nochResolved+1;
                end
            end
               
            %if there is a Formula in BIGG
            if isempty(data.formulae)~=1
                
                %getting all charges
                for k=1:length(data.formulae)
                    if k == length(data.formulae)
                        metBiGG{i,3} = strcat(metBiGG{i,3},data.formulae{k,1});
                    else
                        metBiGG{i,3} = strcat(metBiGG{i,3}, data.formulae{k,1}, ', ');
                        metBiGG{i,5} = length(data.formulae);
                        
                       if isfield(data.database_links,'MetaNetX_MNX_Chemical')
                            metaNetX = webread( data.database_links.MetaNetX_MNX_Chemical.link);
                            regex = regexp(metaNetX, 'InChI=(.+)</td></tr>   <tr><td', 'tokens');
                            if isempty(regex)~=1
                                inchi = regex{1};
                            end
                        end
                        if isempty(inchi)~=1
                             metBiGG{i,6} = inchi;
                             inchiFound = inchiFound +1;
                        end
                        
                        if  isempty(metBiGG{i,6})~=1
                            metBiGG{i,8} = getFormulaFromInChI(metBiGG{i,6});
                            formResolved = formResolved+1;
                        end
                    end
                end
            else
                
                if isfield(data.database_links,'MetaNetX_MNX_Chemical')
                    metaNetX = webread( data.database_links.MetaNetX_MNX_Chemical.link);
                    regex = regexp(metaNetX, 'InChI=(.+)</td></tr>   <tr><td', 'tokens');
                    if isempty(regex)~=1
                        inchi = regex{1};
                    end
                end
                if isempty(inchi)~=1
                     metBiGG{i,6} = inchi;
                     inchiFound = inchiFound +1;
                end
                
                noFormula = noFormula+1;
                metBiGG{i,3} = 'No formula In BiGG';
                
                if  isempty(metBiGG{i,6})~=1
                    metBiGG{i,8} = getFormulaFromInChI(metBiGG{i,6});
                    noformResolved = noformResolved+1;
                end
            end
            %checking if display flag is set to ON
            if disp == 1
                fprintf('\n%d. Looking for Formula and Charge for %s: ',i,all_biggid{i,1}(1:end-3));
                if isempty(data.formulae)~=1
                    fprintf('\nFormula: ');
                    for k=1:length(data.formulae)
                        fprintf('%s ',data.formulae{k,1});
                        if k==length(data.formulae)
                            fprintf('\n');
                        end
                    end
                end
                if isempty(data.charges)~=1
                    fprintf('\nCharge(s): ');
                    for k=1:length(data.charges)
                        fprintf('%s ',num2str(data.charges(k,1)));
                         if k==length(data.charges)
                            fprintf('\n');
                        end
                    end
                end

            end
        
  

    %case where webread returns error    
    catch
        if disp == 1
            notInbigg = notInbigg+1;
            metBiGG{i,2}= 'Not found';
            metBiGG{i,3}= 'Not found';
            fprintf('\n%d. Looking for Formula and Charge for %s: ',i,all_biggid{i,1}(1:end-3));
            fprintf('\n%s.', 'Not found')
        end
    end
if i == length(all_biggid)   
fprintf('\nSummary:\n');
fprintf('%d Metabolites were not in the BiGG database\n',notInbigg);
fprintf('%d Metabolites did not have charge specified in the BiGG database (%d Resolved using InChi)\n',noCharge, nochResolved);
fprintf('%d Metabolites did not have a formula specified in the BiGG database (%d Resolved using InChi)\n',noFormula, noformResolved);
fprintf('%d Charges Resolved using InChi ID\n',chResolved);
fprintf('%d Formulas Resolved using InChi ID\n',formResolved);
end
end
toc;