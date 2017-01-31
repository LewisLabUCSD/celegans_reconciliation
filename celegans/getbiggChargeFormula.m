function metBiGG = getbiggChargeFormula(all_bigg, disp)
% metBiGG = getbiggChargeFormula(all_bigg, 1)
% INPUT: 
% all_bigg: list of Bigg compound ids
% disp: flag - 1 to show outpput - 0 - not to show 
% OUTPUT:
% metBiGG: list of Chemical Formulas and Charges from
%          Bigg 
tic;
metBiGG = cell(length(all_bigg),5);
options = weboptions('Timeout',30, 'RequestMethod', 'GET');
prev='';
notInbigg = 0;
noCharge=0;
noFormula=0;
%check all Bigg IDs
for i=1:length(all_bigg)
    metBiGG{i,1}= all_bigg{i,1};
    try
        %if not same id as previous one
        if strcmp(prev, all_bigg{i,1}(1:end-3))~=1
           
            prev = all_bigg{i,1}(1:end-3);
            data= webread(strcat('http://bigg.ucsd.edu/api/v2/universal/metabolites/', all_bigg{i,1}(1:end-3)),options);
            
            %if there is a charge in BIGG
            if isempty(data.charges)~=1
                
                %getting all charges
                for k=1:length(data.charges)
                    if k == length(data.charges)
                        metBiGG{i,2} = strcat(metBiGG{i,2}, num2str(data.charges(k,1)));
                    else
                        metBiGG{i,2} = strcat(metBiGG{i,2}, num2str(data.charges(k,1)), ', ');
                        metBiGG{i,4} = length(data.charges);
                    end
                end 
            else
                noCharge = noCharge+1;
                metBiGG{i,2} = 'No Charge In BiGG';
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
                    end
                end
            else
                noFormula = noFormula+1;
                metBiGG{i,3} = 'No formula In BiGG';
            end
            %checking if display flag is set to ON
            if disp == 1
                fprintf('\n%d. Looking for Formula and Charge for %s: ',i,all_bigg{i,1}(1:end-3));
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
        
        %ID is same as previous one    
        else
            metBiGG{i,2}= metBiGG{i-1,2};
            metBiGG{i,3}= metBiGG{i-1,3}; 
        end
    %case where webread returns error    
    catch
        if disp == 1
            notInbigg = notInbigg+1;
            metBiGG{i,2}= 'Not found';
            metBiGG{i,3}= 'Not found';
            fprintf('\n%d. Looking for Formula and Charge for %s: ',i,all_bigg{i,1}(1:end-3));
            fprintf('\n%s.', 'Not found')
        end
    end
if i == length(all_bigg)   
fprintf('\nSummary:\n');
fprintf('%d Metabolites were not in the BiGG database\n',notInbigg);
fprintf('%d Metabolites did not have charge specified in the BiGG database\n',noCharge);
fprintf('%d Metabolites did not have a formula specified in the BiGG database\n',noFormula);
end
end
toc;