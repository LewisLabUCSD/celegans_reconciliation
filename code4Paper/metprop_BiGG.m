function prop = metprop_BiGG(all_bigg,info_type)
% prop = metprop_BiGG(all_bigg,info_type)
% finds a given info_type or all into_types
%
% INPUT:
% all_bigg: BiGG ids
% info_type: type of information needed
% (values = 'name','charge','formula','BiGG','KEGG')
% written by chintan joshi on 2/8/2017 at ucsd
if nargin < 2
    info_type = 'all';
end
options = weboptions('time',10,'RequestMethod','GET');
all_bigg = regexprep(all_bigg,'\[c\]|\[e\]|\[m\]|\[n\]','');
for i=1:length(all_bigg)
%     fprintf('%d. Finding in BiGG: %s...',i,all_bigg{i,1});
    try
        data = webread(strcat('http://bigg.ucsd.edu/api/v2/universal/metabolites/',all_bigg{i,1}),options);
%         fprintf('Found\n');
    catch
%         fprintf('Not Found\n');
        data.name = 'not found';
        data.charges = NaN;
        data.formulae = {'not found'};
        data.database_links.KEGGCompound = 'not found';
        data.bigg_id = 'not found';
        data.source = 'not found';
    end
    if strcmpi(info_type,'all') || strcmpi(info_type,'name')
        name{i,1} = data.name;
    end
    if strcmpi(info_type,'all') || strcmpi(info_type,'charge')
        if isempty(data.charges)
%             fprintf('Charge not found in BiGG, finding in MetaNetX...');
            charges{i,1} = getchargeinfo_MNX(data.database_links.MetaNetX_MNX_Chemical.link);
%             fprintf('Found\n');
        elseif length(data.charges) > 1
%             fprintf('More than one charges found in BiGG, finding in MetaNetX...');
data.database_links.MetaNetX_MNX_Chemical.link
            charges{i,1} = getchargeinfo_MNX(data.database_links.MetaNetX_MNX_Chemical.link);
%             fprintf('Found\n');
        else
            charges{i,1} = data.charges;
        end
    end
    if strcmpi(info_type,'all') || strcmpi(info_type,'formula')
        if isempty(data.formulae)
%             fprintf('Formula not found in BiGG, finding in MetaNetX...');
            formula{i,1} = getformulainfo_MNX(data.database_links.MetaNetX_MNX_Chemical.link);
%             fprintf('Found\n');
        elseif size(data.formulae,1) > 1
%             fprintf('More than one formulas found in BiGG, finding in MetaNetX...');
            formula{i,1} = getformulainfo_MNX(data.database_links.MetaNetX_MNX_Chemical.link);
%             fprintf('Found\n');
        else
            formula{i,1} = data.formulae{1,1};
        end
    end
    if strcmpi(info_type,'KEGG')
        if isfield(data.database_links,'KEGGCompound')
            eval(strcat(info_type,'{i,1} = data.database_links.KEGGCompound;'));
        else
            eval(strcat(info_type,'{i,1} = ''not found'';'));
        end
    elseif strcmpi(info_type,'all')
        if isfield(data.database_links,'KEGGCompound')
            KEGG{i,1} = data.database_links.KEGGCompound;
        else
            KEGG{i,1} = 'not found';
        end
    end
    if strcmpi(info_type,'BiGG')
        eval(strcat(info_type,'{i,1} = data.bigg_id;'));
    elseif strcmpi(info_type,'all')
        BiGG{i,1} = data.bigg_id;
    end
end
if strcmpi(info_type,'all')
    prop.names = name;
    prop.charge = charges;
    prop.formula = formula;
    prop.KEGGID = KEGG;
    prop.BiGG = BiGG;
    
%     fprintf('%d metabolites have more than one possible charge info.\n',length(find(cellfun(@length,prop.charge)==2)));
%     fprintf('%d metabolites have more than one possible formula info.\n',length(find(cellfun(@length,prop.formula)==2)));
%     K = length(find(strcmp(prop.BiGG,'not found')));
%     fprintf('BiGG ids could not be found for %d metabolites.\n',K);
%     L = length(find(strcmp(prop.KEGGID(~strcmp(prop.BiGG,'not found')),'not found')));
%     fprintf('KEGG ids could not be found for %d out of %d BiGG ids found.\n',L,length(all_bigg)-K);
else
    eval(strcat('prop.',info_type,'=',info_type));
end