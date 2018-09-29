function metBiGG = getbigg_KEGGID(all_keggid)
% metBiGG = getbigg_KEGGID(all_keggid)
% INPUT: 
% all_keggid: list of KEGG compound ids
% OUTPUT:
% metBiGG: list of IDs Id's retrieved from KEGG

metBiGG = cell(length(all_keggid),1);
options = weboptions('Timeout',100);

for i=1:length(all_keggid)

   %check for proper formatting  
   if strcmp(all_keggid{i,1},' ')~=1 && length(all_keggid{i,1})==6
       if ~strcmp(all_keggid{i,1}(1:2), 'CC') 
           
            fprintf('KEGG:  %s\n', all_keggid{i,1})
            
            %read and edit kegg data
            data = webwrite('http://bigg.ucsd.edu/advanced_search_external_id_results','database_source','kegg.compound', 'query', all_keggid{i,1}, options);
            id = regexp(data, '/models/universal/metabolites/([^"]+)', 'match');
            display(length(id));
            if length(id)>30
                id = id{1}(31:length(id{1}));
            end
            metBiGG{i,1} = id;
            fprintf('BiGG:');
            display(metBiGG{i,1});
       end
   end
end
