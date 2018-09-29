function [genes,gene_entry,header] = getgeneinfo_WormBase(filename_local,gene_entry,disp_flag)
% [genes,header,gene_entry] = getallgenes_WormBase(filename,gene_entry)
% gets information on all the genes querried in gene_entry, else all genes
% in WormBase
%
% INPUT:
% filename: provide the location of the file on disk, or 'online' (optional)
% gene_entry: list of genes that need querrying (optional), genes_not_found
%
% OUTPUT:
% genes: list of genes
% format of genes = {'gene_id','public_name','molecular_name','gene_class_description','info_type_in_model'}
% column 5 will only exist if gene_entry exists
% gene_entry: list of genes not found
% twice_present: list of genes that are duplicates
% header: format of genes
%
% Written by Chintan Joshi
if nargin < 3
    disp_flag = false;
end
if strcmp(filename_local,'online') % this if loop goes to the ftp site and downloads the file and data
    db = ftp('ftp.wormbase.org');
    cd(db,'pub/wormbase/species/c_elegans/annotation/functional_descriptions/');
    mget(db,'c_elegans.canonical_bioproject.current.functional_descriptions.txt.gz');
    close(db);
    gunzip('c_elegans.canonical_bioproject.current.functional_descriptions.txt.gz');
    filename_local = 'c_elegans.canonical_bioproject.current.functional_descriptions.txt';
end
fid = fopen(filename_local,'r'); % open the file for reading
A = textscan(fid,'%s','delimiter','\n'); % read the file
A = A{1,1};
header = A(4); % assign header
A(1:4) = []; % remove header and file description from file
genes = cell(length(A),8);
header = regexp(header,' ','split');
header = header{1,1};
for i=1:length(A) % got to each line
    tindex = regexp(A{i,1},'\t');  % identify tab-delimiters
    if ~isempty(tindex) % if found
        for j=1:(length(tindex)+1) % put each tab-delimited value into a cell
            if j==1
                genes{i,j} = strtrim(A{i,1}(1:tindex(j)));
            elseif j==length(tindex)+1
                genes{i,j} = strtrim(A{i,1}(tindex(j-1)+1:end));
            else
                genes{i,j} = strtrim(A{i,1}(tindex(j-1)+1:tindex(j)));
            end
        end
    else
        genes{i,1} = A{i,1}; % if empty, leave the gene in first column
        genes(i,2:8) = {' '};
    end
end

if nargin>=2
    genes(:,4:7) = []; % remove unwanted columns (is likely to change depending upon the file
    header(4:7) = [];
    rem_genes = {'NA';'ND';'Unknown';'TBD'};
    
    %     making sure that gene_entry has not ignored iso forms in
    %     molecular_name (WB)
    iso = find(~ismember(gene_entry,genes(:,3)));
    iso_bool = false(length(genes(:,3)),1);
    gentry_iso_bool = false(length(gene_entry),1);
    for i=1:length(iso)
        iso_bool = iso_bool | ~cellfun(@isempty,regexp(genes(:,3),strcat(gene_entry(iso(i)),'[a-z]'),'match'));
        if sum(~cellfun(@isempty,regexp(genes(:,3),strcat(gene_entry(iso(i)),'[a-z]'),'match'))) ~= 0
            gentry_iso_bool(iso(i)) = true;
%             gene_entry(iso(i)) = genes(~cellfun(@isempty,regexp(genes(:,3),strcat(gene_entry(iso(i)),'[a-z]'),'match')),3);
            genes(~cellfun(@isempty,regexp(genes(:,3),strcat(gene_entry(iso(i)),'[a-z]'),'match')),3) = gene_entry(iso(i));
        end
    end
    genes_col2 = genes;
    genes_col3 = genes;
    %     accounting perfect matches with molecular_name (WB)
    mem_bool3 = ismember(genes(:,3),gene_entry);
    mem_bool2 = ismember(genes(:,2),gene_entry);
    if sum(ismember(gene_entry,rem_genes))~=0
        gentry_iso_bool(ismember(gene_entry,rem_genes)) = [];
        gene_entry(ismember(gene_entry,rem_genes)) = [];
    end
    col2 = repmat({'public_name (WB)'},sum(mem_bool2),1); % whether the gene in the list is recognized by public_name
    genes_col2(~mem_bool2,:) = [];
%     gene_entry(ismember(gene_entry,genes(:,2))) = [];
    col3 = repmat({'molecular_name (WB)'},sum(mem_bool3 | iso_bool),1); % whether the gene in the list is recognized by molecular_name
    genes_col3(~(mem_bool3|iso_bool),:) = [];
    genes3 = genes_col3; genes2 = genes_col2;
    index = ismember(gene_entry,genes3(:,3)) | gentry_iso_bool;
    index = index | ismember(gene_entry,genes2(:,2));
    gene_entry(index) = [];
    genes3(:,5) = col3;
    genes2(:,5) = col2;
    [~,ia,ib] = intersect(genes3(:,3),genes2(:,3));
    genes2(ib,:) = [];
    genes = [genes3;genes2];
    header{1,5} = 'info_type_in_model';
    twice_present = [];
else
    genes(:,5:7) = []; % remove unwanted columns (is likely to change depending upon the file
    header(5:7) = [];
    gene_entry = [];
    twice_present = [];    
end

if disp_flag
    fprintf('%d genes found in WormBase.\n',size(genes3,1));
    fprintf('%d genes in the list have been annotated based on public_name.\n',size(col2,1));
    fprintf('%d genes in the list have been annotated based on molecular_name.\n',size(col3,1));
    fprintf('%d genes have duplicates.\n',size(twice_present,1));
    fprintf('%d genes were not found in WormBase.\n',size(gene_entry,1));
    fprintf('%d unique genes present in the list.\n',size(genes,1)+size(gene_entry,1));
end

fclose('all');
if strcmp(filename_local,'online')
    delete('c_elegans.canonical_bioproject.current.functional_descriptions.txt');
    delete('c_elegans.canonical_bioproject.current.functional_descriptions.txt.gz');
end