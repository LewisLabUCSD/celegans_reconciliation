function [genes,gene_entry,twice_present,header] = getgeneinfo_WormBase(filename_local,gene_entry,disp_flag)
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
    genes_col2 = genes;
    genes_col3 = genes;
    rem_genes = {'NA';'ND';'Unknown';'TBD'};
    if sum(ismember(gene_entry,rem_genes))~=0
        gene_entry(ismember(gene_entry,rem_genes)) = [];
    end
    col2 = repmat({'public_name (WB)'},sum(ismember(genes(:,2),gene_entry)),1); % whether the gene in the list is recognized by public_name
    genes_col2(~ismember(genes(:,2),gene_entry),:) = [];
    gene_entry(ismember(gene_entry,genes(:,2))) = [];
    col3 = repmat({'molecular_name (WB)'},sum(ismember(genes(:,3),gene_entry)),1); % whether the gene in the list is recognized by molecular_name
    genes_col3(~ismember(genes(:,3),gene_entry),:) = [];
    twice_present = genes_col2(ismember(genes_col2(:,2),genes_col3(:,2)),:);
    col2(ismember(genes_col2(:,2),genes_col3(:,2))) = [];
    genes_col2(ismember(genes_col2(:,2),genes_col3(:,2)),:) = [];
    genes = [genes_col2;genes_col3];
    cols = [col2;col3];
    gene_entry(ismember(gene_entry,genes(:,3))) = [];
    genes(:,5) = cols;
    header{1,5} = 'info_type_in_model';
    twice_present(:,[1,4]) = [];
%     genes(:,[1 4 5]) = [];
else
    genes(:,5:7) = []; % remove unwanted columns (is likely to change depending upon the file
    header(5:7) = [];
    gene_entry = [];
    twice_present = [];    
end

if disp_flag
    fprintf('%d genes found in WormBase.\n',size(genes,1));
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