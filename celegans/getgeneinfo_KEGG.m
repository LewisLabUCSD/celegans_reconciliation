function [genes,genes_not_found] = getgeneinfo_KEGG(org_code,gene_entry,web_options,disp_flag)
% genes = getgeneinfo_KEGG(org_code,gene_entry,web_options)
% gets all the info on all genes in gene_entry list from KEGG
% user-defined functions used:
% 1. getallgenes_KEGG
%
%
% INPUT:
% org_code: 3-4 letter organism code
% gene_entry: list of genes to be querried
% web_options: options for webread
% 
% OUTPUT:
% genes: list of genes with following format
% column1-gene entry, column2-KO, column3-gene name, and
% column4-definition
% 
% Written by Chintan Joshi

if nargin < 3
    web_options = weboptions('time',10);
end
if nargin < 4
    disp_flag = false;
end
rem_genes = {'NA';'ND';'Unknown';'TBD'};
gene_entry(ismember(gene_entry,rem_genes)) = [];
genes = cell(length(gene_entry),4);
h = waitbar(0,'Finding genes in KEGG...');
for i=1:length(gene_entry)
    A = webread(strcat('http://rest.kegg.jp/find/',org_code,'/',gene_entry{i,1}),web_options);
    A = getallgenes_KEGG(org_code,false,A);
    A = A(strcmp(A(:,1),gene_entry{i,1}),:);
    if isempty(A)
        genes(i,:) = repmat({' '},1,4);
    else
        genes(i,:) = A;
    end
    if disp_flag
        fprintf('%d.Query: %s....\n',i,gene_entry{i,1});
    end
    waitbar(i/length(gene_entry));
end
close(h);

genes(strcmp(genes(:,1),' '),:)=[];
genes_not_found = setdiff(gene_entry,genes(:,1));
if disp_flag
    fprintf('%d genes found in KEGG.\n',size(genes,1));
    fprintf('%d genes not found in KEGG.\n',size(genes_not_found,1));
end