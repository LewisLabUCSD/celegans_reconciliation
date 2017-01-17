function [commonGeneRxnsModel1, commonGeneRxnsModel2, rxnsModel1, rxnsModel2] = rxnsFromGenes(model1, model2)
    genesModel1Not2 = setdiff(model1.genes, model2.genes);
    genesModel2not1 = setdiff(model2.genes, model1.genes);
    genesPresentInBoth = intersect(model1.genes, model2.genes);

    %reactions assosiated with same genes Kaleta
    rxnsBoth2S = findRxnsFromGenes(model2, genesPresentInBoth);

    %reactions assosiated with same genes Kaleta
    rxnsBoth1S = findRxnsFromGenes(model1, genesPresentInBoth);

    %reactions unique to each model according to their genes
    rxns1 = findRxnsFromGenes(model1, genes1Not2);


    rxns2 = findRxnsFromGenes(model2, genes2Not1);

    commonGeneRxnsModel2 = fieldnames(rxnsBoth2S)
    for  i = 1:numel(commonGeneRxnsModel2)
        gene = rxnsBoth2S.(commonGeneRxnsModel2{i});
        str='';
        for j = 1:size(gene,1)
            rxns(1,j) = {gene{j,1}};
        end 
        rxns = [rxns{:}];
        str = strjoin(rxns, ', ');
        commonGeneRxnsModel2(i,2) = {str};
        rxns(1,:)=[];
    end

    commonGeneRxnsModel1 = fieldnames(rxnsBoth1S)
    for  i = 1:numel(commonGeneRxnsModel1)
        gene = rxnsBoth1S.(commonGeneRxnsModel1{i});
        str='';
        for j = 1:size(gene,1)
            rxns(1,j) = {gene{j,1}};
        end 
        rxns = [rxns{:}];
        str = strjoin(rxns, ', ');
        commonGeneRxnsModel1(i,2) = {str};
        rxns(1,:)=[];
    end

    rxnsModel1 = fieldnames(rxns1)
    for  i = 1:numel(rxnsModel1)
        gene = rxns1.(rxnsModel1{i});
        str='';
        for j = 1:size(gene,1)
            rxns(1,j) = {gene{j,1}};
        end 
        rxns = [rxns{:}];
        str = strjoin(rxns, ', ');
        rxnsModel1(i,2) = {str};
        rxns(1,:)=[];
    end

    rxnsModel2 = fieldnames(rxns2)
    for  i = 1:numel(rxnsModel2)
        gene = rxns2.(rxnsModel2{i});
        str='';
        for j = 1:size(gene,1)
            rxns(1,j) = {gene{j,1}};
        end 
        rxns = [rxns{:}];
        str = strjoin(rxns, ', ');
        rxnsModel1(i,2) = {str};
        rxns(1,:)=[];
    end
end



