
# coding: utf-8

# In[30]:
import xml.etree.ElementTree as ET
from string import punctuation
import re 


#pulling both model files
kaletaModel = ET.parse('C:\Users\user\Documents\Lewis_LAB\SBML files\Celegans_ecoli.xml').getroot()
icellMod = ET.parse('C:\Users\user\Documents\Lewis_LAB\SBML files\iCEL1273.xml').getroot()

root = kaletaModel.getchildren()

#pulling lists of rections for both models
kaletaListOfRxns = root[0].getchildren()[2]
icellListOfRxns = icellMod.getchildren()[0].getchildren()[4]

genesK = []
genesI = []

#creating regular expression string to match the format of the genes in both models
iGeneFormat = re.compile('(?:(?:^GENE_LIST: )|(?:; ))([^\s;]+)')
kGeneFormat = re.compile('(?:GENE ASSOCIATION: \(| or \(|GENE_ASSOCIATION: \(| and | and \(| or |GENE ASSOCIATION: |GENE_ASSOCIATION: )((\w+\.\w+)|(\w+-\w+\.\w+)|(\w+.\w+))')

#iterating through the Icel list of reactions
for rxnI in icellListOfRxns:
    if len(rxnI)>1:
        rxnNotesI = rxnI[0] #getting the reaction notes
        
        for x in range(len(rxnNotesI)): #find the genes notes
            geneList = re.findall(iGeneFormat,str(rxnNotesI[x].text))
            if geneList:
                break
        if len(geneList)>0:
            for gene in geneList:
                #add gene if not in the list
                if gene not in genesI:
                    genesI.append(gene)
                    
#write all the genes into the txt file                    
iGenesTxt = open('iGenesList.txt', 'w')
for gene in genesI:
    iGenesTxt.write(str(gene) +'\n')

#repeat the process for the kaleta model    
for rxnK in kaletaListOfRxns:
    rxnNotesK = rxnK[0][0]
    length =  len(rxnNotesK)
    for x in range(length):
        string = str(rxnNotesK[x].text)
        if re.search(r'\bgene association\b', string, re.I) or re.search(r'\bgene_association\b', string, re.I):
            #print string
            geneListK = re.findall(kGeneFormat,string)
            if geneListK:
                break
    if len(geneListK)>0:
        for gene in geneListK:
            #print str(gene)    
            
            if gene not in genesK:
                genesK.append(gene)   

kGenesTxt = open('kGenesList.txt', 'w')
for gene in genesK:
    
    kGenesTxt.write(str(gene[0]) +'\n')