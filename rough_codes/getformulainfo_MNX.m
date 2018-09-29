function formula = getformulainfo_MNX(link)

options = weboptions('time',50);
data = webread(link,options);
if ~isempty(strfind(data,'formula'))
    data = data(strfind(data,'formula'):end);
    ix = strfind(data,'</td');
    formula = data(ix(1)+9:ix(2)-1);
    formula = regexprep(regexprep(formula,'</sub>',''),'<sub>','');
else
    formula = ' ';
end