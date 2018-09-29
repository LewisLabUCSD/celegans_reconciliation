function formula = getformulainfo_MNX(link)
% gets formula info for a metabolite usingg MetNetX webaddress (link)
% written by chintan joshi on 2/8/2017 at ucsd

options = weboptions('time',50,'CertificateFilename','');
% link = strrep(link,'http','https')
data = webread(link,options);
if ~isempty(strfind(data,'formula'))
    data = data(strfind(data,'formula'):end);
    ix = strfind(data,'</td');
    formula = data(ix(1)+9:ix(2)-1);
    formula = regexprep(regexprep(formula,'</sub>',''),'<sub>','');
else
    formula = ' ';
end