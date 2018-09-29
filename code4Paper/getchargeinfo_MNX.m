function charge = getchargeinfo_MNX(link)
% gets charge info from MetaNetX
% link is webaddress for metabolite page on MetaNetX
% written by chintan joshi on 1/25/2017 at ucsd

options = weboptions('time',60,'CertificateFilename','');
data = webread(link,options);
if ~isempty(strfind(data,'charge'))
    data = data(strfind(data,'charge'):end);
    ix = strfind(data,'</td');
    charge = str2double(data(ix(1)+9:ix(2)-1));
else
    charge = NaN;
end