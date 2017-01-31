function charge = getchargeinfo_MNX(link)

options = weboptions('time',50);
data = webread(link,options);
if ~isempty(strfind(data,'charge'))
    data = data(strfind(data,'charge'):end);
    ix = strfind(data,'</td');
    charge = str2double(data(ix(1)+9:ix(2)-1));
else
    charge = NaN;
end