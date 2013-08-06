%% Local function 5: CheckDomain

function er=CheckDomain(model)
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

er=0;
[nA,mA]=size(model.A);
[nD,mD]=size(model.Domain);
if max(nD, mD)~=mA
    fprintf('Domain of all variables should be defined.\n');
    er=1;
end
end

