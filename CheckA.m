%% Local function 1: CheckA

function er=CheckA(model);
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%
er=0;
[nA,mA]=size(model.A);
if nA~=mA
    fprintf('Matrix A should be square.\n\n');
    er=1;
end
end


