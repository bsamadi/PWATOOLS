%% Local function 4: Checkxcl
function er=Checkxcl(model)
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

er=0;
[nA,mA]=size(model.A);
[nx, mx]=size(model.xcl);

if mx~=1
    fprintf('xcl should be a vector of size %i \times 1.\n', nA);
    er=1;
end
if  nx~=nA
    fprintf('sizes of A and xcl do not match.\n');
    er=1;
end
end


