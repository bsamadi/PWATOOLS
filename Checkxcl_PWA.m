%% Local function 4: Checkxcl
function er=Checkxcl_PWA(model)
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

er=0;
n=model.dim(1);

[nx, mx]=size(model.xcl);

if mx~=1
    fprintf('xcl should be a vector of size %i \times 1.\n', nA);
    er=1;
end
if  nx~=n
    fprintf('sizes of A and xcl do not match.\n');
    er=1;
end



