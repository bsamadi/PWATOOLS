%% Local function 3: CheckBx

function er=CheckB_PWA(model)
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

er=0;
L=model.NR;
n=model.dim(1);
m=model.dim(2);

for i=1:L
    [nB,mB]=size(model.B{i});
    if (nB~=n)
        fprintf('Sizes of A and B do not match.\n');
        er=1;
        return;
    end
    if (mB~=m)
        fprintf('Sizes of B do not match the number of inputs.\n');
        er=1;
        return;
    end
    
end





