%% Local function 3: CheckBx

function er=CheckK_PWA(model)
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

er=0;
L=model.NR;
n=model.dim(1);
m=model.dim(2);

if isfield(model, 'K') & isfield(model, 'k')
    for i=1:L
        [nKi,mKi]=size(model.K{i});
        [nki,mki]=size(model.k{i});
        if (nKi~=m | mKi~=n)
            fprintf('Sizes of K and the order of the system/number of inputs do not match.\n');
            er=1;
            return;
        end
        if (nki~=m)
            fprintf('Sizes of k do not match the number of the inputs.\n');
            er=1;
            return;
        end
        if (mki~=1)
            fprintf('All k should be a vector of size %i x 1 .\n', m);
            er=1;
            return;
        end
    end
    
end

end



