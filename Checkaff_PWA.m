%% Local function 6: Checkxcl
function er=Checkaff_PWA(model)
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

er=0;
L=model.NR;

n=model.dim(1);

for i=1:L
    [na,ma]=size(model.a{i});
    if all([na ma])
        if (ma~=1 | na~=n)
            fprintf('Nonzedro affine terms should be a vector of size %i \times 1.\n', n);
            er=1;
        end
    end
end
