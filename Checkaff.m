%% Local function 6: Checkxcl
function er=Checkaff(model)
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

er=0;
[nA,mA]=size(model.A);
[na, ma]=size(model.aff);

if (ma~=1 | na~=nA)
    fprintf('Affine term should be a vector of size %i \times 1.\n', nA);
    er=1;
end

end

