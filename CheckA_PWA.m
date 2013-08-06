%% Local function 1: CheckA

function er=CheckA_PWA(model);
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

er=0;
L=model.NR;
n=model.dim(1);
size_bank=[];

for i=1:L
    [nA,mA]=size(model.A{i});
    size_bank=[size_bank; [nA mA]];
    if nA~=mA
        fprintf('Mtrices A should be square.\n\n');
        er=1;
    end
end

V=reshape(size_bank, 2*L, 1);
[I, J]=find(V);
V=V(I);   % all nonzero elements
if any(V-n)
    fprintf('All non-zero matrices A should have dimension %i x %i.\n\n', n,n);
    er=1;
end