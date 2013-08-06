
function [error, index]=IsC1Function(pwasys)
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

error=0;
index=[];
Tol=1e-2;
RX=pntrgn(pwasys);
[NR NS] = size(pwasys.Abar);

A = pwasys.A;
a = pwasys.a;
B = pwasys.B;
E = pwasys.E;
e = pwasys.e;
F = pwasys.F;
f = pwasys.f;
pwatype=pwasys.type;
if strcmp(pwatype, 'lower-envelope')
    col_index=[1];
elseif strcmp(pwatype, 'pwadi')
    col_index=[1 2];
elseif strcmp(pwatype, 'null')
    col_index=[];
end
n=size(A{1},1);

for i=1:NR
    for j=i+1:NR
        if ~isempty(F{i,j})
            for h=1:3
                x=F{i,j}*rand(n-1, 1)+f{i,j};
                for k=1:col_index
                    Z_i=A{i,k}*x+a{i,k};
                    Z_j=A{j,k}*x+a{j,k};
                    if norm(Z_i-Z_j)>Tol
                        index=[index; i j];
                        error=1;
                        return
                    end
                end
            end
        end
    end
end
