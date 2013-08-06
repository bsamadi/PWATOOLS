
function RX=pntrgn(pwasys)
% This function detects the edge points in each region from given region
% eqations E and e
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

X=computex(pwasys);
E = pwasys.E;
e = pwasys.e;

pnt_num=size(X, 2);  % number of the points
NR=pwasys.NR;
rln_num=size(E{1}, 1);

for i=1:NR
    Z=[];
    for j=1:pnt_num
        Y=E{i}*X(:,j)+e{i};
        Y=round(Y*10^10)./10^10;
        if min(sign(Y))>=0
            Z=[Z X(:,j)];
        end
    end
    RX{i}=Z;
end

    
