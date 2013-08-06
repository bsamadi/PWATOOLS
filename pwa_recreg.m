function R = pwa_recreg(X,GridInfo)

% Input:
% GridInfo.Domain   Domain of interest
% GridInfo.GR       Grid resolution
% X                 Data points
%
% Output
% R                 Regions corresponding to X 

Domain = GridInfo.Domain;
GR = GridInfo.GR;

n = length(Domain);

if size(X,1)==n && size(X,2)~=n,
    X = X';
elseif size(X,2)~=n,
    if n>1,
        error(['X should have ' num2str(n) ' columns.']);
    else
        error(['X should have ' num2str(n) ' column.']);
    end
end

p = size(X,1);                      % Number of data points
Ind = pwa_index(GR);
NR = size(Ind,2);                   % Number of regions

GridPoints = cell(n,1);
for i=1:n,
    GridPoints{i}=linspace(Domain{i}(1),Domain{i}(end),GR(i)+1);
end

for i=1:p,
    R{i} = [];
    for r=1:NR,
        GP = Ind(:,r);
        Yes = zeros(n,1);
        for j=1:n,
            Xjbounds = GridPoints{j}([GP(j),GP(j)+1]);
            if Xjbounds(1)<= X(i,j) && X(i,j) <= Xjbounds(2),
                Yes(j) = 1;
            end
        end
        if all(Yes),
            R{i} = [R{i}; r];
        end
    end
end