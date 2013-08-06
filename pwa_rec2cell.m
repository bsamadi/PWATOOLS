function CellInfo = pwa_rec2cell(GridInfo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Create the cell descriptions and boundary descriptions using          %
%   information about the rectangular grid.                               %
%                                                                         %
%   Input:                                                                %
%       GridInfo.Domain: Domain                                           %
%       GridInfo.GR: Grid resolution                                      %
%                                                                         %
%   Output:                                                               %
%       CellInfo.Ebar: Cell descriptions                                  %
%       CellInfo.Fbar: Boundary descriptions                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GR = GridInfo.GR;                   % Grid resolution
Domain = GridInfo.Domain;           % Domain

n = length(Domain);                 % Number of variables

Ind = pwa_index(GR);                % Indice of the regions
NR = size(Ind,2);                   % Number of regions

for i=1:n,
    GridPoints{i}=linspace(Domain{i}(1),Domain{i}(end),GR(i)+1);
end

Ebar = cell(NR,1);
for i=1:NR,
    GP = Ind(:,i);
    for j=1:n,
        % GridPoints{j}(GP(j))<= x(j) <= GridPoints{j}(GP(j)+1)
        Ebar{i} = [Ebar{i}; zeros(1,j-1) 1 zeros(1,n-j) -GridPoints{j}(GP(j))];
        Ebar{i} = [Ebar{i}; zeros(1,j-1) -1 zeros(1,n-j) GridPoints{j}(GP(j)+1)];
    end
end

CellInfo.Ebar = Ebar;

Fbar = cell(NR,NR);

xb = nan*ones(n,1);

% Create boundary descriptions
for i = 1:NR-1,     % Region i
    GPi = Ind(:,i);
    for j = i+1:NR  % Region j
        GPj = Ind(:,j);
        for k=1:n,
            GP = intersect(GPi(k),GPj(k));
            if ~isempty(GP),
                xb(k) = GridPoints{k}(GP);
            end
        end
        db = sum(isnan(xb));
        if db<n,
            ind = find(isnan(xb));
            Fbar{i,j} = zeros(n+1,length(ind)+1);
            Fbar{i,j}(ind,end) = xb(ind);
            notind = sort(setdiff(1:n,ind));
            Fbar{i,j}(notind,notind) = eye(length(notind));
        end
    end
end

CellInfo.Fbar = Fbar;