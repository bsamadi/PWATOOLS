function CellInfo = pwa_cell(pwafun)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Create the cell descriptions and boundary descriptions using          %
%   information about the vertices and the triangulation.                 %
%                                                                         %
%   Input:                                                                %
%       pwafun.X: Vertices                                                %
%       pwafun.T: Triangulation                                           %
%                                                                         %
%   Output:                                                               %
%       CellInfo.Ebar: Cell descriptions                                  %
%       CellInfo.Fbar: Boundary descriptions                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Obtain dimensions of Triangulation
[m n] = size(pwafun.T);
% Create cell descriptions
for i = 1:m,
    CellInfo.Ebar{i} = [];
    for j = 1:n,
        Index = 1:n;
        Index(j)=[];
%         M = [pwafun.X(pwafun.T(i,Index),:) ones(n-1,1)];
%         Ebar = [1 -(M(:,2:end)\M(:,1))'];
        M = [pwafun.X(pwafun.T(i,Index),:) ones(n-1,1)];
        Ebar = null(M)';
        xbar = [pwafun.X(pwafun.T(i,j),:) 1]';
        if all(Ebar*xbar < 0),
            Ebar = -Ebar;
        end
        CellInfo.Ebar{i} = [CellInfo.Ebar{i}; Ebar];
    end
end

% Create boundary descriptions
for i = 1:m-1,
    for j = i+1:m
        Index = intersect(pwafun.T(i,:),pwafun.T(j,:));
        if ~isempty(Index) && length(Index)==n-1,
%             M = [pwafun.X(Index,:) ones(n-1,1)];
%             if rank(M(:,2:end))==length(M(:,2:end)),
%                 CellInfo.Fbar{i,j} = [(M(:,2:end)\M(:,1))';eye(n-1)];
%             else
%                 error('Computation of Fbar should be updated!');
%             end
            M = [pwafun.X(Index,:) ones(n-1,1)]';
            CellInfo.Fbar{i,j} = M*[eye(n-2,n-1);-ones(1,n-2) 1];
            CellInfo.Fbar{j,i} = CellInfo.Fbar{i,j};
        end
    end
end