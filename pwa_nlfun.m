function Y = pwa_nlfun(nlfun,X,Param)

% Evaluate Function
% Y = pwa_nlfun(nlfun,X,Param)

% Define variables locally
if isfield(nlfun,'Handle'),
    f = nlfun.Handle;
elseif isfield(nlfun,'NonlinearFunction'),
    f = nlfun.NonlinearFunction;
end

if isfield(nlfun,'Dim'),
    Dim = nlfun.Dim;
    Index = nlfun.Index;
else
    Dim = size(X,2);
    Index = 1:Dim;
end

%
if nargin == 2,
    Param = [];
end

m = size(X,1);

%
Y = [];
for i = 1:m,
    %Xi = zeros(1,Dim);
    try,
        Xi = nlfun.xcl;
    catch,
        Xi = nlfun.xstar;        
    end
    Xi(Index) = X(i,:);

    if isempty(Param),
        Yi = feval(f,Xi);
    else
        Yi = feval(f,Xi,Param);
    end

    if isfield(nlfun,'YIndex'),
        Yi = Yi(nlfun.YIndex);
    end
    Y = [Y; reshape(Yi,1,prod(size(Yi)))];
end
