function PWABounds = pwa_bounds_single(nlfun, Option)

% Nonlinear function: nlfun.Handle
% Vertices: nlfun.X
% Domain: nlfun.Domain
% Resolution: nlfun.Resolution
% Parameters: nlfun.Param
% The linearization point: nlfun.xstar
% The linear approximation: nlfun.AbarLin

if nargin == 1,
    Option.ObjFun = 'L2';
    Option.Type = 'AP'; % 'LB','UB','AP'
    Option.SafetyBound = 0.0; % 
else
    if ~isfield(Option,'ObjFun'),
        Option.ObjFun = 'L2';
    end
    if ~isfield(Option,'Type'),
        Option.Type = 'AP'; % 'LB','UB','AP'
    end
    if ~isfield(Option,'SafetyBound'),
        Option.SafetyBound = 0.0;
    end
end

f = nlfun.Handle;                     % Specify the handle of the nonlinear function
a=[];
b=[];
for i=1:length(nlfun.Domain),
    a=[a; nlfun.Domain{i}(1)];        % Set start point of the domain
    b=[b; nlfun.Domain{i}(end)];      % Set start point of the domain
end
N = nlfun.Resolution;                 % The resolution of sampling
X = nlfun.X;                          % Vertices of the partition

if isfield(nlfun,'Param'),
    Param = nlfun.Param;                  % Parameters of the nonlinear function
elseif isfield(nlfun,'Parameters'),
    Param = nlfun.Parameters;                  % Parameters of the nonlinear function
else
    Param = [];
end

[m n] = size(X);

if n == 1,
    X = sort(X);
end

for i = 1:n,
    if a(i) >= b(i),                  % Check ?
        eval(['error(''nlfun.domain(1,' num2str(i) ') should be less than nlfun.domain(2,' num2str(i) ').'');'])
    end
end

try
    Y = pwa_nlfun(nlfun,a,Param);
catch
    error(['The function can not be evaluated.']);
end

% Create fine mesh
W = pwa_grid(a,b,N);
Z = pwa_nlfun(nlfun,W,Param);


p = size(W,1);
q = size(Z,2);

Ver = ver('MATLAB');
if str2num(Ver.Version(1)) == 7,
    T = delaunayn(X,{'Qt','Qbb','Qc','Qz'});
else
    T = delaunayn(X);
end

numt = size(T,1);                       % Determine number of triangulations (cells)

constraints = set([]);
for i = 1:numt,
    Abar{i} = sdpvar(q,n+1,'full');     % PWA approximation parameters
end

% The approximation is exact at xstar
if isfield(nlfun,'xstar'),
%    if n == 2,
%        t = tsearch(X(:,1),X(:,2),T,nlfun.xstar(1),nlfun.xstar(2));
%    else
        t = tsearchn(X,T,nlfun.xstar');
%    end
    constraints = constraints + set(Abar{t}*[nlfun.xstar;1] == pwa_nlfun(nlfun,nlfun.xstar',Param)',['Abar' num2str(t) '=AbarLin']);
end

% Continuity constraints
if numt > 1,
    for i = 1:m,
        Ti = find(any((T==i)'));
        if length(Ti) > 1,
            for j = 1:length(Ti)-1,
                k = j+1;
                constraints = constraints + set((Abar{Ti(k)}-Abar{Ti(j)})*[X(i,:) 1]'==0,['(Abar' num2str(Ti(k)) '-Abar' num2str(Ti(j)) ')*[X(' num2str(i) ',:) 1]''=0']);
            end
        end
    end
end

% Error
% if n == 2,
%    t = tsearch(X(:,1),X(:,2),T,W(:,1),W(:,2));
% else
    t = tsearchn(X,T,W);
% end

% Objective
if ~isfield(Option,'ObjFun'),
    Option.ObjFun = 'L2';
end

Err=[];
if strcmp(Option.ObjFun,'L2'),
    Obj = 0;
    for i = 1:p,
        if strcmp(Option.Type,'LB'),
            if W(i,:)>nlfun.xstar,
                constraints = constraints + set(Z(i,:)'-Option.SafetyBound*sat(abs(W(i,:)-nlfun.xstar))>Abar{t(i)}*[W(i,:) 1]','Z>AbarWbar');
            elseif W(i,:)<nlfun.xstar,
                constraints = constraints + set(Z(i,:)'+Option.SafetyBound*sat(abs(W(i,:)-nlfun.xstar))<Abar{t(i)}*[W(i,:) 1]','Z<AbarWbar');                
            end
        elseif strcmp(Option.Type,'UB'),
            if W(i,:)>nlfun.xstar,
                 constraints = constraints + set(Z(i,:)'+Option.SafetyBound*sat(abs(W(i,:)-nlfun.xstar))<Abar{t(i)}*[W(i,:) 1]','Z<AbarWbar');
            elseif W(i,:)<nlfun.xstar,
                 constraints = constraints + set(Z(i,:)'-Option.SafetyBound*sat(abs(W(i,:)-nlfun.xstar))>Abar{t(i)}*[W(i,:) 1]','Z>AbarWbar');
            end
        end
        Error = Z(i,:)'-Abar{t(i)}*[W(i,:) 1]';
        Obj = Obj + sum(sum(Error.*Error));
        Err = [Err; Error'];
    end
end

% Solve the optimization problem
solvesdp(constraints,Obj,sdpsettings('usex0',1,'shift',1e-7)); %,sdpsettings('solver','SEDUMI'));
%solution = solvesdp(constraints,Obj,sdpsettings('solver','SEDUMI'));
checkset(constraints);

for i = 1:length(Abar),
    Abar{i} = double(Abar{i});
end
Obj = double(Obj);
Err = double(Err);

Y = [];
for i = 1:m,
    J = find(any((T==i)'));
    Y = [Y; (Abar{J(1)}*[X(i,:) 1]')'];
end

PWABounds.Abar = Abar;
PWABounds.Obj = Obj;
PWABounds.X = X;
PWABounds.Y = Y;
PWABounds.W = W;
PWABounds.Z = Z;
PWABounds.T = T;
PWABounds.Err = Err;

function y = sat(x)
if abs(x)<1,
    y=x;
else,
    y = sign(x);
end