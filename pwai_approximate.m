function pwadi = pwadi_approximate(nlsys)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Create a piecewise-affine (PWA) differential inclusions pwadi from    %
%   the nonlinear system nlsys. First, the nonlinear function nlfun is    %
%   extracted from nlsys. Then, nlfun is bounded by PWA functions. And    %
%   finally, the overall PWA differential inclusions pwadi is created.    %
%                                                                         %
%   Input:                                                                %
%       nlsys                                                             %
%                                                                         %
%   Output:                                                               %
%       pwadi                                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if  ~isfield(nlsys,'NonlinearFunction'),
    nlsys.NonlinearFunction = inline('0*x;');
end

% Specify the handle of the nonlinear function
nlfun.Handle = nlsys.NonlinearFunction;

% Define the domain of nonlinearity and an index of variables in the domain
nlfun.Domain = [];
nlfun.Index = [];
Domain = nlsys.Domain;

if  isfield(nlsys,'NonlinearDomain'),
    index = find(nlsys.NonlinearDomain==0);
    for i = index',
        nlsys.Domain{i}=[];
    end
end
j=1;
for i = 1:length(nlsys.Domain),
    if ~isempty(nlsys.Domain{i}),
        nlfun.Index = [nlfun.Index i];
        nlfun.Domain{j} = nlsys.Domain{i};
        j=j+1;
    end
end

if ~isfield(nlsys,'xcl'),
    error('Closed loop equilibrium point (nlsys.xcl) is not defined.');
else
    nlfun.xcl = nlsys.xcl(nlfun.Index);
end

% Set the dimension of the nonlinear function
nlfun.Dim = length(nlfun.Domain);
% Obtain indices of state variables described by nonlinear functions
YIndex = find(nlsys.NonlinearStateEquations);
% Obtain indices of state variables present in domain of nonlinearity
XIndex = nlfun.Index;
nlfun.YIndex = YIndex;

% Set the system parameters
if isfield(nlsys,'Parameters')
    nlfun.Parameters = nlsys.Parameters;
end
% Set the sampling resolution
nlfun.Resolution = nlsys.Resolution;
% Set the objective (L2 or Linf)
if isfield(nlsys,'ObjFun')
    nlfun.ObjFun = nlsys.ObjFun;
end

if ~isfield(nlsys,'Method'),
    nlsys.Method = 'Rectangular';
end

% Set the resolution of the uniform gird
if isfield(nlsys,'GR'),
    nlfun.GR = nlsys.GR;
else
    nlfun.GR = 1;
end

% Determine approximation method to use and compute the PWA approximation
if strcmpi(nlsys.Method,'Rectangular'),
    % Set the resolution of the gird
    % Compute the approximation
    pwafun = pwai_rectangular(nlfun);
end

% Obtain the cell information
CellInfo = pwa_rec2cell(pwafun);

% Define the number of state variables
n = size(nlsys.A,1);

% Determine if B(x) is a matrix
try
    nlsys.B = 0*nlsys.B+nlsys.B;
    [mB nB] = size(nlsys.B);
    ismatrix = 1;
catch
    ismatrix = 0;
end

% Convert the approximation to the system form

% Obtain indices of state variables not in the domain of nonlinearity
NotXIndex = setdiff(1:n,XIndex);

NR = prod(nlfun.GR);             % Number of regions
nE = size(CellInfo.Ebar{1},1);

for i = 1:NR,
    for j = 1:2,
        if isempty(YIndex),
            pwasys.Abar{i,j} = [nlsys.A nlsys.a; zeros(1,n+1)];
        else
            % Create the auxiliary matrices for A and a
            Afun = zeros(n+1,n+1);
            Afun(YIndex,XIndex) = pwafun.Abar{i,j}(:,1:nlfun.Dim);
            Afun(YIndex,end) = pwafun.Abar{i,j}(:,end);
            pwasys.Abar{i,j} = [nlsys.A nlsys.a; zeros(1,n+1)]+Afun;
        end
    end
    % Define augmented cell description matrices
    pwasys.Ebar{i} = zeros(nE,n+1);
    pwasys.Ebar{i}(:,XIndex) = CellInfo.Ebar{i}(:,1:nlfun.Dim);
    pwasys.Ebar{i}(:,end) = CellInfo.Ebar{i}(:,end);
    pwasys.Ebar{i} = [pwasys.Ebar{i};zeros(1,n) 1];

    % Define augmented boundary description matrices
    for j = 1:NR,
        if isempty(CellInfo.Fbar{i,j}),
            pwasys.Fbar{i,j} = [];
        else
            pwasys.Fbar{i,j} = zeros(n+1,n);
            pwasys.Fbar{i,j}(end,end) = 1;
            pwasys.Fbar{i,j}(NotXIndex,1:length(NotXIndex))=eye(length(NotXIndex));
            pwasys.Fbar{i,j}([XIndex end],length(NotXIndex)+(1:size(CellInfo.Fbar{i,j},2))) = CellInfo.Fbar{i,j};
        end
    end

    % Create B matrix
    if ismatrix,
        % If B(x) is already a matrix, use "as is"
        pwasys.Bbar{i} = [nlsys.B;zeros(1,size(nlsys.B,2))];
    else
        % If B(x) is not a matrix, find the Chebychedf center of the
        % polytopic region described by P={x| Ebar*xbar>0}
        % [Page 148 Convex optimization book by Boyd]

        yalmip('clear');
        r = sdpvar(1);
        xc = sdpvar(n,1);
        xcbar = [xc;1];
        constraints = set([]);
        for j = 1:nT,
            constraints = constraints + set(pwasys.Ebar{i}*xcbar > r*norm(pwasys.Ebar{i}(j,:)));
        end
        objective = -r;
        solvesdp(constraints,[],objective,sdpsettings('verbose',0,'warning',0));
        xcheb = double(xc);
        NaN = find(isnan(xcheb));
        xcheb(NaN) = 0*NaN;
        if isfield(nlsys,'Parameters')
            B{i} = feval(nlsys.B,xcheb,nlsys.Parameters);
        else
            B{i} = feval(nlsys.B,xcheb);
        end
        pwasys.Bbar{i} = [B{i}; zeros(1,size(B{i},2))];
    end
end

% Pass local values to pwasys
pwasys.Abar = pwasys.Abar;
pwasys.Bbar = pwasys.Bbar';
pwasys.Ebar = pwasys.Ebar';
pwasys.Fbar = pwasys.Fbar';
pwasys.Index = XIndex;
pwasys.Rw = pwafun.Rw;
pwasys.W = pwafun.W;
pwasys.Z = pwafun.Z;
pwasys.Err = pwafun.Err;
pwasys.Ind = pwafun.Ind;
pwasys.NR = pwafun.NR;
pwasys.GridPoints = pwafun.GridPoints;
pwasys.xcl = nlsys.xcl;
pwasys.GR = nlsys.GR;
if isfield(nlsys,'Domain'),
    pwasys.Domain = Domain;
end
if isfield(nlsys,'NonlinearDomain'),
    pwasys.NonlinearDomain = nlsys.NonlinearDomain;
end
pwasys.istar = pwa_recreg(pwasys.xcl,pwasys);
pwadi = pwasys;

% Extract A and a from Abar
for i=1:size(pwadi.Abar,1),
    for j=1:2,
        pwadi.A{i,j} = pwadi.Abar{i,j}(1:n,1:n);
        pwadi.a{i,j} = pwadi.Abar{i,j}(1:n,end);
    end
        pwadi.B{i} = pwadi.Bbar{i}(1:n,:);
end

% Compute Cqi and dqi
%
% Input: pwadi.Index         Variables in the domain of nonlinearity
%        pwadi.GridPoints    Uniform grid, breaking points
%        pwadi.Ind           Index of regions
%        pwadi.NR            Number of regions
Index = pwadi.Index;
GridPoints = pwadi.GridPoints;
Ind = pwadi.Ind;
NR = pwadi.NR;

for i=1:NR
    for j=1:length(Index)
        % GridPoints{j}(Ind(j,i))<= x(Index(j)) <= GridPoints{j}(Ind(j,i)+1)
        Cq{i,j} = 2*[zeros(1,Index(j)-1) 1 zeros(1,n-Index(j))]/(GridPoints{j}(Ind(j,i)+1)-GridPoints{j}(Ind(j,i)));
        dq{i,j} = -(GridPoints{j}(Ind(j,i)+1)+GridPoints{j}(Ind(j,i)))/(GridPoints{j}(Ind(j,i)+1)-GridPoints{j}(Ind(j,i)));
    end
end

pwadi.Cq = Cq;
pwadi.dq = dq;