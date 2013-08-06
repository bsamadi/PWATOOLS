function pwainc = pwa_rectangular(nlfun)

% pwa_rectangular Piecewise-affine inclusion for a rectangular grid.
%
% pwainc = pwa_rectangular(nlfun)
%
% nlfun.Handle = @function_name;
% nlfun.Domain = [xmin xmax];
% nlfun.GR = Uniform grid resolution;
% nlfun.Resolution = Fine grid resolution; % (Optional)
% nlfun.ObjFun = 'L2'; % or 'Linf' (Optional)
%
% pwainc.X = Griding points
% pwainc.Y = The value of the pwainc at griding points
% pwainc.T = Triangles
% pwainc.Abar = Linear coefficients; % Y = Abar{i}*X
% pwainc.Err = Approximation error for the fine mesh
% pwainc.W = Fine mesh points
% pwainc.Z = The value of pwainc at W
% pwainc.Obj = The value of the objective function
%
% Copyright: Behzad Samadi, Amirkabir University of Technology, Oct. 2008

% Argument checking
f = nlfun.Handle;                     % Specify the handle of the nonlinear function
n = length(nlfun.Domain);             % Number of variables in the domain of nonlinearities

a=[];
for i=1:length(nlfun.Domain),
    a=[a; nlfun.Domain{i}(1)];    % Set start point of the domain
end

if isfield(nlfun,'Parameters'),
    Param = nlfun.Parameters;             % Set the system parameters
else
    Param = [];
end

if isfield(nlfun,'ObjFun'),
    Option = nlfun.ObjFun;                % Specify the objective (L2 or Linf)
end

%
for i = 1:n,
    if any(diff(nlfun.Domain{i})<0),
        eval(['error(''nlfun.domain{' num2str(i) '} should be increasing.'');'])
    end
end

%
try
    Y = pwa_nlfun(nlfun,a,Param);
catch
    error(['The function can not be evaluated.']);
end

% The Approximation

% Create coarse mesh
X = pwa_grid(nlfun.Domain,nlfun.GR);
if n == 1,
    X = sort(X);
end

% Evaluate nlfun using coarse mesh points
Y = pwa_nlfun(nlfun,X,Param);

Ver = ver('MATLAB');
if str2num(Ver.Version(1))==7,
    T = delaunayn(X,{'Qt','Qbb','Qc','Qz'});
else
    T = delaunayn(X);
end
[mT nT] = size(T);

% Create Abar{i}
for i = 1:mT,
    y = [];
    x = [];
    for j = 1:nT,
        y = [y Y(T(i,j),:)'];
        x = [x [X(T(i,j),:)'; 1]];
    end
    Abar{i} = y/x;
    if isempty(Abar{i}),
        Abar{i} = 0;
    end
end

% Pass local values to pwainc
pwainc.X = X;
pwainc.Y = Y;
pwainc.T = T;
pwainc.Abar = Abar;

if ~isempty(nlfun.Resolution),
    % Construct the fine mesh
    % Z = f(W)
    % W = [w1; w2; ... ; wm] is m times n
    % Z = [z1; z2; ... ; zm] is m times p
    %
    % zi = f(wi)

    % Create the fine mesh
    for i=1:length(nlfun.Domain), 
        Domain{i} = nlfun.Domain{i}([1 end]);
    end
    W = pwa_grid(Domain,nlfun.Resolution);
    % Evaluate nlfun using fine mesh points
    Z = pwa_nlfun(nlfun,W,Param);

    % Define error function
    p = size(W,1);
    Err = [];
    for i = 1:p,
        if ~isempty(Z),
            Err = [Err;Z(i,:)-pwa_eval(X,Y,W(i,:))];
        end
    end
    pwainc.Err = Err;
    pwainc.W = W;
    pwainc.Z = Z;
end

if exist('Option'),
    % Set objective function
    if strcmp(Option,'L2'),
        Obj = sum(sum(Err.*Err));
    elseif strcmp(Option,'Linf'),
        Obj = max(max(abs(Err)));
    end
    pwainc.Obj = Obj;
end