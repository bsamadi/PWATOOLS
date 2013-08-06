function pwafun = pwa_split(nlfun)

% PWA_Split Piecewise-affine approximation.
%
% pwafun = pwa_split(nlfun) computes a piecewise-affine (PWA) approximation
% of a nonlinear function,
%
% nlfun.Handle = @function_name;
% nlfun.Domain = {[xmin xmax]};  % It's a cell vector
% nlfun.TNR = Target number of regions;
% nlfun.RStar = The center region; % Optional
% nlfun.Parameters = Parameter vector; % If required
% nlfun.Resolution = Fine grid resolution;
% nlfun.ObjFun = 'L2'; % or 'Linf'
% nlfun.AP = Anchoring point; % Optional
% nlfun.AbarLin = Linear approximation of the function at the anchoring point;  % Optional
% nlfun.Index = Index of the variables in the domain of nonlinearity
% nlfun.YIndex = Index of the nonlinear entries
%
% pwafun.X = Griding points
% pwafun.Y = The value of the pwafun at griding points
% pwafun.T = Triangles
% pwafun.Abar = Linear coefficient for each triangle; % Y = Abar{i}*X
% pwafun.Err = Approximation error for the fine mesh
% pwafun.W = Fine mesh points
% pwafun.Z = The value of pwafun at W
% pwafun.Obj = The value of the objective function
%
% Copyright: Behzad Samadi, Concordia University June 2005, April 2007


% Argument checking

f = nlfun.Handle;                     % Specify the handle of the nonlinear function

a=[];
b=[];
for i=1:length(nlfun.Domain),
    a=[a; nlfun.Domain{i}(1)];    % Set start point of the domain
    b=[b; nlfun.Domain{i}(end)];  % Set start point of the domain
end

if isfield(nlfun,'Parameters'),
    Param = nlfun.Parameters;             % Set the parameters
else
    Param = [];
end

TNR = nlfun.TNR;                      % The target number of regions
Option.ObjFun = nlfun.ObjFun;         % Specify the objective: L2 or Linf
if isfield(nlfun,'xcl') && isfield(nlfun,'AbarLin'),
    Option.xcl = nlfun.xcl;                % Specify the anchoring point
    Option.AbarLin = nlfun.AbarLin;      % Specify the linear approximation at the anchoring point
end
n=length(a);

try
    Y = pwa_nlfun(nlfun,a,Param);
catch
    error(['The function can not be evaluated.']);
end

% Construct the fine mesh
%
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
Z = pwa_nlfun(nlfun,W,Param);

% Initial Approximation Mesh
X = pwa_grid(Domain,1);
if isfield(nlfun,'Rstar'),
   X = [X;nlfun.Rstar];
end
[Y Err Obj T Abar] = pwa_estim(W,Z,X,Option);

while size(T,1)<TNR,
    Y = pwa_eval(W,Z,X);
    S = pwa_eval(X,Y,W);
    Err = Z-S;
    ErrMag = pwa_mag(Err);
    if isfield(nlfun,'Rstar'),
        
        tri = find(1-isnan(tsearchn(nlfun.Rstar,1:length(nlfun.Rstar),W)));
        ErrMag(tri) = 0*tri;
    end
    [Max i]=max(ErrMag);
    
    NewX = union(X, W(i,:),'rows');
    
    if size(NewX,1)==size(X,1),
        break; % If splitting does not increase the accuracy.
    else,
        X = NewX;
    end

    Ver = ver('MATLAB');
    if str2num(Ver.Version(1))==7,
        T = delaunayn(X,{'Qt','Qbb','Qc','Qz'});
    else
        T = delaunayn(X);
    end

    if n == 1,
        X = sort(X);
    end
    [Y Err Obj T Abar] = pwa_estim(W,Z,X,Option);
end

pwafun.X = X;
pwafun.Y = Y;
pwafun.T = T;
pwafun.Err = Err;
pwafun.Obj = Obj;
pwafun.W = W;
pwafun.Z = Z;
pwafun.Abar = Abar;

% pwa_split