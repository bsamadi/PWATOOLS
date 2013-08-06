function PWP = pwp_approximate(NLFcn)

% Usage: PWP = pwp_approximate(NLFcn)
% Example:
% NLFcn.Handle = @nlfcn;
% NLFcn.Resolution = Resolution;
% NLFcn.Domain = Domain;
% NLFcn.Order = nl;
% 
% PWP = pwp_approximate(NLFcn);
% coeff = PWP.coeff;
% X = PWP.X;
% NLFun = PWP.NLFun;
% PWPFun = PWP.PWPFun;
% 
% figure(100);
% line(X,NLFun,'Color',[0 0 0])
% line(X,PWPFun,'Color',0.7*[1 1 1],'LineWidth',2,'LineStyle','--');
% xlabel('$x$');
% ylabel('$f(x)$');
% legend('$f(x)$','$\hat f(x)$',2);
% set(gca,'Box','on','XTick',Domain,'XLim',[Xmin Xmax],'XGrid','on');


Domain = NLFcn.Domain;
Resolution = NLFcn.Resolution;
nl = NLFcn.Order;
xstar = NLFcn.xstar;

NR = length(Domain)-1;
X = sort(union(linspace(Domain(1),Domain(end),Resolution),Domain));

F = [];
x = sdpvar;
[p,c,v] = polynomial(x,nl);
P = sdpvar(1,length(v));
Funx = char(sdisplay(P*v));
FunX = strrep(Funx,'x','X(In)');
FunX = strrep(FunX,'^','.^');

% Approximation Error
for i = 1:NR,
    In = find(X>=Domain(i)&X<=Domain(i+1));
    NLFun(In) = feval(NLFcn.Handle,X(In));
    coeff{i} = sdpvar(1,length(v));
    P = coeff{i};
    PWPFun(In) = eval(FunX);
end

Err = PWPFun-NLFun;
Obj = sdpvar(1);
F = [F, [Obj Err;Err' eye(length(Err))]>0];

% Continuity
for i=1:NR-1,
    x = Domain(i+1);
    P = coeff{i};
    PWPFuni = eval(Funx);
    P = coeff{i+1};
    PWPFunip1 = eval(Funx);
    F = [F, PWPFuni==PWPFunip1];
end

Istar = min(find(Domain>=xstar))-1;
x = xstar;
P = coeff{Istar};
PWPFuni = eval(Funx);
NLFunc = feval(NLFcn.Handle,x);

F = [F, PWPFuni==NLFunc];

sol = solvesdp(F,Obj);

for i = 1:NR,
    In = find(X>=Domain(i)&X<=Domain(i+1));
    NLFun(In) = feval(NLFcn.Handle,X(In));
    coeff{i} = double(coeff{i});
    P = coeff{i};
    PWPFun(In) = eval(FunX);
end

PWP.coeff = coeff;
PWP.X = X;
PWP.NLFun = NLFun;
PWP.PWPFun = PWPFun;
PWP.NR = NR;