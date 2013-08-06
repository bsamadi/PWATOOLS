function [Y,Err,Obj,T,Abar] = pwai_estim(W,Z,X,Option)

% PWA Estimation
% [Y Err Obj T Abar] = pwa_estim(W,Z,X,Option,Abarstar)

[m n] = size(X);
p = size(W,1);
q = size(Z,2);

Ver = ver('MATLAB');
if str2num(Ver.Version(1))==7,
    T = delaunayn(X,{'Qt','Qbb','Qc','Qz'});
else
    T = delaunayn(X);
end

numt = size(T,1);

constraints=set([]);

for i=1:numt,
    Abar{i} = sdpvar(q,n+1,'full'); % PWA approximation parameters
end

% Linear approximation
if isfield(Option,'xcl') && isfield(Option,'AbarLin'),
    if n == 2,
        t = tsearch(X(:,1),X(:,2),T,Option.xcl(1),Option.xcl(2));
    else
        t = tsearchn(X,T,Option.xcl');
    end
    if any(ismember(X,Option.xcl','rows')),
        i = find(ismember(X,Option.xcl','rows'));
        [t,j]=find(T==i);
    end
    for i=1:length(t),
        constraints = constraints + set(Abar{t(i)} == Option.AbarLin);
    end
end

% Continuity constraints
if numt>1,
    for i = 1:m,
        Ti=find(any((T==i)'));
        if length(Ti)>1,
            for j=1:length(Ti)-1,
                k = j+1;
                constraints=constraints+set((Abar{Ti(k)}-Abar{Ti(j)})*[X(i,:) 1]'==0,['Abar' num2str(Ti(k)) '*x=Abar' num2str(Ti(j)) '*x']);
            end
        end
    end
end

% Error
if n == 2,
    t = tsearch(X(:,1),X(:,2),T,W(:,1),W(:,2));
else
    t = tsearchn(X,T,W);
end

Err = [];
for i = 1:p,
    Err = [Err;Z(i,:)-(Abar{t(i)}*[W(i,:) 1]')'];
end

% Objective
if ~isfield(Option,'ObjFun'),
    Option.ObjFun='L2';
end
if strcmp(Option.ObjFun,'L2'),
%     Obj=sum(sum(Err.*Err));
        Obj = sdpvar(1,1);
        ErrVec = reshape(Err,numel(Err),1);
        Schur = [Obj ErrVec';ErrVec eye(length(ErrVec))];
        constraints=constraints+set(Schur>0,'Shur>0');
elseif strcmp(Option.ObjFun,'Linf'),
    Obj = sdpvar(1,1);
    constraints=constraints+set(Obj>0,'Obj>0');
    constraints=constraints+set(Err>-Obj,'Error>-Obj');
    constraints=constraints+set(Err<Obj,'Error<Obj');
end

% Solve the optimization problem
solution=solvesdp(constraints,Obj);

for i=1:length(Abar),
    Abar{i} = double(Abar{i});
end
Obj = double(Obj);
Err = double(Err);

Y = [];
for i=1:m,
    J=find(any((T==i)'));
    Y=[Y; (Abar{J(1)}*[X(i,:) 1]')'];
end