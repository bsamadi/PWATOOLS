
function pwasys=MinElp(pwasys)
% This function computes the minimum volume ellipsoid for slab and polytopic
% regions. 
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

E=pwasys.E;
n=size(E{1},2);
NR=pwasys.NR;
RX=pntrgn(pwasys);  % points at the boundaries 
dE1=size(E{1}, 1);  % number of equation in E

if dE1==2  % slab regions
    for i=1:NR
        X=RX{i};
        direct=any(X'); % direction at which region is defiend; equivalent to NonlinearDomain
        Y=direct*X;
        pwasys.EpA{i}=2*direct/(Y(2)-Y(1));
        pwasys.Epb{i}=-(Y(1)+Y(2))/(Y(2)-Y(1));
    end
elseif dE1>2   % polytopic
    for i=1:NR
        X=RX{i};
        [X, Izero]=compactx(X);
        n=size(X,1);
        X_dim_2=size(X, 2);
        EpA=sdpvar(n,n);
        Epb=sdpvar(n,1);
        constraints=[];
        for k=1:X_dim_2
            constraints=constraints + set([eye(n) EpA*X(:,k)+Epb;(EpA*X(:,k)+Epb)'  1] >0);
        end
        solvesdp(constraints, -geomean2(EpA))
        EpA=double(EpA);
        Epb=double(Epb);
        [EpA, Epb]=expandx(EpA, Epb, Izero);
        pwasys.EpA{i}=EpA;
        pwasys.Epb{i}=Epb;
    end
    
end
end
%% local function 1: Compact X
function [X, Izero]=compactx(X)
% this function remobes zero from the rows which are zero and decreases the
% state space dimension

SysemOrder=size(X, 1);

Izero=[];
for i=1:SysemOrder
    RowZero=all(~X(i, :));
    if RowZero
     Izero=[Izero i];
    end
end
X(Izero, :)=[];

end

%% local function 2: Expand X
function [EpA, Epb]=expandx(A, b, Izero)
% this function returns matrices A and b to the original state space

SystemOrder=size(A, 1)+length(Izero);
Icomplement=setdiff([1:SystemOrder], Izero);

EpA=zeros(SystemOrder);
Epb=zeros(SystemOrder, 1);


EpA(Icomplement, Icomplement)=A;
Epb(Icomplement)=b;

end
