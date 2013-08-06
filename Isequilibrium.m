function Isxcl=Isequilibrium(pwasys, xcl, tol)
% this function determines a given point xcl is an open/closed-loop
% equilibrium point. If so, Isxcl accordingly will be one.
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

Isxcl.open=0;
Isxcl.closed=0;

NR=pwasys.NR;
A=pwasys.A;
B=pwasys.B;
a=pwasys.a;
K=pwasys.K;
k=pwasys.k;

if strcmp(pwasys.type, 'pwadi')
    col_index=[1 2];
elseif strcmp(pwasys.type, 'lower-envelope')
    col_index=[1];
end

I=[];  % indices of the region which contain xcl
for i=1:NR
    reg_Er=pwasys.E{i}*xcl+pwasys.e{i};
    reg_Rr=round(reg_Er*10^8)/10^8;
    if min(sign(reg_Er))>= 0
        I=[I i];
    end
end


Y=[];
Z=[];
for i=I
    for j=col_index
        % open loop
        y=A{i,j}*xcl+a{i,j};
        Y=[Y y];
        % closed loop
        z=(A{i,j}+B{i,j}*K{i})*xcl+a{i,j}+B{i,j}*k{i};
        Z=[Z z];
    end
end
if norm(Y)<tol
    Isxcl.open=1;   % one if xcl is an open-loop equilibrium point
end

if norm(Z)<tol
    Isxcl.closed=1; % one if xcl is a closed loop equilibrium point
end




