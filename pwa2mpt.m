function sysStruct = pwa2mpt(pwasys,Ts)

N = length(pwasys.Abar);
n = length(pwasys.Abar{1})-1;
m = size(pwasys.Bbar{1},2);

for i=1:N,
    sysStruct.A{i} = eye(n)+Ts*pwasys.Abar{i}(1:n,1:n);
    sysStruct.B{i} = Ts*pwasys.Bbar{i}(1:n,:);
    sysStruct.f{i} = Ts*pwasys.Abar{i}(1:n,n+1);
    sysStruct.C{i} = eye(n);
    sysStruct.D{i} = zeros(n,m);
    sysStruct.guardX{i} = -pwasys.Ebar{i}(1:n,1:n);
    sysStruct.guardC{i} = pwasys.Ebar{i}(1:n,n+1);
end