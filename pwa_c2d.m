function pwadsys = pwac2d(pwasys,Ts)

pwadsys = pwasys;
pwadsys.Ts = Ts;

N = length(pwasys.Abar);
n = length(pwasys.Abar{1})-1;
m = size(pwasys.Bbar{1},2);

for i=1:N,
    pwadsys.A{i} = eye(n)+Ts*pwasys.Abar{i}(1:n,1:n);
    pwadsys.B{i} = Ts*pwasys.Bbar{i}(1:n,:);
    pwadsys.f{i} = Ts*pwasys.Abar{i}(1:n,n+1);
    pwadsys.C{i} = eye(n);
    pwadsys.D{i} = zeros(n,m);
    pwadsys.guardX{i} = -pwasys.Ebar{i}(1:n,1:n);
    pwadsys.guardC{i} = pwasys.Ebar{i}(1:n,n+1);
    pwadsys.Abar{i} = [pwadsys.A{i} pwadsys.f{i}; zeros(1,n) 0];
    pwadsys.Bbar{i} = [pwadsys.B{i}; zeros(1,m) ];
end