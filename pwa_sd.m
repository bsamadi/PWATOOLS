function SD = pwa_sd(pwasys,pwactrl)

tau = sdpvar(1);
gamma = sdpvar(1);
eta = sdpvar(1);

M = length(pwasys.Abar);
n = length(pwasys.Abar{1})-1;
m = size(pwasys.Bbar{1},2);

P = sdpvar(n,n);
X = sdpvar(n,n);
R = sdpvar(n,n);
I = eye(n);
Im = eye(m);

F = set(P>0)+set(R>0)+set(X>0)+set(eta>gamma)+set(gamma>10);
F = F+set(P(:)<1e4)+set(R(:)<1e4)+set(X(:)<1e4);

for i=1:M,
    if ismember(i,pwasys.istar),
        A = pwasys.Abar{i}(1:n,1:n);
        B = pwasys.Bbar{i}(1:n,:);
        K = pwactrl.Kbar{i}(:,1:n);
        Fbar{i} = [A B*K];
        N{i} = sdpvar(2*n,2);

        Psi{i} = [P;zeros(n,n)]*Fbar{i}+Fbar{i}'*[P zeros(n,n)]-[I;-I]*X*[I -I]-2*N{i}*[I -I]-2*[I;-I]*N{i}'+eta*[I zeros(n,n);zeros(n,n) I];
        M1{i} = [I;-I]*X*Fbar{i}+Fbar{i}'*X*[I -I];
        M2{i} = Fbar{i}'*R*Fbar{i};

        F = F+set([Psi{i}+tau*M1{i} [P;zeros(n,n)]*B+tau*[I;-I]*X*B; B'*[P zeros(n,n)]+tau*B'*X*[I -I] -gamma*Im]<0);
        F = F+set([Psi{i}+tau*M2{i} tau*Fbar{i}'*R*B+[P;zeros(n,n)]*B tau*N{i}; tau*B'*R*Fbar{i}+B'*[P zeros(n,n)] tau*B'*R*B-gamma*Im zeros(m,n); tau*N{i}' zeros(n,m) -tau*R/2]<0);

%       Linear systems
%         F = F+set(Psi{i}+tau*M1{i}<0);
%         F = F+set([Psi{i}+tau*M2{i} tau*N{i};tau*N{i}' -tau*R/2]<0);

    else,
        A = pwasys.Abar{i}(1:n,1:n);
        a = pwasys.Abar{i}(1:n,n+1);
        B = pwasys.Bbar{i}(1:n,:);
        K = pwactrl.Kbar{i}(:,1:n);
        k = pwactrl.Kbar{i}(:,n+1);
        L = pwasys.L{i};
        l = pwasys.l{i};
        
        Fbar{i} = [A B*K B*k+a];
        N{i} = sdpvar(2*n,2);
        lambda{i} = sdpvar(1);


        Psi{i} = [P;zeros(n+1,n)]*Fbar{i}+Fbar{i}'*[P zeros(n,n+1)]-[I;-I;zeros(1,n)]*X*[I -I zeros(n,1)]-2*[N{i};zeros(1,n)]*[I -I zeros(n,1)]-2*[I;-I;zeros(1,n)]*[N{i}' zeros(n,1)]+[eta*[I zeros(n,n);zeros(n,n) I] zeros(2*n,1);zeros(1,2*n) -lambda{i}]+lambda{i}*[L';zeros(n,1);l]*[L zeros(1,n) l];
        M1{i} = [I;-I;zeros(1,n)]*X*Fbar{i}+Fbar{i}'*X*[I -I zeros(n,1)];
        M2{i} = Fbar{i}'*R*Fbar{i};
        
        F = F+set(lambda{i}<0);
        F = F+set([Psi{i}+tau*M1{i} [P;zeros(n+1,n)]*B+tau*[I;-I;zeros(1,n)]*X*B; B'*[P zeros(n,n+1)]+tau*B'*X*[I -I zeros(n,1)] -gamma*Im]<0);
        F = F+set([Psi{i}+tau*M2{i} tau*Fbar{i}'*R*B+[P;zeros(n+1,n)]*B tau*[N{i}; zeros(1,n)]; tau*B'*R*Fbar{i}+B'*[P zeros(n,n+1)] tau*B'*R*B-gamma*Im zeros(m,n); tau*[N{i}' zeros(n,1)] zeros(n,m) -tau*R/2]<0);
    end
end

sol = solvesdp(F,-tau);
checkset(F);
sol

SD.TauMax = double(tau);
SD.gamma = double(gamma);
SD.eta = double(eta);
SD.P = double(P);
SD.R = double(R);
SD.X = double(X);
