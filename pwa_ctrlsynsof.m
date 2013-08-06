function pwactrl = pwa_ctrlsynsof(pwasys, pwactrl,Param)

% State Feedback Controller Synthesis
% pwactrl = pwa_ctrlsynsf(pwasys,Param)

% pwasys
% xbarDot = Abar{i}*xbar + Bbar{i}*u
%
% pwasys.Abar = {Abar1, ..., Abarn}
% pwasys.Bbar = {Bbar1, ..., Bbarn}
% pwasys.Cybar = {Cybar1, ..., Cybarn}
% pwasys.Czbar = {Czbar1, ..., Czbarn}
% pwasys.Dzw = {Dzw1, ..., Dzwn}
% pwasys.Dzu = {Dzu1, ..., Dzun}
%
% Polytopic region i: {x| Ebar{i}*xbar > 0 }
%
% Boundary information table:
% Boundary between Ri and Rj : {x| xbar=Fbar{i,j}*sbar }
%
% Param.Lyapunov = 'SOS'
% Param.Order = Order of SOS function (should be even)
% Param.LimitC = Limit for polynomial coefficients

pwasys = center_region(pwasys);
x = sdpvar(pwasys.n,1);
xbar = [x;1];

% Polynomial Lyapunov
[V,c] = polynomial(x,Param.Order);
F = set(-Param.LimitC<c(:)<Param.LimitC);

% Controller
for i = 1:pwasys.n,
    pwactrl.Kbar{i} = sdpvar(1,pwasys.n+1);
    F = F + set(-Param.LimitK<pwactrl.Kbar{i}<Param.LimitK);
end

% Positive definite
eps1 = sdpvar(Param.Order/2,1);
F = F + set(eps1(:)>0);
LB = 0;
for i=1:Param.Order/2,
    LB = LB+eps1(i)*(x'*x)^i;
end
F = F + set(sos(V-LB));    % Lower Bound


for i = pwasys.istar,
    A{i} = pwasys.Abar{i}(1:pwasys.n,1:pwasys.n);
    B{i} = pwasys.Bbar{i}(1:pwasys.n,:);
end
pwactrl.Lin2PWA = 1;
if isfield(Param,'Poles'),
    for i = pwasys.istar,
        K{i} = -acker(A{i},B{i},Param.Poles);
    end
elseif isfield(Param,'KLin'),
    for i = pwasys.istar,
        K{i} = Param.KLin;
    end
elseif isfield(Param,'QLin'),
    for i = pwasys.istar,
        K{i} = -lqr(A{i},B{i},Param.QLin,Param.RLin);
    end
else
    pwactrl.Lin2PWA = 0;
end

if pwactrl.Lin2PWA,
    for i=pwasys.istar,
        k = sdpvar(pwasys.m,1);
        pwactrl.Kbar{i} = [K{i} k];
    end
end

for i=pwasys.istar,
    F = F+set( (pwasys.Abar{i}+pwasys.Bbar{i}*pwactrl.Kbar{i})*[pwasys.xcl; 1] == 0,['(Abar+Bbar*Kbar)*xclbar=0']);
end

for i=1:pwasys.NR,
    for j = i:pwasys.NR,
        if ~isempty(pwasys.Fbar{i,j}),
            % Continuity of the control input
            F = F + set((pwactrl.Kbar{i}-pwactrl.Kbar{j})*pwasys.Fbar{i,j}== 0,['Ki*F' num2str(i) '_' num2str(j) 's=Kis*F' num2str(i) '_' num2str(j) 's']);
        end
    end
end

for i = 1:pwasys.NR,
    % Description of the regions
    S = pwasys.Ebar{i}*xbar;
    S = S(1:end-1);

    % S procedure parameters
    for j = 1:length(S),
        [L(j),cs{j}] = polynomial(x,Param.Order-2);
        F = F + set(-Param.LimitC<cs{j}(:)<Param.LimitC);
        F = F + set(sos(L(j)));
    end

    % Stability constraint dV/dt < -alphaV
    xDot = (pwasys.Abar{i}+pwasys.Bbar{i}*pwactrl.Kbar{i})*xbar;
    F = F + set(sos(-jacobian(V,x)*xDot(1:end-1)-L*S-Param.alpha*V));
end

[sol,m,B] = solvesos(F,[]);
checkset(F)

for i=1:pwasys.NR,
    pwactrl.Kbar{i}=double(pwactrl.Kbar{i});
end