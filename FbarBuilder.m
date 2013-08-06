function F=FbarBuilder(pwasys)
% This function synthesizes Fbar from Ebar
% it computes the normal direction and calculates the first effective
% variable in terms of other variables.
% effective and noneffective variables are defined below
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

NR=pwasys.NR;

RX=pntrgn(pwasys); % points in each region

% DI: effective variables
% DIC: noneffective variables
[DI, DIC]= EffectiveVariable(pwasys);

n=size(pwasys.A{1,1},1);

S=[];  % for noneffective variables
for h=DIC
    temp=zeros(n,1);
    temp(h)=1;
    S=[S temp];
end

for i=1:NR % compacting state space
    RX{i}(DIC, :) =[];
end

F{1,1}=[];
for i = 1:NR,
    for j = i+1:NR
            F{i,j}=[];
            F{j,i}=[];
        pointscommon = intersect(RX{i}', RX{j}', 'rows');
        L=size(pointscommon, 1);
        if ~isempty(pointscommon) && L==length(DI),
            normaldirection=null([pointscommon ones(L,1)]);
            normaldirection=-normaldirection/normaldirection(1); % finding the first effective variable in terms ogf other
            affine=zeros(n,1);
            affine(DI(1))=normaldirection(end);   % 
            G=[];
            index=1;
            for y=DI(2:end)
                index=index+1;
                temp=zeros(n,1);
                temp(DI(1))=normaldirection(index);  % 
                temp(y)=1;
                G=[G temp];
            end
            F{i,j}=[S G affine; zeros(1, n-1) 1];
            F{j,i}=F{i,j};
        end
    end
end


