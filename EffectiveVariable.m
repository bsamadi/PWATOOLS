function [eff, noneff]= EffectiveVariable(pwasys)
% this function determines which variables are used in regions equations.
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

NR=size(pwasys.A, 1);
n=size(pwasys.A{1}, 1);


S=[];
for i=1:NR
    S=[S; pwasys.E{i}];
end

eff=find(any(S));
noneff=setdiff([1:n], eff);


