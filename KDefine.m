function [pwasys, sysloop]=KDefine(pwasys, Gain)
% This function ensures the controller gains K and k are either set or
% zero with proper sizes!
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%

A=pwasys.A;
a=pwasys.a;
B=pwasys.B;
n=size(pwasys.Abar{1,1},1)-1;
m=size(pwasys.Bbar{1,1},2);
NR=size(A,1);

if nargin==2
    pwasys.Kbar=Gain;  
end

if isfield(pwasys, 'Kbar')
    for i=1:NR
        pwasys.K{i}=pwasys.Kbar{i}(:, 1:end-1);
        pwasys.k{i}=pwasys.Kbar{i}(:,end);
    end
else
    for i=1:NR
        try  % ensuring K is defined properly
            if isempty(pwasys.K{i})
                pwasys.K{i}=zeros(m,n);
            end
        catch
            pwasys.K{i}=zeros(m,n);
        end
        %-------------------------------------
        try  % ensuring k is defined properly
            if isempty(pwasys.k{i})
                pwasys.k{i}=zeros(m,1);
            end
        catch
            pwasys.k{i}=zeros(m,1);
        end
        pwasys.Kbar{i}=[pwasys.K{i}  pwasys.k{i}];
    end
end

Norm=0;

for i=1:NR
    Norm=Norm+norm(pwasys.K{i})+norm(pwasys.k{i});
end

if Norm==0
    sysloop='open';
else
    sysloop='closed';
end


