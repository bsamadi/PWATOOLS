function pwasys = center_region(pwasys)

% Returns pwasys.istar that contains the index of all regions that include
% the desired equilibrium point: pwasys.xcl

[pwasys.NR pwasys.NS] = size(pwasys.Abar);    % number of regions, number of systems per region
pwasys.n = size(pwasys.Abar{1},1)-1;          % number of state variables
pwasys.m = size(pwasys.Bbar{1},2);            % number of inputs
pwasys.p = size(pwasys.Ebar{1},1)-1;          % number of constraints for each region

pwasys.istar = [];

for i=1:pwasys.NR,
    xcl_is_inside_Ri = all(pwasys.Ebar{i}*[pwasys.xcl; 1]>=0-1e-7);
    if xcl_is_inside_Ri,
        pwasys.istar = union(pwasys.istar,i);
    end
end