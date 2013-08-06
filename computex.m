function X=computex(pwasys)
% this function computs X from given region equations E and e
% each column contains a point which is the intersection of two boundaries 
% Copyright: Mohsen Zamani Fekri, Concordia University April 2011
%
%


E = pwasys.E;
e = pwasys.e;
NR=pwasys.NR;
dE1=size(E{1}, 1);  % number of equation in E



if dE1==2    % slab regions
    X=[];
    for i=1:NR
        for j=1:2
            pnt=-E{i}(j,:)\e{i}(j);
            X=[X pnt];
        end
    end
    X=round(X*10^10)./10^10;
    Y=unique(X', 'rows');
    X=Y';
elseif dE1>2 % polytopic regions
    X=[];
    index=combntns(1:dE1, 2);  %all possible systems of equatons resulted from E
    L=size(index, 1);
    for i=1:NR
        for j=1:L
            EA=E{i}(index(j,:), :);
            Ee=e{i}(index(j,:));
            pnt=-EA\Ee;
            Er=EA*pnt+Ee;
            if norm(Er)<1e-6
                X=[X pnt];
            end
        end
        X=round(X*10^10)./10^10;
        Y=unique(X', 'rows');
        X=Y';
    end
else         % unbounded regions
    fprintf('I cannot handle unbounded regions.\n')
    return;
end


