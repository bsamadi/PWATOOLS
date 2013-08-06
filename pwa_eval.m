function S = pwa_eval(X,Y,W)

% PWA function evaluation
% S = pwa_eval(X,Y,W)

[m n] = size(X);
[mW nW] = size(W);

if n==mW & n~=nW,
    W = W';
end

if n > 1,
    T = DelaunayTri(X);
    t = pointLocation(T,W);
    if isnan(t),
        error('The given point is out of the function domain.');
    end
    for i = 1:length(t),
        A = inv([X(T(t(i),:),:) ones(n+1,1)])*Y(T(t(i),:),:);
        S(i,:) = [W(i,:) 1]*A;
    end
else
    [X i] = sort(X);
    Y = Y(i,:);
    S = [];
    for i = 1:size(Y,2),
        S = [S interp1(X,Y(:,i),W,'linear')];
    end
end