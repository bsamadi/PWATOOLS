function M = pwa_mag(A)

% Customized magnitude of a matrix
% M = pwa_mag(A)

A = A'.*A';
if size(A,1)>1,
    A = sum(A);
end
M = sqrt(A);