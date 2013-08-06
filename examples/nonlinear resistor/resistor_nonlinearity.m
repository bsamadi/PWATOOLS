function F = resistor_nonlinearity(X)

x=X(2);

if (-2e4 <= x & x< .2)
    F=[0; -50*(1e-3/.2*x)];
elseif (.2 <= x & x<.6)
    F=[0; -50*(-2e-3*x+1.4e-3)];
elseif (.6 <= x & x <= 2e4) 
    F=[0; -50*(4e-3*x-2.2e-3)];
end


