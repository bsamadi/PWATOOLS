function B=pendul_gain_nonlinearity(X)


l=.2;
mp=1/3;
g=9.8;
mc=1;

x2=X(2);

M = mc+mp*(1-3/4*[cos(x2)]^2);
b3=1/M;
b4=-3/(2*l*M)*cos(x2);

B=[0;0; b3; b4];

