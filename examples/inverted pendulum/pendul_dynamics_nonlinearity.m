function F=pendul_dynamics_nonlinearity(X)


l=.2;
mp=1/3;
g=9.8;
mc=1;


x2=X(2);
x4=X(4);


M = mc+mp*(1-3/4*[cos(x2)]^2);
f1= 1/M*(mp*l/2*x4^2*sin(x2)+3*mp*g/8*sin(2*x2));
f2= -3/(2*l)*(g*sin(x2)+ f1*cos(x2));
F=[0;0;f1;f2];



