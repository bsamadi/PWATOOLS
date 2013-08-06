% PWA Approximation
% This example was not tried. 

clear nlsys pwasys DesignParam


nonlinear_function='F=[0; -.5*x1^2*x2];';
ANL = [0 1;-1 .5];
BNL = [0;1];
xcl = [0;0];
nlsys.A = ANL;
nlsys.B = BNL;
nlsys.a = [0;0];
nlsys.NonlinearStateEquations = [0 ; 1];
nlsys.NonlinearFunction = @NLFun_Extension_Example_6_3;
nlsys.Domain = {[-30 30],[-60 60]};
nlsys.NonlinearDomain = [1;1];
nlsys.Parameters = nonlinear_function;
nlsys.Resolution = 9;
nlsys.Method = 'Uniform';          % Uniform grid
nlsys.UGR =3;                     % Unifrom grid resolution
%nlsys.Method = 'OptimalUniform';    % Optimal approximation using a uniform grid
%nlsys.UGR = 3;                      % Unifrom grid resolution
%nlsys.Method = 'Multiresolution'; % Multiresolution grid
%nlsys.TNR =3;                    % Target number of regions
nlsys.Rstar = [-2; 0];
nlsys.ObjFun = 'L2';                % ObjFun: L2 or Linf
nlsys.x0 = [0; .1];
nlsys.AbarLin = [ANL nlsys.a;zeros(1,3)];   % Linear approximation of the nonlinear system
nlsys.xcl = xcl;

% PWA approximation
pwasys = pwa_approximate(nlsys);
% PWA Envelope
nlsys.X = pwasys.X;
nlsys.Method = 'UpperBound';
pwasysU = pwa_approximate(nlsys);

nlsys.Method = 'LowerBound';
pwasysL = pwa_approximate(nlsys);

pwasys = pwasysU;

pwasys.Abar = [pwasysL.Abar  pwasysU.Abar];
pwasys.Bbar = [pwasysL.Bbar  pwasysU.Bbar];

pwasys.X = pwasysL.X;
pwasys.W = pwasysL.W;

pwasys.Z = {pwasysL.Z, pwasysU.Z};
pwasys.Y = {pwasysL.Y, pwasysU.Y};

