% PWA Approximation

clear nlsys pwasys DesignParam

model_parameters;

ANL = [0 0 1 0;
    0 0 0 1;
    -M\(Ko+Ku) -M\(Co+Cu)];

BNL = [0 0; 0 0;muu*inv(M)*B];

xcl = [0;0;0;0];

nlsys.A = ANL;
nlsys.B = BNL;
nlsys.a = [0;0;0;0];
nlsys.NonlinearStateEquations = [0;0;1;1];
nlsys.NonlinearFunction = @NLFun;
nlsys.Domain = {[-1 1],[-1 1],[-1 1],[-5 5]};
nlsys.NonlinearDomain = [0;1;0;0];
%nlsys.Parameters = qx2;
 nlsys.Resolution =7;
nlsys.Method = 'Uniform';          % Uniform grid
nlsys.UGR =[4];                     % Unifrom grid resolution
% nlsys.Method = 'OptimalUniform';    % Optimal approximation using a uniform grid
% nlsys.UGR = 5;                      % Unifrom grid resolution
% nlsys.Method = 'Multiresolution'; % Multiresolution grid
% nlsys.TNR =[5 3];                    % Target number of regions
% nlsys.Rstar = [ -.1; .1];
nlsys.ObjFun = 'L2';                % ObjFun: L2 or Linf
nlsys.x0 = xcl;
nlsys.AbarLin = [ANL zeros(4,1);zeros(1,5)];   % Linear approximation of the nonlinear system
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




