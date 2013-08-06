
l=.2;
mp=1/3;
g=9.8;
mc=1;

ANL = [0 0 1 0;
       0 0 0 1;
       0 0 0 0;
       0 0 0 0];
   
model.A = ANL;
model.aff = [0;0;0;0];
model.Bx = @pendul_gain_nonlinearity;
model.fx = @pendul_dynamics_nonlinearity;
model.Domain = {[-100 100],[3*pi/4 5*pi/4], [-100 100], [-2 2]};
model.NonlinearDomain = [0;1;0;1];
%nlsys.Parameters = nonlinear_function;
model.mtd = 'Uniform';          % Uniform grid
model.NR =3;                     % Unifrom grid resolution
model.xcl=[0;pi;0;0];

%                        END OF MODEL DEFINITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        END OF MODEL DEFINITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------------------------------------------------------------
% Save the model and run it by pressing (F5) to obtian the PWA pproximation
% ------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %%%
%             THIS PART SHOULD NOT BE MODIFIED BY THE USER             %%%
%                                                                      %%%
rehash;
[Err, model]=NonCheckModel(model);
if Err==0
    [pwainc, pwasys, nlsys]=PWAComp(model);
    fprintf('\nNonlinear system saved in nlsys.\n\n');
    fprintf('PWA Differential Inclusion approximation saved in pwainc.\n\n');
    fprintf('PWA approximation saved in pwasys.\n\n');
elseif Err==1
    fprintf('\nI quit because model has error.\n\n');
end
