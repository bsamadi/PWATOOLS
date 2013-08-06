clear nlsys pwasys pwainc

% parameters are assigned
model_parameters;
model.xcl = [0;0;0;0];
model.A = [0     0      1 0;
           0     0      0 1;
           -M\(Ko+Ku) -M\(Co+Cu)];
model.Bx = [0 0; 0 0;muu*inv(M)*B];
model.aff = [0;0;0;0];
model.fx = @Flutter_Nonlinearity;
model.Domain = {[-1 1],[-1 1],[-1 1],[-5 5]};
model.NonlinearDomain = [0;1;0;0];
model.mtd = 'Uniform';          % Uniform grid
model.NR =6;                     % Unifrom grid resolution

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
