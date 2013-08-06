
model.A = [0   1  0;
           0 -.01 0;
           0   0  0];
model.Bx = [0;1; 0];
model.aff = [0;0;0];
model.fx = @XY_Cart_Nonlinearity;
model.Domain = {[-3*pi/5 3*pi/5],[-20 20], [-10 10]};
model.NonlinearDomain = [1;0;0];
model.mtd = 'Uniform'; 
model.NR =4;                    
model.Rstar = [-pi/15; pi/15];
model.xcl = [0;0;0;];

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
    [pwainc, pwasys, model]=PWAComp(model);
    fprintf('\nNonlinear system saved in model.\n\n');
    fprintf('PWA Differential Inclusion approximation saved in pwainc.\n\n');
    fprintf('PWA approximation saved in pwasys.\n\n');
elseif Err==1
    fprintf('\nI quit because model has error.\n\n');
end