% PWA Approximation

clear nlsys pwasys DesignParam


ANL = [0 1; 0 -.1];
BNL = [0; 1];

model.xcl = [0;0];
model.A = ANL;
model.Bx = BNL;
model.aff = [0;0];
model.fx = @NLFun_extension_6_2;
model.Domain = {[-4 4],[-2 2]};
model.NonlinearDomain = [1;0];
model.Resolution = 9;
model.mtd = 'Multiresolution'; % Multiresolution grid
model.NR =5;                    % Target number of regions
model.Rstar = [-2; -1; 1; 2];
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
    [pwainc, pwasys, model]=PWAComp(model);
    fprintf('\nNonlinear system saved in model.\n\n');
    fprintf('PWA Differential Inclusion approximation saved in pwainc.\n\n');
    fprintf('PWA approximation saved in pwasys.\n\n');
elseif Err==1
    fprintf('\nI quit because model has error.\n\n');
end
