%                        START OF MODEL DEFINITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear model; 

model.A=[-30 -20; .05  0];               %Enter the matrix A (2 x 2) here
model.aff=[24;  0];                     %Enter the affine term a (2 x 1) here
model.Bx=[20;   0];                     %Enter the matrix B(x) (2 x 1) as a function of x here
model.Domain={[-20 20],[-2e4 2e4]};     %Enter the domains of 2 variables here
model.fx=@resistor_nonlinearity;        %Enter the function handler here to call the function
                                        %that generates f(x) (2 x 1)
model.NonlinearDomain=[0 ;1];  %Change the elements whose corresponding variable
                              %appear in f(x) to one
                              %For example for a 2nd order system,
                              %if only x(1) appears in f(x), then model.
                              %NonlinearDomain=[1;0]
model.xcl=[0.371428571428570 ;
           0.642857142857146];                %Enter the desired equilibrium point here
model.mtd='Multiresolution';                %Enter the griddig type here:
                             %Choose from {'Uniform','Optimaluniform','Multiresolution'}
model.NR=3;                 %Enter the desired number of partitions
model.Rstar=[.2 .6];            %This is OPTIONAL and is in effect only if you choose 'Multiresolution' mode


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

