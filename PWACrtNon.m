clear model pwasys
echo on;
clc
% 
% *************************************************************************
%  CREATING THE MODEL (NONLINEAR)
% *************************************************************************
% 
% PWATOOL starts with importing a nonlinear or piecewise-affine model. If 
% the model is nonlinear PWATOOL approximates it with a PWA model.
% 
% References
% (1) B. Samadi and L. Rodrigues. Extension of local linear controllers to
%     global piecewise affine controllers for uncertain non-linear systems.
%     International Journal of Systems Science, 39(9):867-879, 2008.
% 
pause; %strike any key to continue
%
% Consider a nonlinear model 
%
%                     dx/dt=A x + a + f(x) + B(x) u 
%
% which consists of a linear dynamics A, an affine term a, a nonlinear
% dynamics f(x) and input gain B(x).
%
% In order to show how PWATOOL imports the above nonlinear model, we use an
% an example from (L. Rodrigues and S. Boyd, 2005); we model a nonlinear
% resistor and later approximate it by PWADIs.
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% | Creating a m-file for writing the data of the nonlinear model         |
% |-----------------------------------------------------------------------|
%
% The nonlinear resistor example from (L. Rodrigues and S. Boyd, 2005) is 
% of the order two (n=2) and has one input (m=1). 
%
% Therefore, to start creating the nonlinear model we type
%
%            pwacreate(2, 1, 'nonlinear_resistor_non.m')
%
pause; %strike any key to continue
%
% or in general 
%
%            pwacreate(n, m, 'myfilename.m')
%
% where myfilename.m is the m-file which is to store the nonlinear
% model. By doing so, 'myfilename.m' will be created (if does not exist
% before) or rewritten (if existed before).
%
% Now open 'nonlinear_resistor_non.m' which we created previously. We show
% how it has been filled.
%
pause; %strike any key to continue
clc
%
            open nonlinear_resistor_non.m
%            
%
% nonlinear_resistor_non.m has two sections of which, only, Section 1, 
% should be filled by the user. Section 2 should remain as is to let the
% code be executed after completion.
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% |  model.A : linear dynamics                                            |
% |-----------------------------------------------------------------------|
%
% model.A is a matrix of size (n x n). For the nonlinear resistor example
% we have 
%
           model.A=[-30 -20
                    .05  0];
%
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% |  model.aff : affine term                                                |
% |-----------------------------------------------------------------------|
%
% model.aff is column vector of size (n x 1). For the nonlinear resistor
% example we have 
%
           model.aff=[24
                       0];
%
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% |  model.Bx : input gain                                               |
% |-----------------------------------------------------------------------|
%
% model. Bx is a matrix of size (n x m). It can be a variable of class
% double or a function handle (please hit number 1 in the main menu of
% PWATOOL for more details). For the nonlinear example we have
%
           model.Bx=[20
                      0];
                  
%
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% |  model.Domain : domain of the variables                                               |
% |-----------------------------------------------------------------------|
%
% model.Domain is a cell of (1 x n). For the nonlinear resistor example we
% have
%
           model.Domain={[-20 20],[-2e4 2e4]};
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% |  model.fx : nonlinear dynamics                                               |
% |-----------------------------------------------------------------------|
%
% model.fx denotes the function handle which generates the nonlinearity in
% the model. For instance, for the nonlinear resistor example, we write 
%
            model.fx=@resistor_nonlinearity;
%
% to mean that 'resistor_nonlinearity.m' file contains the code which for
% each value of X=[x(1); x(2);...;x(n)] produces a column vector of size
% (n x 1).
%
pause; %strike any key to continue
%
% In our example 'resistor_nonlinearity.m' is coded as follows:
%
% 
%           ---------------------------------------
%          | function F = resistor_nonlinearity(X) |
%          | x=X(2);                               |
%          | if (-2e4 <= x & x< .2)                |  
%          |     F=[0; -50*(1e-3/.2*x)];           |  
%          | elseif (.2 <= x & x<.6)               |   
%          |     F=[0; -50*(-2e-3*x+1.4e-3)];      |     
%          | elseif (.6 <= x & x <= 2e4)           |  
%          |     F=[0; -50*(4e-3*x-2.2e-3)];       |   
%          | end                                   | 
%           ---------------------------------------
%
% Note that the output of the function should be a column vector of size 
% (n x 1). Hence some elements of vector 'F' above are set to zero since
% they are not set by the nonlinearity dynamics.
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% |  model.NonlinearDomain : domain of the nonlinear dynamics             |
% |-----------------------------------------------------------------------|
%
% model.NonlinearDomain is a row vector of size (1 x n) which shows which 
% variables appear in the nonlinear dynamics. For every variable x(j) which
% appear in the nonlinear dynamics model.fx we put a '1' at the j-th place
% of vector model.NonlinearDomain.
%
% For the nonlinear resistor example write
%
          model.NonlinearDomain=[0;1];
%
% to mean that model.fx is a function of only x(2).
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% |  model.xcl: equilibrium point                                         |
% |-----------------------------------------------------------------------|
% 
% model.xcl is a column vector of size (n x 1) and shows the "desired"
% equilibrium point which we need in analysis or synthesis. For the
% nonlinear resistor example we have
 
            model.xcl=[0.371428571428570
                       0.642857142857146];
%
%
pause; %strike any key to continue
clc
%
% |-----------------------------------------------------------------------|
% |  model.NR: number of regions                                          |
% |-----------------------------------------------------------------------|
%
% model.NR is a vector of size (1 x ND) where ND is the number of the
% variables which appear in model.fx. Elements of model.NR show the number
% of the regions which we like to have at each direction x(j) which appear
% in model.fx. 
%
% For the nonlinear resistor example we have
%
 
            model.NR=3;
% 
% If model.NR is given as a scalar, all of the variables x(j) which
% appear in model.fx will have the same number of regions in their
% corresponding direction.
%
pause; %strike any key to continue
clc
%
% |-----------------------------------------------------------------------|
% |  model.mtd: gridding type                                             |
% |-----------------------------------------------------------------------|
%
% model.mtd is a string which takes a value from 
% {'Uniform', 'OptimalUniform', 'Multiresolution'}. It shows the method
% that PWATOOL will use to grid the domain. For the nonlinear resistor
% exaplme we write
%
         model.mtd='Multiresolution';
%
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% |  model.Rstar: Rstar region                                            |
% |-----------------------------------------------------------------------|
%
% model.Rstar is matrix of size (* x ND) where ND is the number of
% variables which appear in model.fx and * denotes an optional number. For
% the nonlinear resistor example that ND=1 we write
%
           model.Rstar=[.2 ; .6];
%
%
% By specifying above Rstar and since the total number of the regions will
% be 3, the regions which we obtain will be R1=[-2e4 .2], R2=[.2 .6]
% and R3=[.6 2e4]
% 
%
%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                %       MODELING IS COMPLETE        % 
%                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% |  Checking and saving the data
% |-----------------------------------------------------------------------|
% 
% Now, save the changes in the file 'nonlinear_resistor_non.m' and run it
% by pressing F5 or typing 'run nonlinear_resistor_pwa'.
%
% PWATOOL imports the data that user has entered and checks for dimension
% consistency between the model parameters by calling the function
%
% [Err, model]=NonCheckModel(model);
%
pause; %strike any key to continue
 
              [Err, model]=NonCheckModel(model);
 
%               
pause; %strike any key to continue
%
% If there is no error (Err=0), function 'PWAComp' is called to build the
% PWADI approximation. 
%
%             [pwainc, pwasys, nlsys]=PWAComp(model);
%
%
pause; %strike any key to continue
 
              [pwainc, pwasys, nlsys]=PWAComp(model);
              
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% |  Outputs and results
% |-----------------------------------------------------------------------|
% 
% If PWATOOL accepts the data that user has entered, it computes the PWADI
% approximation for the nonlinear system and prints a message for the user.
%
% PWATOOL, generates three outputs. 
%
pause; %strike any key to continue
clc
% First, the original nonlinear system is saved in a variable called
% 'nlsys'. The following table shows the fields of 'nlsys' which, in fact,
% user has entered.
%
% |------------------------------------------------------------------------
% |PARAMETER             VARIABLE
% |------------------------------------------------------------------------
% |A                     model.A  
% |a                     model.aff
% |f(x)                  model.fx
% |B(x)                  model.Bx
% |xcl                   model.xcl
% |Domain                model.Domain
% |Gridding Type         model.mtd
% |Rstar                 model.Rstar
% |------------------------------------------------------------------------
%
pause; %strike any key to continue
clc
% Second, PWADI approximation is stored in a variable called 'pwainc'. A
% similar PWA approximation (without envelops) is also stored in a variable
% called 'pwasys'. The table below shows the fields in 'pwainc' or 'pwasys'
% 
% |------------------------------------------------------------------------
% |PARAMETER    STRUCTURE          REPRESENTRS             Size of Each 
% |                                                        Matrix / Vector
% |------------------------------------------------------------------------
% |A            NR x 2  cell       Linear dynamics         n x n
% |a            NR x 2  cell       affine term             n x 1
% |B            NR x 2  cell       Input gain              n x m
% |E            NR x 1  cell       Regions equation        * x n
% |e            NR x 1  cell       Regions equations       * x 1 
% |F            NR x NR cell       Boundary equations      n x (n-1) 
% |f            NR x NR cell       Boundary equations      n x 1 
% |Abar         NR x 2  cell       Linear dynamics         (n+1) x (n+1)
% |Bbar         NR x 2  cell       Input gain              (n+1) x m
% |Ebar         NR x 1  cell       Regions equation          *   x (n+1)
% |Fbar         NR x NR cell       Boundary equations      (n+1) x (n+1) 
% |NR           1  x ** array      Number of regions       1 x  **
% |------------------------------------------------------------------------
% *: The number of rows in matrix E and vector e depends on the gridding 
% type. If the system is slab, then [2, n]=size(E); if the gridding type is
% uniform and the system is not slab, then [n+1, n]=size(E) and
% [n+1, 1]=size(e).
%
% **: NR could be a scalar or a row vector of size 1 x ND where ND is the
% number of variables which appear in f(x).
% 
pause; %strike any key to continue
echo off
 


