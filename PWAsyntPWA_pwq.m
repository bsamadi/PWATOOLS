clear model pwasys pwainc setting
echo on;
clc
% *************************************************************************
%  SYNTHESIZING CONTROLLERS FOR A PWA model (PWQ LYAPUNOV APPROACH)
% *************************************************************************
% 
% PWATOOL can synthesize PWA controllers to stabilize a PWA system around
% a desired equilibrium point xcl. 
%
% Let us explain how PWA controllers are synthesized by PWATOOL. 
%
% We use the "Local Controller Extension" method (see B. Samadi and
% L. Rodrigues, 2008) to design a local controller at the equilibrium point
% region and extend it later to the whole regions.
%
% Reference: B. Samadi and L. Rodrigues. Extension of local linear 
%            controllers to global piecewise affine controllers for uncertain
%            non-linear systems. International Journal of Systems Science,
%            39(9):867-879, 2008.
%
pause; %strike any key to continue
clc
%  SYNTHESIZING PWA CONTROLLERS; A PIECEWISE QUADRATIC (PWQ) LYAPUNOV APPROACH
%
%
% Consider a PWA model 
% 
%          dx/dt = A{i} x+ a{i}+ B{i} u   
%
% where regions R_i is defined as 
%
%             R_i = {x | E{i} x + e{i} >= 0 }.    i=1,2,...,NR.
% 
% Let the common boundary between regions R_i and R_h (if applicable) be
% contained in   
%
%             F{i,h} s + f{i,h}
%
% where s ranges over the vector space R^{n-1}
%
pause; %strike any key to continue
%
% We use Lyapunov functions to ensure the stability of the closed-loop
% PWADIs. 
%
pause; %strike any key to continue
clc
% By PWQ Lyapunov approach in synthesizing PWA controllers u = K{i} x + k{i},
% we mean that matrices P{i}, vectors q{i} and scalars r{i} as well as the
% controller gains K{i} and k{i} should be found such that the 
% following two conditions hold:
%
%    1)  V(x)  = x'*P{i}*x + 2*q{i}'*x +r{i}   >0         
%    2)  dV/dt =   <-alpha * V(x)                        
%
% Here 'alpha' is the decay rate constant that should be set by the user in
% advance. In a sense, the higher alpha is, the faster the closed-loop 
% dynamics is.
%                                                   
pause; %strike any key to continue
clc
% Also let Z{i} and W{i} be matrices of proper sizes with non-negative
% elements which are to be found. Then, conditions (1) and (2) amount to
% the following problem, formulated as BMI:
%
pause; %strike any key to continue
%
% PROBLEM: Find P{i}, q{i}, r{i}, K{i}, k{i}, Z{i} and W{i} for i=1,2,...,NR
% such that we have
% 
% (1) In regions R_i which contain the equilibrium point xcl
%
%   I:  P{i}*(A{i}+B{i}*K{i})+(A{i}+B{i}*K{i})'*P{i}+alpha*P{i}  < 0
%    
%  II:  [P{i}       q{i}
%        q{i})      r{i}] > 0
%
%
pause; %strike any key to continue
%
% (2) In regions R_i which do NOT contain the equilibrium point xcl
%
%   I: [P{i}*(A{i}+B{i}*K{i})+(A{i}+B{i}*K{i})'*P{i}+alpha*P{i}+E{i}'*Z{i}*E{i}     P{i}*(a{i}+B{i}*k{i})+E{i}'*Z{i}*e{i}+(A{i}+B{i}*K{i})'*q{i}
%      (P{i}*(a{i}+B{i}*k{i})+E{i}'*Z{i}*e{i}+(A{i}+B{i}*K{i})'*q{i})'                                  e{i}'*Z{i}*e{i}+2*q{i}'*(a{i}+B{i}*k{i})]  < 0
%
%
%  II:  [P{i}-E{i}'*W{i}*E{i}          -E{i}'*W{i}*e{i}+q{i}
%        (-E{i}'*W{i}*e{i}+q{i})'      -e{i}'*W{i}*e{i}+r{i}] > 0
%    
%
pause; %strike any key to continue
% 
% (3) Continuity of the Lyapunov function
% 
%      [F{i,h}'*(P{i}-P{h}*F{i,h}            F{i,h}'*(P{i}-P{h})*f{i,h}
%       (F{i,h}'*(P{i}-P{h})*f{i,h})'        f{i,h}'*(P{i}-P{h})*f{i,h}] = 0      for F{i,h} not empty 
%
pause; %strike any key to continue
%
% (4) Continuity of the control outputs
% 
%     ((A{i}+B{i}*K{i})-(A{h}+B{h}*K{h})) * F{i,h}=0       for F{i,h} not empty
%     ((a{i}+B{i}*k{i})-(a{h}+B{h}*k{h})) * f{i,h}=0
%
%
pause; %strike any key to continue
clc
% NOTE 1: PWATOOL uses Yalmip to solve above BMIs. A solution, then, may or
%         may not be found by the solver, depending on the size of the 
%         problem. 
%
pause; %strike any key to continue
%
% NOTE 3: The regions can also be approximated with ellipsoids which in the
%         case of slab regions they will be degenerate ellipsoids. This 
%         yields to a different set of BMIs.
%
% References: (1) L. Rodrigues and S. Boyd. Piecewise-affine state feedback 
%             for piecewise-affine slab systems using convex optimization. 
%             Systems and Control Letters, 54:835-853, 200
%             
%             (2) A. Hassibi and S. Boyd. Quadratic stabilization and
%             control of piecewise-linear systems. Proceedings of the
%             American Control Conference, 6:3659-3664, 1998
% 
%             (3) L. Vandenberghe, S. Boyd, and S.-P. Wu. Determinant
%             maximization with linear matrix inequality constraints. SIAM
%             Journal on Matrix Analysis and Applications,9(2):499-533,1998
%
pause; %strike any key to continue
%
% To obtain the controllers with ellipsoidal approximation, let 'Elip_i' be
% the ellipsoid that approximates region R_i:
% 
%                Elip_i={x | || EL{i} x + eL{i}|| < 1}
% 
pause; %strike any key to continue
clc
% Then, for i=1,2,...,NR the following problem, formulated as BMI, gives
% sufficient conditions for stabilizing the PWA model around xcl.
%
% PROBLEM: Find P{i}, q{i}, r{i} and scalars miu{i} <0 and bita{i}<0 such
% that the constraints below hold for i=1,2,...,NR and j=1,2:
%
pause; %strike any key to continue
%
% (1) In regions R_i which contain the equilibrium point xcl
% 
%
%   I:  P{i}*(A{i}+B{i}*K{i})+(A{i}+B{i}*K{i})'*P{i}+alpha*P{i}  < 0
%    
%  II:  [P{i}       q{i}
%        q{i})      r{i}] > 0
%
%
pause; %strike any key to continue
%
% (2) In regions R_i which do NOT contain the equilibrium point xcl
%
%   I:  [Q*(A{i}+B{i}*K{i})+(A{i}+B{i}*K{i})'*Q+alpha*Q+miu{i}*EL{i}'*EL{i}     Q*(a{i}+B{i}*k{i})+miu{i}*EL{i}'*eL{i}+(A{i}+B{i}*K{i})'*q{i}
%       (Q*(a{i}+B{i}*k{i})+miu{i}*EL{i}'*eL{i}+(A{i}+B{i}*K{i})'*q{i})'                   -miu{i}*(1-eL{i}'*eL{i})+2*q{i}'*(a{i}+B{i}*k{i})] < 0
%
%  II:  [P{i}-bita{i}*E{i}'*E{i}               -bita{i}*E{i}'*e{i}+q{i}
%        (-bita{i}*E{i}'*e{i}+q{i})'        bita{i}*(1-e{i}'*e{i})+r{i}] > 0
%
%
pause; %strike any key to continue
% 
% (3) Continuity of the Lyapunov function
% 
%      [F{i,h}'*(P{i}-P{h}*F{i,h}            F{i,h}'*(P{i}-P{h})*f{i,h}
%       (F{i,h}'*(P{i}-P{h})*f{i,h})'        f{i,h}'*(P{i}-P{h})*f{i,h}] = 0      for F{i,h} not empty 
%
pause; %strike any key to continue
%
% (4) Continuity of the control outputs
% 
%     ((A{i}+B{i}*K{i})-(A{h}+B{h}*K{h})) * F{i,h}=0       for F{i,h} not empty
%     ((a{i}+B{i}*k{i})-(a{h}+B{h}*k{h})) * f{i,h}=0
%
%
pause; %strike any key to continue
clc
%             USING POWTOOL FOR SYNTHESIZING PWA CONTROLLERS
% 
% We use the Cart example (L. Rodrigues and S. Boyd, 2005) to show that 
% how PWATOOL synthesizes PWA controllers for a PWA model using a PWQ
% Lyapunov function.
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% | STEP 0: APPROXIMATING THE NONLINEAR MODEL WITH PWADI
% |-----------------------------------------------------------------------|
% 
% Please hit number 2 in the main menu of PWATOOL for "CREATING A PWA MODEL"
%
% For the Cart example we load the variable 'cart_pwa'
% 
        load cart_pwa
%
% by doing so, a variable called 'pwasys' will contain the PWA model of the
% Cart.
% 
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% | STEP 1: SYNTHESIZING PWA CONTROLLERS
% |-----------------------------------------------------------------------|
%
% PWA controllers can be synthesized by PWATOOL by the following command:
% 
% ctrl = pwasynth(pwasys, x0, setting)
%
% where 'pwasys' is the PWADI approximation of the nonlinear model, x0 is
% the initial point for simulation and 'setting' an array structure.
%
% We, first, choose the initial point
%
      x0=[1.2; 3; 2];
%
% Note that most of the fields of 'setting', if not set by the user, will
% take default values by PWATOOL. We explain the fields of 'setting' in the
% following.
%
pause; %strike any key to continue
clc
% setting.Lyapunov: shows whether 'global' or 'pwq' (piecewise quadratic)
%                   Lyapunov function should be used in the synthesis
%
% For the Cart example we set 
%
      setting.Lyapunov='pwq';
%
% This means that PWATOOL uses only global Lyapunov function in controller
% synthesis.
%
pause; %strike any key to continue
clc
% setting.ApxMeth:  shows whetehr 'ellipsoidal' or 'quadratic' curve
%                   method should be used to approximate the regions.
%
% For the Cart example we set 
%
       setting.ApxMeth='quadratic';  
% 
% so that only quadratic curve method is used for approximating the
% regions.
%
pause; %strike any key to continue
clc
% setting.SynthMeth: shows whether 'LMI' or 'BMI' methods should be used.
%
% For the Cart example we set
%
      setting.SynthMeth='bmi';
%
pause; %strike any key to continue
clc
% setting.QLin  :  positive definite matrix of size n x n where n is the 
%                  order of the system. It is used as the LQR parameter
%                  for the region which contains the equilibrium point
%
% The Cart example is of degree 3. We set QLin as follows
%
       setting.QLin = [ 14.9748    7.4636    9.6918
                         7.4636    4.0334    5.1840
                         9.6918    5.1840   10.7943];
% 
%    
pause; %strike any key to continue
clc
% setting.RLin  :  positive definite matrix of size m x m where m is the 
%                  number of the inputs to the system. It is used as the 
%                  LQR parameter for the region which contains the 
%                  equilibrium point
%
% The Cart example has one input. Therefore, we set 
%
       setting.RLin = [5.3309];
% 
%
pause; %strike any key to continue
clc
% setting.RandomQ: if set, shows setting.QLin is set randomly by PWATOOL.
%                  If it is zero, PWATOOL keeps the value that user enters
%                  for setting.QLin
%
% For the Cart example we reset this parameters:
%
     setting.RandomQ=0;
%
pause; %strike any key to continue
%
% setting.RandomR: if set, shows setting.RLin is set randomly by PWATOOL.
%                  If it is zero, PWATOOL keeps the value that user enters
%                  for setting.RLin
%
% For the Cart example we reset this parameters:
%
      setting.RandomR=0;
%
pause; %strike any key to continue
clc
% setting.alpha :  specifies the decay rate of the closed-loop system.
%                  Generally, the bigger alpha is the faster the dynamics
%                  of the system will be. However, bigger values of alpha
%                  make the convergence of the controller design harder.
%
% For the Cart example dynamics, we set
%
      setting.alpha=.1;
%
pause; %strike any key to continue
clc
% setting.NormalDirectionOnly : specifies if only the normal directions at
%                  the boundaries should meet the continuity constraints.
%                  By default, it would be zero so that the whole
%                  directions are included in the continuity constraints.
%                  In that case, if convergence is obtained, the results are
%                  usually smoother than those obtained by
%                  NormalDirectionOnly set to 1. Please see the PWATOOL
%                  document for a detailed explanation.
%
%
pause; %strike any key to continue
clc
% setting.xcl   :  shows the desired equilibrium point which user sets.
%                  If it is not set by the user, PWATOOL uses the
%                  pwainc.xcl (if exists) as the equilibrium point.
%                  PWATOOL issues an error and stops running if at xcl
%                  we cannot satisfy the equilibrium point equations. 
%                  analytically.
%
pause; %strike any key to continue
clc
% setting.StopTime: shows the simulation time in seconds. If not set by the
%                   user, it will be set to 10 sec.
%
pause; %strike any key to continue
%
% setting.IterationNumber: shows the maximum number of times that PWATOOL 
%                  tries to synthesize PWA controller if it fails in doing
%                  so in the first try. The default value is 5.
%
%
pause; %strike any key to continue
clc
% Table below shows a summary of setting's fields and their default values
%
% |------------------------------------------------------------------------
% |FIELD            STRUCTURE       SIZE             DEFAULT VALUE
% |------------------------------------------------------------------------
% |ApxMeth          Cell           (1 x 2)           {'ellipsoidal', 'quadratic'}
% |Lyapunov         Cell           (1 x 2)           {'global', 'pwq'}
% |SynthMeth        Cell           (1 x 2)           {'lmi', 'bmi'}
% |QLin             Array          (n x n)           Random QLin>0
% |RLin             Array          (m x m)           Random RLin>0
% |RandomQ          scalar         (1 x 1)           1
% |RandomR          scalar         (1 x 1)           1
% |alpha            scalar         (1 x 1)           0.1
% |xcl              Array          (n x 1)           pwasys.xcl (if exists)
% |StopTime         scalar         (1 x 1)           10 sec
% |IterationNumber  scalar         (1 x 1)           5
% |------------------------------------------------------------------------
%
pause; %strike any key to continue
clc
% |-----------------------------------------------------------------------|
% | SSTEP 2: RESULTS AND DUISCUSSIONS
% |-----------------------------------------------------------------------|
%
% Based on the settings that user has chosen for the fields 'Lyapunov', 
% 'ApxMeth' and 'SynthMeth', PWATOOL synthesizes PWA controllers by
% different methods. The current options are:
%
% 1) LMI approach with ellipsoidal approximation     (golabl Lyapunov)
% 2) BMI approach with ellipsoidal approximation     (golabl Lyapunov)
% 3) BMI approach with quadratic curve approximation (golabl Lyapunov)
% 4) BMI approach with ellipsoidal approximation     (pwq Lyapunov)
% 5) BMI approach with quadratic curve approximation (pwq Lyapunov)
%
pause; %strike any key to continue
clc
% A variable named 'pwacont' stores the result of the synthesis if at
% least one of the methods among the possibilities converge. In this case
% another variable 'pwasetting' stores the setting that PWATOOL has used
% for the synthesis.
%
% However, if none of the methods, among the possibilities, converge in the
% first try, PWATOOL tries to a maximum of 'IterationNumber' times which by
% default is set to 5. 
%
pause; %strike any key to continue
clc
% Now let's run pwasynth to synthesize PWA controllers for the nonlinear
% resistor model. We type
% 
% ctrl=pwasynth(pwasys, x0, setting);
% 
pause; %strike any key to continue
%
ctrl=pwasynth(pwasys, x0, setting);
 
pause; %strike any key to continue
clc
% The varable called 'pwacont' contains the PWA controllers:
%
    pwacont
%
pause; %strike any key to continue
%
% 'pwacont' will contain several fields where each field is itself an 
% array cell, containing the information about he whole methods which have
% converged.
%
% We explain the fields of 'pwacont' in the following:
%
pause; %strike any key to continue
clc
% pwacont.Gian : contains the controller which has been designed (Kbar), in
%                the aggregated vector space. For example pwacont.Gain{1}
%                will be a cell of size (1 x n) where each element would
%                contain the controller Kbar{i}=[K{i} k{i}] for region R_i
%                (i=1,2,...,NR). 
%
pause; %strike any key to continue
clc
% pwacont.AppMethod:    shows the corresponding approximation method for
%                       pwacont.Gain which would be 'ellipsoidal' or
%                       'quadratic curve'.
%
pause; %strike any key to continue
clc
% pwacont.SynthesisType: shows the corresponding synthesis type for the Gain 
%                        which would be LMI or BMI
%
pause; %strike any key to continue
clc
% pwacont.Lyapunov :     shows the corresponding Lyapunov function type for
%                        the Gain, which would be 'global' or 'pwq'.
%
pause; %strike any key to continue
clc
% pwacont.Sprocedure1 and 
% pwacont.Sprocedure2  : shows the corresponding S-procedure coefficients
%                        which PWATOOL has used for each method.
%
pause; %strike any key to continue
clc
% Table below shows a summary of 'pwacont' fields:
%
% |------------------------------------------------------------------------
% |FIELD             STRUCTURE       SIZE          COMMENT
% |------------------------------------------------------------------------
% |Gain              cell           (1 x *)        controller gain
% |AppMethod         cell           (1 x *)        approximation method
% |SynthesiMethod    cell           (1 x *)        synthesis method
% |Lyapunov          cell           (1 x *)        Lyapunov type
% |------------------------------------------------------------------------
%
% *: The size of the cells depend on the number of the methods that
%    converge.
%
pause; %strike any key to continue
clc
pwacont
% we see that one method has converged:
%
% 1) BMI approach with quadratic curve approximation (pwq Lyapunov)
% 
% This information is obtained from the first four fields of 'pwacont'.
%
pause; %strike any key to continue
clc
% 'pwacont' shows the quadratic approximation, by BMI approach and a 
% PWQ Lyapunov function has yielded a PWA controller for regions R_i
% i=1,2,...,NR in the aggregated vector space y=[x;1] as follows:
%
ctrl.Gain{1}{:}
%
pause; %strike any key to return to the main menu
echo off
 


