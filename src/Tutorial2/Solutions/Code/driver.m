%% Question 1

% Clear
clear all, close all, clc

% Drawing semin-analytic bifurcation diagram
bd = figure();
fimplicit( @(lam,u) lam.*u + u.^3 - u.^5,[-0.5 1 -1.5 1.5]);
xlabel('lambda'); ylabel('C');

% Preparing vector of lambda values (and norm)
lambdaVals = [0.7 -0.4];
l2NormVals = zeros(size(lambdaVals));

% Setup differentiation matrices
nx = 100
[x,~,Dxx] = NeumannDiffMat([-5,5],nx);

% For each value of lambda
for i = 1:length(lambdaVals) 

  % Set parameters and timestep
  p = [1; lambdaVals(i); 0; 1; 1];
  u0 = 1./cosh(x).^2;
  rhs = @(t,u) AllenCahn(u,p,Dxx);
  tSpan = [0 50];
  [t,UHist] = ode15s(rhs,tSpan,u0);
  uFinal = UHist(end,:)';

  % Compute and display l2 norm of final state
  l2NormVals(i) = ComputeL2Norm(uFinal,x);

end
figure(bd); hold on; plot(lambdaVals,l2NormVals,'*'); hold off;

%% Question 2
% The continuation follows the upper branch of steady states until the saddle-node
% bifurcation is reached. To the left of the saddle-node bifurcation, the only
% branch of stable steady states is the trivial state $u(x) \equiv 0$, hence the
% continuation jumps to that branch and sticks with it. In general, one should expect
% that, when the branch has an instability, the brute-force continuation will jump to
% a nearby stable branch. This is because, via time simulation, the only attainable
% steady states are the linearly stable ones: when one initialises a time simulation
% from an unstable steady state (or nearby) the solution departs from it, by
% definition of linear instability. As we shall see, however, this must be taken with
% a pinch of salt: sometimes the dynamics stays nearby an unstable equilbirium for
% long times before departing, and this can lead one to think that an equilibrium
% found via direct simulation is stable, when in fact it is unstable. There is no way
% to asess stability with brute-force simulation of the original PDE.

% Values of lambda and l2 norms
lambdaVals = 0.7:-0.02:-0.5;
l2NormVals = zeros(size(lambdaVals));
resVals   = zeros(size(lambdaVals));

disp('Branch of homogeneous steady states');
fprintf('----------------------------------------\n');
fprintf('%12s %12s %12s\n','lambda','l2 Norm', '||F(U)||');
fprintf('----------------------------------------\n');

for i = 1:length(lambdaVals) 

  % Time Step
  p = [1; lambdaVals(i); 0; 1; 1];
  if i == 1
    u0 = 1./cosh(x).^2;
  else
    %pause
    u0 = uFinal;
  end
  rhs = @(t,u) AllenCahn(u,p,Dxx);
  tSpan = [0 50];
  [t,UHist] = ode15s(rhs,tSpan,u0);

  % Plot initial and final state
  uFinal = UHist(end,:)';

  % Compute and display l2 norm of final state
  l2NormVals(i) = ComputeL2Norm(uFinal,x);
  resVals(i) = norm(AllenCahn(uFinal,p,Dxx),'inf');
  fprintf('%12.4e %12.4e %12.4e \n',lambdaVals(i),l2NormVals(i),resVals(i));

end
fprintf('----------------------------------------\n');

figure(bd), hold on;
plot(lambdaVals,l2NormVals,'*');
hold off;

%% Question 3
% A branch of "bump" solutions exixsts. The branch is stable until it undergoes a
% saddle-node bifurcation for $\lambda \in [-0.14, -0.12]$. We also observe a jump
% towards the trivial homogeneous steady state $u(x) \equiv 0$.

% Values of lambda and l2 norms
lambdaVals = 0.7:-0.02:-0.5;
l2NormVals = zeros(size(lambdaVals));
resVals   = zeros(size(lambdaVals));

disp('Branch of patterned steady states');
fprintf('----------------------------------------\n');
fprintf('%12s %12s %12s\n','lambda','l2 Norm', '||F(U)||');
fprintf('----------------------------------------\n');

for i = 1:length(lambdaVals) 

  % Time Step
  p = [1; lambdaVals(i); 0; 1; 1];
  if i == 1
    u0 = 0.4*cos(2*pi/10*x);
  else
    u0 = uFinal;
  end
  rhs = @(t,u) AllenCahn(u,p,Dxx);
  tSpan = [0 50];
  [t,UHist] = ode15s(rhs,tSpan,u0);

  % Plot initial and final state
  uFinal = UHist(end,:)';
  uNorm = norm(AllenCahn(uFinal,p,Dxx),'inf');

  if i == 1
    figure(); plot(x,uFinal); xlabel('x'); ylabel('u(x)');
    title(['Equlibrium for lambda = ' num2str(lambdaVals(i))...
      ', ||F(u)|| = ', num2str(uNorm)]);
  end


  % Compute and display l2 norm of final state
  l2NormVals(i) = ComputeL2Norm(uFinal,x);
  resVals(i) = norm(AllenCahn(uFinal,p,Dxx),'inf');
  fprintf('%12.4e %12.4e %12.4e \n',lambdaVals(i),l2NormVals(i),resVals(i));

end
fprintf('----------------------------------------\n');

figure(bd), hold on;
plot(lambdaVals,l2NormVals,'*');
hold off;


%%  Question 4

% Time step with lambda = 0.1
p = [1; 0.7; 0; 1; 1];
u0 = 0.5*tanh(-x);
rhs = @(t,u) AllenCahn(u,p,Dxx);
tSpan = [0 10];
[t,UHist] = ode15s(rhs,tSpan,u0);

% Plot
[X,T] = meshgrid(x,t);
figure; surf(X,T,UHist); shading interp;
xlabel('x'); ylabel('t'); zlabel('u');

uFinal = UHist(end,:)';
uNorm = norm(AllenCahn(uFinal,p,Dxx),'inf');

figure; plot(x,uFinal); xlabel('x'); ylabel('u(x)');
title(['Equlibrium for lambda = ' num2str(lambdaVals(i)) ...
  ', ||F(u)|| = ', num2str(uNorm)]);

%% Question 5
% The branch becomes unstable at a bifurcation (whose type can not yet be inferred)
% in the interval $\lambda = [0.36,0.38]$. When the instability occurs, the solution
% jumps to the stable branch of homogeneous steady states. It may seem
% counterintuitive, at first, that the solution does not jump to the closer branch of
% patterned bump states (in $S_2$ norm). The fact that the branch appears closer
% than the homogeneous one depends clearly on the norm we selected, and it is
% therefore a projection effect.

% Values of lambda and l2 norms
lambdaVals = 0.7:-0.02:-0.5;
l2NormVals = zeros(size(lambdaVals));
resVals   = zeros(size(lambdaVals));

disp('Branch of patterned steady states');
fprintf('----------------------------------------\n');
fprintf('%12s %12s %12s\n','lambda','l2 Norm', '||F(U)||');
fprintf('----------------------------------------\n');

for i = 1:length(lambdaVals) 

  % Time Step
  p = [1; lambdaVals(i); 0; 1; 1];
  if i == 1
    u0 = 0.5*tanh(-x);
  else
    %pause
    u0 = uFinal;
  end
  rhs = @(t,u) AllenCahn(u,p,Dxx);
  tSpan = [0 50];
  [t,UHist] = ode15s(rhs,tSpan,u0);

  % Plot initial and final state
  uFinal = UHist(end,:)';

  % Compute and display l2 norm of final state
  l2NormVals(i) = ComputeL2Norm(uFinal,x);
  resVals(i) = norm(AllenCahn(uFinal,p,Dxx),'inf');
  fprintf('%12.4e %12.4e %12.4e \n',lambdaVals(i),l2NormVals(i),resVals(i));

end
fprintf('----------------------------------------\n');

figure(bd), hold on;
plot(lambdaVals,l2NormVals,'*');
hold off;

%% Question  6
% Spatial grid
nx = 100;
[x,~,Dxx] = NeumannDiffMat([-5,5],nx);

% For all flags
flagList = {'flat','bump','sigmoid'};
for k = 1:length(flagList)

  % Get flag and setup initial condition
  u0Flag = flagList{k}
  switch u0Flag
    case 'flat'
      u0 = zeros(size(x))+0.01*rand(size(x));
    case 'bump'
      u0 = 0.4*cos(pi/5*x);
    case 'sigmoid'
      u0 = 0.5*tanh(x);
    otherwise 
      error('u0Flag not implemented');
  end

  % Setup problem and time step
  p = [1; 0.7; 0; 1; 1];
  rhs = @(t,u) AllenCahn(u,p,Dxx);
  tSpan = [0 50];
  [t,UHist] = ode15s(rhs,tSpan,u0);

  % Plot surface
  [X,T] = meshgrid(x,t);
  figure; 
  subplot(1,2,1);
  surf(X,T,UHist); shading interp;
  xlabel('x'); ylabel('t'); zlabel('u');

  % Get fina state and estimate its time variation
  uFinal = UHist(end,:)';
  uNorm = norm(AllenCahn(uFinal,p,Dxx),'inf');

  % Plot final state
  subplot(1,2,2);
  plot(x,u0,x,uFinal); xlabel('x'); ylabel('u(x)');
  legend({'u(0)','u(T)'});
  title(['Equlibrium for lambda = ' num2str(p(2)) ...
    ', ||F(u)|| = ', num2str(uNorm)]);

  % Save a file (flag dependent)
  fileName = ['solution_' u0Flag '.mat'];
  u = uFinal;
  save(fileName,'u','p');

end

%% Question 7
% For 7.1 dn 7.2 see AllenCahn.m
 
% Setup problem
sol = load('solution_sigmoid.mat');
u0 = sol.u; p = sol.p;
prob = @(u) AllenCahn(u,p,Dxx);

% Setup and solve BVP
opts = optimset('Display','iter',...
                'TolFun',1e-10,...
                'MaxIter',50,...
                'Jacobian','on');
u = fsolve(prob,u0,opts);

% Plot solution
lims = [-5 5 -2 2];
eqFig = PlotSteadyState(x,u,[],lims,'',[]);

%% Question 8

% Perturb parameter
p(2) = p(2) - 0.3;
prob = @(u) AllenCahn(u,p,Dxx);

% Setup and solve BVP
opts = optimset('Display','iter',...
                'TolFun',1e-10,...
                'MaxIter',50,...
                'Jacobian','on');
u = fsolve(prob,u0,opts);

% Plot solution
eqFig = PlotSteadyState(x,u,eqFig,[-5 5 -2 2],'',false); hold off;

%% Question 9

% Grid and differentiation matrix
[x,~,Dxx] = NeumannDiffMat([-5,5],nx);

% Initial guess and parameters
sol = load('solution_sigmoid.mat'); u0 = sol.u; p0 = sol.p;

% Problem setup
prob = @(u,p) AllenCahn(u,p,Dxx);

% Continuation parameters
icp    = 2;    % Continuation parameter index. icp = 2 means continuation in p(2)
ds     = -0.1; % Parameter stepsize (negative means initially decreasing parameter)
nSteps = 200;  % Number of continuation steps to be taken
opts = optimset('Display','off',... % Options to pass to fsolve
                'TolFun',1e-10,...
                'MaxIter',50,...
                'Jacobian','on');

% Launch continuation and plot bifurcation diagaram
[bd1,sol1] = SecantContinuation(prob,u0,p0,icp,ds,nSteps,opts);
bdFig = figure; plot(bd1(:,2),  bd1(:,3), '.-','DisplayName','sigmoidal');
ylim([-1 14]);
xlabel('lambda'); ylabel('2-Norm');
legend show


%% Question 10

% Load initial guess and parameters
sol = load('solution_bump.mat'); u0 = sol.u; p0 = sol.p;
prob = @(u,p) AllenCahn(u,p,Dxx);

% Continuation 
icp    = 2;   
ds     = -0.1;
nSteps = 100; 
opts = optimset('Display','off', 'TolFun',1e-10,'MaxIter',50,'Jacobian','on');
[bd4,sol4] = SecantContinuation(prob,u0,p0,icp,ds,nSteps,opts);
figure(bdFig), hold on;
plot(bd4(:,2),  bd4(:,3), '.-','DisplayName','bump');
lgd = legend('show'); lgd.Location = 'northwest';

%% Question 11
% The solution presented below uses way too many resources. In principle, we can
% continue homogeneous steady states by continuing in $\lambda$ solutions to the
% scalar problem $0=G(\lambda,C)$, where $C$ is a real number. This would require
% defining a new problem. Here, mostly because I am being lazy and $n$ is small, I
% find the branch by continuing the $n$-dimensional discretised boundary value
% problem $0 = D_xx U + N(U,\lambda)$.

% Continuation of homogeneous equilibrium C ~= 0 
sol = load('solution_flat.mat'); u0 = sol.u; p0 = sol.p;
prob = @(u,p) AllenCahn(u,p,Dxx);

% Continuation
icp    = 2;    
ds     = -0.1;
nSteps = 200; 
opts = optimset('Display','off', 'TolFun',1e-10, 'MaxIter',50,'Jacobian','on');
[bd2,sol2] = SecantContinuation(prob,u0,p0,icp,ds,nSteps,opts);
figure(bdFig), hold on;
plot(bd2(:,2),  bd2(:,3), '.-','DisplayName','u_* = C, C > 0 ');
legend show

% Continuation of homogeneous equilibrium C = 0 
u0 = zeros(size(x)); p0 = sol.p;
prob = @(u,p) AllenCahn(u,p,Dxx);

% Continuation 
icp    = 2;  
ds     = -0.05;
nSteps = 30;  
opts = optimset('Display','off', 'TolFun',1e-10, 'MaxIter',50, 'Jacobian','on');
[bd3,sol3] = SecantContinuation(prob,u0,p0,icp,ds,nSteps,opts);
figure(bdFig), hold on;
plot(bd3(:,2),  bd3(:,3), '.-','DisplayName','u_*=0');
legend show

%% Question 12
% See accompanying pdf

%% Question 13

% Set value of lambda and u around which we linearise
lambda = 0; p0(2) = lambda; u = 0*x;

% Get Jacobian and compute spectrum
[~,J] = AllenCahn(u,p0,Dxx);
[V,D] = eig(full(J)); % Note, this is a bad idea for large matrices! For larger
                      % problems one should use eigs instead.
% Sort the eigenvalues in descending order.
[d,ix] = sort(diag(D),'descend'); 

figure;
for j = 1:5
  subplot(5,2,2*(j-1)+[1 2]);
  plot(x,real(V(:,ix(j))),'DisplayName','Re \psi(x)'); hold on;
  plot(x,imag(V(:,ix(j))),'DisplayName','Im \psi(x)'); hold off;
  title(['lambda = ' num2str(real(d(j))) ' + i ' num2str(imag(d(j))) ]);
  xlabel('x'); ylim([-0.3 0.3]);
  lgd = legend;
  lgd.Location = 'eastoutside';
end

% I fuond it handy to write two functions, ComputeStability, and
% ExploreBifurcationDiagram which compute and display useful information using the
% branch/solution structure outputted by SecantContinuation. So I create here a
% branch structure with just 3 solutions, and use them to provide evidence as
% required in this question
bd0 = [0 0 0;...
       0 (pi/10)^2 0;...
       0 (2*pi/10)^2 0];
sol0 = [0*x'; 0*x'; 0*x'];

%% Question 14  
% We now use this framework for analysing the solutions at $\lambda = 0, (\pi/10)^2,
% (2\pi/10)^2$. When $\lambda = 0$, we see one critical eigenvalue, with a critical
% eigenmode that is constant in space $\psi_0(x) \equiv 0$

bd = bd0; sol = sol0; id = 1;
[V,D] = ComputeStability(bd,sol,p,Dxx,id);
ExploreBifDiag(bd,sol,id,D,V,x);

%% 
% When $\lambda = (\pi/10)^2$, we see a second critical eigenvalue (the former
% critical eigenvalue has now become an unstable eigenvalue), with a critical
% eigenmode $\psi_1(x) = c_1 \sin( \pi x /10)$ (here $c_1 > 0$). From this point a
% he branch of sigmoidal equilibria emanate. In fact, sigmoidal profiles with small amplitude
% (close to the bifurcation) will look similar to the function $u_*(x) + \psi_1(x) =
% \psi_1(x)$. Here we also note that the critical eigenvalue is not exactly
% $0$, and this should be expected because we are looking at eigenvalues of the
% differentiation matrix with $n=100$, not the linear operator.

bd = bd0; sol = sol0; id = 2;
[V,D] = ComputeStability(bd,sol,p,Dxx,id);
ExploreBifDiag(bd,sol,id,D,V,x);

%% 
% When $\lambda = (2\pi/10)^2$, we see a third critical eigenvalue (the former
% critical eigenvalues have now become an unstable eigenvalues), with a critical
% eigenmode $\psi_2(x) = c_2 \sin(\pi x /5)$ (here $c_2 < 0$). Here the emanating
% branch originates "bumps".

bd = bd0; sol = sol0; id = 3;
[V,D] = ComputeStability(bd,sol,p,Dxx,id);
ExploreBifDiag(bd,sol,id,D,V,x);


%% Question 15 
% These are a few considerations on the first 2 open questions. We compute the first
% 5 eigenvalues for all the solutions along the branch of bumps. We then plot $\mu_i(s)$ and
% $\lambda(s)$, where $s$ is a coordinate parametrising the solution branch. One can
% see along the branch there are always at least 2 unstable eigenvalues. At the
% saddle node bifurcation we see one real eigenvalue crossing at nonzero speed (see
% red star on both graphs). When looking at the eigenvalues plot, recall that the
% continuation is from $\lambda = 0.7$, for decreasing $\lambda$. We also report the
% value of $\lambda$ along the branch in the same plot as the eigenvalues, so as to
% have a reference point for the bifurcations. Overall, we conclude that this branch
% is unstable for all values of $\lambda \leq 0.7$.

% Computing eigenvalues along the branch of bumps
bd = bd4; sol = sol4;
nSol = size(bd,1);
eVals = zeros(nSol,5);
for i = 1:size(bd,1)
  id = i;
  [~,D] = ComputeStability(bd,sol,p0,Dxx,id);
  [d,ix] = sort(diag(D),'descend'); 
  eVals(i,:) = d(1:5);
end

% Plotting branch
figure;
subplot(2,2,[1 2]);
plot(bd(:,2),bd(:,3)); 
hold on; ii = 41;
plot(bd(ii,2),bd(ii,3),'r*'); 
xlabel('\lambda'); ylabel('S_0');
title('Branch of bump solutions');

% Plotting eigenvalues
subplot(2,2,[3 4]);
plot(1+bd(:,1),real(eVals),'k-', 1+bd(:,1), bd(:,2),'-'); grid on; 
title('Most unstable eigenvalues along the branch');
xlabel('s'); ylabel('Re \mu_i');
lgd = legend({'\mu_1','\mu_2','\mu_3','\mu_4','\mu_5','\lambda'});
lgd.Location = 'southeast';
hold on;
plot(ii, eVals(ii,3),'r*');
plot(ii, bd(ii,2),'r*');
hold off;

%% 
% A similar statement can be made for the branch of sigmoidal steady states, which
% are unstable for all $\lambda \leq 0.7$. We note that there are several
% bifurcations on this branch (and we should expect branches of solutions emanating
% from them, but we are not pursuing it here). However, this branch emanates from the
% trivial solution branch, at $\lambda = (\pi/10)^2$; it is unstable at onset, and it
% does not re-stabilise at the saddle node bifurcation (another unstable eigenvalue
% is added at the saddle-node). In fact, the branch remains unstable for all $\lambda
% \leq 0.7$.

bd = bd1; sol = sol1;
nSol = size(bd,1);
eVals = zeros(nSol,5);
for i = 1:size(bd,1)
  id = i;
  [~,D] = ComputeStability(bd,sol,p0,Dxx,id);
  [d,ix] = sort(diag(D),'descend'); 
  eVals(i,:) = d(1:5);
end

% Plotting branch
figure;
subplot(2,2,[1 2]);
plot(bd(:,2),bd(:,3)); 
% hold on; ii = 41;
% plot(bd(ii,2),bd(ii,3),'r*'); 
xlabel('\lambda'); ylabel('S_0');
title('Branch of sigmoidal solutions');

% Plotting eigenvalues
subplot(2,2,[3 4]);
plot(1+bd(:,1),real(eVals),'k-', 1+bd(:,1), bd(:,2),'-'); grid on; 
title('Most unstable eigenvalues along the branch');
xlabel('s'); ylabel('Re \mu_i');
lgd = legend({'\mu_1','\mu_2','\mu_3','\mu_4','\mu_5','\lambda'});
lgd.Location = 'southeast';
% hold on;
% plot(ii, eVals(ii,3),'r*');
% plot(ii, bd(ii,2),'r*');
% hold off;

%% 
% Unlike the other branches, the branch of homogeneous solution is initially
% unstable, but it has stable segments. There are several bifurcation in addition to
% the saddle-node bifurcation.
bd = bd2; sol = sol2;
nSol = size(bd,1);
eVals = zeros(nSol,5);
for i = 1:size(bd,1)
  id = i;
  [~,D] = ComputeStability(bd,sol,p0,Dxx,id);
  [d,ix] = sort(diag(D),'descend'); 
  eVals(i,:) = d(1:5);
end

% Plotting branch
figure;
subplot(2,2,[1 2]);
plot(bd(:,2),bd(:,3)); 
% hold on; ii = 41;
% plot(bd(ii,2),bd(ii,3),'r*'); 
xlabel('\lambda'); ylabel('S_0');
title('Branch of homogeneous solutions');

% Plotting eigenvalues
subplot(2,2,[3 4]);
plot(1+bd(:,1),real(eVals),'k-', 1+bd(:,1), bd(:,2),'-'); grid on; 
title('Most unstable eigenvalues along the branch');
xlabel('s'); ylabel('Re \mu_i');
lgd = legend({'\mu_1','\mu_2','\mu_3','\mu_4','\mu_5','\lambda'});
lgd.Location = 'southeast';

%% 
% One consideration is that, when we computed the brute-force bifurcation diagram in
% Figure 1.2 of the tutorial, using a time stepper, the sigmoidal and bump steady
% states both "seemed" stable at $\lambda = 0.7$. On the other hand, the stability
% analysis above tells us taht both are unstable. We can see from the results below
% that the solution spends a long time close to the unstable bump equilibrium, before
% being ejected from it, and being attracted by the homogeneous stable state.
% First let us make a simulation for $t \in [0,50]$: the bump solution would appear
% as a stable. We know this is in contraddiction with the numerical bifurcation
% analysis.

p = [1; 0.7; 0; 1; 1];
u0 = 0.4*cos(2*pi/10*x);
rhs = @(t,u) AllenCahn(u,p,Dxx);
tSpan = [0 50];
[t,UHist] = ode15s(rhs,tSpan,u0);
[X,T] = meshgrid(x,t);
figure; surf(X,T,UHist); shading interp;
xlabel('x'); ylabel('t'); zlabel('u');

%% 
% Let's repeat the simulation for $t \in [0,6000]$: after about 1500 time units the
% bump reveals itself as unstable, and the solution reaches the homogeneous
% stable steady state. This is in accordance with the numerical bifurcation analysis
tSpan = [0 6000];
[t,UHist] = ode15s(rhs,tSpan,u0);
[X,T] = meshgrid(x,t);
figure; surf(X,T,UHist); shading interp;
xlabel('x'); ylabel('t'); zlabel('u');

