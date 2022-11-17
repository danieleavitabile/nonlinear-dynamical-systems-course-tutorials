% Clean
clear all, close all, clc

% Setup
u0 = 3;                     % Initial conditions
p0 = [2,1];                 % Parameters [lambda,beta]
prob = @(u,p) Poly135(u,p); % Problem definition

% Continuation parameter index: 
% icp = 1 means continuation in p(1)
icp = 1;     

% Parameter stepsize: 
% a negative stepsize gives an initially decreasing parameter
ds = -0.05; 

% Number of continuation steps
nSteps = 80;    

% Options to pass to fsolve
opts = optimset('Display','off', 'TolFun',1e-10,...
                'MaxIter',50, 'Jacobian','on');

% Launch continuation
[bd,sol] = SecantContinuation(prob,u0,p0,icp,ds,nSteps,opts);

% Plotting
plot(bd(:,2), bd(:,3), '.-');
xlabel('lambda'); ylabel('2-norm');

% Right-hand side function
function [F,DFDU] = Poly135(u,p);
   F = p(1)*u + p(2)*u^3 - u^5; 
   if nargout > 1
     DFDU = p(1) + 3*p(2)*u^2 -5*u^4; 
   end
end

