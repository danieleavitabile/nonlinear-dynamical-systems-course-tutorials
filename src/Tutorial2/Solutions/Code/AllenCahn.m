function [F,DFDU] = AllenCahn(u,p,Dxx)

  % Rename parameters
  nu     = p(1); 
  lambda = p(2); 
  alpha  = p(3);
  beta   = p(4);
  gamma  = p(5);
  nx     = length(u);

  % Right-hand side
  F = nu*Dxx*u + lambda*u + alpha*u.^2 + beta*u.^3 - gamma*u.^5;

  if nargout > 1
    v = (lambda + 2*alpha*u + 3*beta*u.^2 - 5*gamma*u.^4);
    DFDU = nu*Dxx + spdiags(v,0,nx,nx);
  end

end
