function F = AllenCahn(u,p,Dxx)

  % Rename parameters
  nu     = p(1); 
  lambda = p(2); 
  alpha  = p(3);
  beta   = p(4);
  gamma  = p(5);

  F = nu*Dxx*u + lambda*u + alpha*u.^2 + beta*u.^3 - gamma*u.^5;

end
