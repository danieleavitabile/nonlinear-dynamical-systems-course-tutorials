function [F,DFDZ] = Schnakenberg(z,p,Dxx)

  % Rename parameters
  b  = p(1); 
  d  = p(2); 

  % Ancillary variables and solution split
  nx = length(z)/2; iU = 1:nx; iV = nx+iU;
  u = z(iU); v = z(iV);

  % Function handles for reaction terms, and their derivatives
  f = @(u,v) -u+v.*u.^2; dfdu = @(u,v) -1 +2*u.*v; dfdv = @(u,v)  u.^2;
  g = @(u,v)  b-v.*u.^2; dgdu = @(u,v)    -2*u.*v; dgdv = @(u,v) -u.^2;

  % Right-hand side
  F = zeros(size(z));
  F(iU) =   Dxx*u + f(u,v);
  F(iV) = d*Dxx*v + g(u,v);

  if nargout > 1
    DFDZ = spdiags([],[],2*nx,2*nx);
    DFDZ(iU,iU) =   Dxx + spdiags(dfdu(u,v),0,nx,nx);
    DFDZ(iU,iV) =         spdiags(dfdv(u,v),0,nx,nx);
    DFDZ(iV,iU) =         spdiags(dgdu(u,v),0,nx,nx);
    DFDZ(iV,iV) = d*Dxx + spdiags(dgdv(u,v),0,nx,nx);
  end

end
