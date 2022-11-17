function [x,Dx,Dxx] = PeriodicDiffMat(xSpan,nx)

      % Gridpoints
      a = xSpan(1); b = xSpan(2);
      hx = (b-a)/nx;
      x = a+[0:nx-1]'*hx;

      % Auxiliary vecor
      e = ones(nx,1);

      % First order differentiation matrix
      Dx = spdiags([-e e],[-1 1],nx,nx);
      Dx(1,nx) =  -1; Dx(nx,1) = 1;
      Dx = Dx/(2*hx);
 
      % Second order differentiation matrix
      Dxx = spdiags([e -2*e e],-1:1,nx,nx);
      Dxx(1,nx) = 1; Dxx(nx,1) = 1; 
      Dxx = Dxx/(hx^2);

end
