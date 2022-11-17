function [x,Dx,Dxx] = NeumannDiffMat(xSpan,nx)

      % Gridpoints
      a = xSpan(1); b = xSpan(2);
      x = linspace(a,b,nx)';
      hx = (b-a)/(nx-1);

      % Auxiliary vecor
      e = ones(nx,1);

      % First order differentiation matrix
      Dx = spdiags([-e e],[-1 1],nx,nx);
      Dx(1,:) = 0; Dx(nx,:) = 0;
      Dx = Dx/(2*hx);
 
      % Second order differentiation matrix
      Dxx = spdiags([e -2*e e],-1:1,nx,nx);
      Dxx(1,2) = 2; Dxx(nx,nx-1) = 2; 
      Dxx = Dxx/(hx^2);

end
