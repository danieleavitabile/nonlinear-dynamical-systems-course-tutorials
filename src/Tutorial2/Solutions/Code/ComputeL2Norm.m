function S = ComputeL2Norm(u,x)

  % Rename parameters
  nx = length(x); hx = x(2) - x(1); Lx = x(end)-x(1);

  % Integration weights
  w = ones(size(x)); w([1 nx]) = 0.5;  w = w*hx;

  % Compute square of the integral of |u(x)|^2,
  % And scale it by Lx
  S = sqrt(w'*u.^2)/sqrt(Lx);

end
