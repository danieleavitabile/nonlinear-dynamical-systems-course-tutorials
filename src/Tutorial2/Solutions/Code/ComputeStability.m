function [V,D] = ComputeStability(bd,sol,p0,Dxx,id)

  % Branch and solution identifier
  display('Linear stability of solution') 
  fprintf('ID = %4d,    p(2) = %5.4e,    ||u|| = %5.4e\n',id, bd(id,2),bd(id,3));

  % Compute Jacobian at the solution
  p0(2) = bd(id,2); u = sol(id,:)';
  [~,J] = AllenCahn(u,p0,Dxx);

  % Compute eigenvalues and eigenvectors
  [V,D] = eig(full(J));

end
