function plotHandle = ExploreBifDiag(bd,sol,id,D,V,x)

  % Bifurcation digaram
  plotHandle = figure;
  plotHandle.Position = [2280 313 1144 1024];
  subplot(5,5,[1:3]); 
  plot(bd(:,2),bd(:,3),'-'); hold on; plot(bd(id,2),bd(id,3),'*');
  xlabel('lambda'); ylabel('2-norm');

  % Solution profile
  u = sol(id,:)'; p(2) = bd(id,2);
  subplot(5,5,[6:8 11:13]); plot(x,u,'.-');
  xlabel('x');  ylabel('u_*'); title(['Steady state, lambda=' num2str(p(2))]);

  % Spectrum
  [d,ix] = sort(diag(D),'descend'); 
  subplot(5,5,[16:18 21:23]); plot(real(d),imag(d),'*'); hold on; xline(0,'red');
  xlabel('real part');  ylabel('imaginary part'); title('Eigenvalues');
  grid on;

  % Top 5 eigenmodes
  for j = 1:5
    if j == 1
      subplot(5,5,[4 5]);
    else
      subplot(5,5,5*(j-1)+3+[1 2]);
    end
    plot(x,real(V(:,ix(j))),'DisplayName','Re \psi(x)'); hold on;
    plot(x,imag(V(:,ix(j))),'DisplayName','Im \psi(x)'); hold off;
    title(['lambda = ' num2str(real(d(j))) ' + i ' num2str(imag(d(j))) ]);
    xlabel('x'); ylim([-0.3 0.3]);
    lgd = legend;
    lgd.Location = 'eastoutside';
  end

end
