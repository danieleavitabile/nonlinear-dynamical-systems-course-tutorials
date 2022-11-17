function plotHandle = PlotHistory(x,t,U,p,parentHandle)

  %% Extract number of components
  numComp = 2;
  nx = size(U,2)/2;

  %% Assign solution label
  solLabel(1).name = "U";
  solLabel(2).name = "V";

   %% Position and eventually grab figure
   if isempty(parentHandle)
     %scrsz = get(0,'ScreenSize');
     % plotHandle = figure('Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/4 scrsz(4)/4]);
     plotHandle = figure();
     parentHandle = plotHandle;
    else
      plotHandle = parentHandle;
   end
   figure(parentHandle);

   %% Grid
   [T,X] = meshgrid(t,x);

   %% Plots
   for k = 1:numComp
     subplot(1,numComp,k)
     % pcolor(X,T,U(:,idx(:,k))'); shading interp; view([0 90]);
     surf(X,T,U(:,nx*(k-1)+[1:nx])'); shading interp; view([0 90]);
     title(solLabel(k).name);
     xlabel('x'); ylabel('t');
   end

   %% Save
   % print -dtiff history.tiff

end

