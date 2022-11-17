function plotHandle = PlotSteadyState(x,u,parentHandle,lims,aTitle,replace)

   % Eventually grab an older figure (parentHAndle)
   if isempty(parentHandle)
     plotHandle = figure();
    else
      plotHandle = parentHandle;
   end
   figure(plotHandle);

   % Determine if we want to replace the figure or write on it
   if isempty(replace)
     replace = true;
   end
   if ~replace
     hold on;
   end

   % Plot
   plot(x,u);
   xlabel('x'); ylabel('u_*(x)');

   % Optionally add title and axes
   if ~isempty(aTitle);
     title(aTitle);
   end

   % Optionally change limits
   if ~isempty(lims);
     axis(lims);
   end


   if replace

end

