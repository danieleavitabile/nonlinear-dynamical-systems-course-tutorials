fig = gcf;
figName = 'bd';
%ylim([1e-6 1e2]);
%xlim([1 2e3]);
%legend('NumColumns',2);
set(0,'defaulttextinterpreter','latex')
%xlabel(''); ylabel('$\displaystyle{\max_k \| u(t_k) - u_n(t_k) \|_\infty}$');

% define figure properties
opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = '~/Desktop';
opts.width      = 8;
opts.height     = 6;
opts.fontType   = 'Times';
opts.fontSize   = 9;
opts.gridLineStyle = '-'; 
opts.minorGridLineStyle = '-';
opts.gridAlpha = 0.05; 
opts.minorGridAlpha = 0.05;

fig.Units               = 'centimeters';
fig.Position(3)         = opts.width;
fig.Position(4)         = opts.height;

% set text properties
set(fig.Children, ...
    'FontName',     'CMU Serif', ...
    'FontSize',     10);

% remove unnecessary white space
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02),'GridLineStyle',opts.gridLineStyle)
grid on; set(gca,'GridLineStyle',opts.gridLineStyle,...
                 'MinorGridLineStyle',opts.minorGridLineStyle,...
		 'GridAlpha',opts.gridAlpha,...
		 'MinorGridAlpha',opts.minorGridAlpha);

% export to png
fig.PaperPositionMode   = 'auto';
print([opts.saveFolder '/' figName '.png'], '-dpng', '-r600')
saveas(gcf,[opts.saveFolder '/' figName '.pdf'],'pdf')
