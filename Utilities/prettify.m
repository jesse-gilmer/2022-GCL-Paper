function [] = prettify(h)
h1 = get(h,'children');

%% Make the axes figure normal (bottom and left only)
set(h1,'box','off')

%% Change the figure qualities
set(h,'color','w');
% set(h1,'fontsize',16);
set(gca,'TickDir','out');
set(gcf, 'Position', [100 100 840 740]);
end