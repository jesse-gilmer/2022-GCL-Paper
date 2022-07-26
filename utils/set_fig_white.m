% Set current figure white for printing.

fig_hand = gcf;
fig_hand.Color = [1 1 1];
fig_hand.InvertHardcopy = 'off';
set(fig_hand,'PaperPositionMode','auto')
set(fig_hand,'renderer','painter');