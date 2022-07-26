function [ F t ] = cat_rate_curves( st,t_event,window,g_width,bin )
% [ F t ] = cat_rate_curves( st,t_event,window,g_width,bin )
% Given a cell array st of spike times, this function extracts the firing
% rates within [window(1) window(2)] of the event times t_event. The spike
% trains are smoothed with a gaussian kernel of width g_width, and a bin
% size of bin. If t_event is a vector, the same event times are used for
% all spike trains. If it is a cell array with the same length as st,
% different event times are used for each train.

t = window(1):bin:window(2);

if ~iscell(t_event);
    foo = t_event;
    clear t_event
    for i = 1:length(st);
        t_event{i} = foo;
    end
end

n_samp = round(diff(window)/bin)+1;

F = cell(1,length(st));
%for i = 1:length(st);
%    [nbr i_nbr] = nndt(t_event{i}-window(1),t_rate);
%    [foo min_ind] = min(abs(nbr),[],1);
%    min_ind = sort([i_nbr(1,min_ind == 1) i_nbr(2,min_ind == 2)]);
%    ind = [];
%    for j = 1:length(t_event{i});
%        ind = [ind min_ind(j):(min_ind(j)+n_samp-1)];
%    end
%    F{i} = reshape(rate(i,ind),n_samp,length(t_event{i}))';
%end

for i = 1:length(st)
    %t_samp = [];
    t_samp = zeros(1,n_samp*length(t_event{i}));
    [rate t_rate] = firing_rate(st(i),g_width,bin,[],'gaussian');
    F{i} = zeros(length(t_event{i}),n_samp);
    for j = 1:length(t_event{i})
        %t_samp = [t_samp (t_event{i}(j) + (window(1):bin:window(2)))];
        t_samp((1+(j-1)*n_samp):(j*n_samp)) = t_event{i}(j) + (window(1):bin:window(2));
    end
    F{i} = reshape(interp1(t_rate,rate,t_samp),n_samp,length(t_event{i}))';
    %display(i);
end

end