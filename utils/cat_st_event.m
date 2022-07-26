function [ st_event ind ] = cat_st_event( st,t_event,window )
% [ st_event ind ] = cat_st_event( st,t_event,window )
% Given a cell array st of spike times, a vector of event times, and upper
% and lower limits for a peri-event window, this function returns a cell
% array of spike times centered on each event.

for i = 1:length(st)
    for j = 1:length(t_event)
        ind{i,j} = st{i}>t_event(j)+window(1) & ...
            st{i}<t_event(j)+window(2);
        st_event{i,j} = st{i}(ind{i,j})-t_event(j);
    end
end
        

end

