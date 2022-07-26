function [ rates t_rates ] = firing_rate( st,window,bin,tlim,smooth )
% [ rates t_rates ] = firing_rate( st,window,bin,tlim,smooth )
% This function returns an array of firing rates. Time limits
% should be in gtms. Smoothing options are 'rect', 'gaussian', and
% 'causal', where the window value specifies the length of the boxcar
% function and the width of the gaussian kernel, respectively (both values
% in ms). The option 'causal' uses a truncated half-gaussian kernel that
% only uses spikes in the past, and 'acausal' uses only the future.

% Set default time limits.
if isempty(tlim)
    minlim = inf;
    maxlim = 0;
    for i = 1:length(st)
        if ~isempty(st{i})
        minlim = min(minlim,min(st{i}));
        maxlim = max(maxlim,max(st{i}));
        end
    end
    tlim = [minlim maxlim];
end

% Iterate over cells.
for i = 1:length(st)
    if sum(st{i}>tlim(1) & st{i}<tlim(2)) > 0
        train = histc(st{i},tlim(1):bin:tlim(2));
    else
        train = zeros(1,length(tlim(1):bin:tlim(2)));
    end
    if strcmp(smooth,'rect') == 1
        kernel = (1000/window)*ones(1,round(window/bin));
    elseif strcmp(smooth,'gaussian') == 1
        kernel = 1000*(1/(window*sqrt(2*pi)))*exp((-1/(2*(window^2)))*((-(window*5):bin:(window*5)).^2));
	elseif strcmp(smooth,'causal') == 1
        kernel = 1000*(1/(window*sqrt(2*pi)))*exp((-1/(2*(window^2)))*((-(window*5):bin:(window*5)).^2));
        kernel(1:floor(length(kernel)/2)) = 0;
        kernel = 2*kernel;
	elseif strcmp(smooth,'acausal') == 1
        kernel = 1000*(1/(window*sqrt(2*pi)))*exp((-1/(2*(window^2)))*((-(window*5):bin:(window*5)).^2));
        kernel(ceil(length(kernel)/2):end) = 0;
        kernel = 2*kernel;
    end
    rates(i,:) = conv(train,kernel,'same');
end

t_rates = tlim(1) + bin*(0:(size(rates,2)-1));

end