function [z nn] = nndt(x,y)

%NNDT   Spike nearest neighbor time differences.
%
%   Z = NNDT(POST, PRE) given spike times POST and PRE returns the 
%   time differences POST-PRE between each spike in POST and the 
%   spikes in PRE immediately preceding and following.
%  
%   [Z I] = NNDT(POST, PRE) also returns the indices I of the
%   spikes in PRE used in the time difference computation.

x = x(:); y = y(:)';
edge = [-inf y inf];
[foo bin] = histc(x, edge);
y(end+1) = NaN; 
nn = [bin-1 bin];
nn(nn<1) = length(y);
nn(nn>length(y)) = length(y);
z = [x x] - y(nn); 
z = z'; 
nn = nn';
% z = z(:); 
% z(isnan(z)) = [];
