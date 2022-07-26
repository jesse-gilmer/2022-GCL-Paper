function [ q h ] = fdr_ben_hoch( p,q_thresh )
%function [ q h ] = fdr_ben_hoch( p,q_thresh )

[p_sorted order] = sort(p,'ascend');
q = (length(p)./(1:length(p))).*p_sorted;

[foo resort] = sort(order);
q = q(resort);
h = zeros(size(p));
h(q < q_thresh) = 1;

q(q>1) = 1;

end