function [ ib id iu ] = find_closest( a, b )
% [ ib id iu ] = find_closest( a, b )
% Given vectors a and b, this function returns, for each element of 
% a, the index of the closest element of b (ib), as well as the closest 
% lesser element (id) and greater element (iu). It also  The algorithm used
% has complexity (m+n)log(m+n); code adapted from Roger Stafford:
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/243878 

 m = size(a,2); n = size(b,2);
 [c,p] = sort([a,b]);
 q = 1:m+n; q(p) = q;
 t = cumsum(p>m);
 r = 1:n; r(t(q(m+1:m+n))) = r;
 s = t(q(1:m));
 id = r(max(s,1));
 iu = r(min(s+1,n));
 [d,it] = min([abs(a-b(id));abs(b(iu)-a)]);
 ib = id+(it-1).*(iu-id);

end