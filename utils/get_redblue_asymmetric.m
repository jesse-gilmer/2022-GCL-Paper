function [ col ] = get_redblue_asymmetric( nr,nb )
% [ col ] = get_redblue_asymmetric( nr,nb )

nn = 1001;
nn_half = .5+nn/2;

foo = get_redblue(nn);
Mr = foo(1:nn_half,:);
Mb = foo(nn_half:end,:);

for k = 1:nr
    colr(k,:) = Mr(round(1+(k-1)*(nn_half-1)/(nr-1)),:);
end
for k = 1:nb
    colb(k,:) = Mb(round(1+(k-1)*(nn_half-1)/(nb-1)),:);
end

col = [colr; colb];


end