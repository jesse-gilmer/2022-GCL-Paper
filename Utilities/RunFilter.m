function [output] = RunFilter(input,target,gc,options)

nPk = size(target,1);
nPts = size(target,2);
nMF = size(input,1);
nGC = size(gc,1);

nTr = options.nTrials;
lrnRate = options.lrnRate;

% generate connection matrix to Purkinjie cells
W2 = rand(nPk, nGC);
% scale weight matrix so that mean response = mean target
response = W2 * gc;
W2 = W2 .* mean(target,2) ./ mean(response,2);
% W2 = W2 ./ nGC;
net.W2_init = W2;

%%
% training loop
mse = zeros(nPk, nTr+1);
for k = 1:nTr
    
    % compute non-plastic response
    pk = W2 * gc;
    mse(:,k) = mean((pk-target).^2, 2);
    
    % implement learning
    for i = 1:nPts
        error = W2 * gc(:,i) - target(:,i);
        W2 = W2 - (error .* gc(:,i) .* lrnRate)';
    end

    W2(isnan(W2)) = 0;
end

% final performance
pk = W2 * gc;
mse(:, nTr+1) = mean((pk-target).^2, 2);

output.encode = pk;
output.mse = mse;
output.W = W2;

end

