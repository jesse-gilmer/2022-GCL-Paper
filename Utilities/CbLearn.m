function [mse, net, mf, target] = CbLearn(mf, gc, target, train, noise, dt)

% generate mossy fiber rate functions
if isstruct(mf)
    nPts = round(mf.T/dt);
    nMF = mf.N;
    mu = mf.mean;
    sigma = mf.sd;
    tau = mf.tau;
    r0 = randn(nMF,1)*sigma;
    mf = zeros(nMF, nPts) + r0;
    for i = 2:nPts
        mf(:,i) = mf(:,i-1)*exp(-dt/tau) + sigma*sqrt(1-exp(-2*dt/tau))*randn(nMF,1);
    end
    mf = mf + mu; 
else
    nMF = size(mf,1);
    nPts = size(mf,2);
end

% generate target rate functions
if isstruct(target)
    nPk = target.N;
    mu = target.mean;
    sigma = target.sd;
    tau = target.tau;
    r0 = randn(nPk,1)*sigma;
    target = zeros(nPk, nPts) + r0;
    for i = 2:nPts
        target(:,i) = target(:,i-1)*exp(-dt/tau) + sigma*sqrt(1-exp(-2*dt/tau))*randn(nPk,1);
    end
    target = target + mu;
else
    nPk = size(target,1);
end

% generate connection matrix to GCs
nGC = gc.N;
gcConv = gc.conv;
if nGC > 0
    W1 = zeros(nGC, nMF);
    for k = 1:nGC
        sel = randsample(nMF, gcConv);
        W1(k,sel) = 1;
    end
else
    nGC = nMF;
    W1 = eye(nMF);
    gc.thr = nan(1);
    gc.alphaZ = nan(1);
end
net.W1 = W1;

% set GC threshold
if ~isnan(gc.alphaZ)
    mu = gcConv * mean(mean(mf));
    sigma = sqrt(gcConv) * mean(std(mf,[],2));
    thr = mu + gc.alphaZ * sigma;
else
    thr = gc.thr;
end

% unpack remaining parameters
nTr = train.nTr;
lrnRate = train.rate;
batchLearn = train.batch;
normW = train.normalize;
sigma = noise.sd;
tau = noise.tau;

% compute granule cell response
gc = W1 * mf;   % pre-threshold activity
if ~isnan(thr)
    gc(gc < thr) = 0;
    gc(gc > 0) = gc(gc > 0) - thr;  % thresholded activity
end

% generate connection matrix to Purkinjie cells
W2 = rand(nPk, nGC); 
if normW
    % scale weight matrix so that mean response = mean target
    response = W2 * gc;
    W2 = W2 .* mean(target,2) ./ mean(response,2);
end
net.W2_init = W2;

% training loop
mse = zeros(nPk, nTr+1);
for k = 1:nTr
    
    % compute non-plastic response
    pk = W2 * gc;
    mse(:,k) = mean((pk-target).^2, 2);
    
    % initialize noise
     noise = sigma * randn(nPk, 1);
    
    % implement learning
    for i = 1:nPts
        noise = noise*exp(-dt/tau) + sigma*sqrt(1-exp(-2*dt/tau))*randn(nPk,1);
        if batchLearn
            error = pk(:,i) + noise - target(:,i);
        else
            error = W2 * gc(:,i) + noise - target(:,i);
        end
        W2 = W2 - error .* gc(:,i)' * lrnRate;
    end
    
end

% final performance
pk = W2 * gc;
mse(:, nTr+1) = mean((pk-target).^2, 2);

net.W2_final = W2;
net.thr = thr;