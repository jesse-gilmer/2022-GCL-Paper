clear all
close all
rng('shuffle','philox')
for T = [25]
    clearvars -except T
    mRes = struct();

nMF = 100;
nGC = 6000;
conv = 4;

mu = 1;
sigma = .33;

options.mu = mu;
options.sigma = sigma;
Li  = 0;
for L = [.001]
    Li = Li + 1;
    figure(Li)
    title(L)
    hold on
options.lrnRate = L;
options.nTrials = 1000;

t = imread('cat_g.png');
t = t(:,:,1);
tel = numel(t);

nPts = tel;

ntr = 1;
i = 1;
    
    
%     for A = -4:.25:2.5
for mtr = 1
    AA = 0;
% for A = linspace(-1,1,6)
for A = [-2]
    AA = AA + 1;
    for q = 1:ntr
%       try  
%% Project params:
alpha = A;
mRes(i).options = options;
options.tau = T;
mRes(i).alpha = alpha;
mRes(i).conv = conv;
mRes(i).nGC = nGC;
mRes(i).nMF = nMF;


%% Make input signal:
input = MakeSignal(nPts,nMF,1,options);
nPts = size(input,2);
mRes(i).nPts = nPts;

%% Make target signal:

options2 = options;
options2.tau = 100;
target = MakeSignal(nPts,1,1,options2);
% target = rand(size(target));
tt = 1:nPts;
target = cos(10*2*pi*(tt/1000));
target = (target-min(target))./(max(target)-min(target));

target = double(t(:)');


%% Make GCs:
W1 = zeros(nGC, nMF);
for k = 1:nGC
    sel = randsample(nMF, conv);
    W1(k,sel) = 1/conv;
end
mu = mean(mean(input));
sigma = mean(std(input,[],2));
thr = mu + alpha * sigma;

gc = W1 * input;   
gc = gc - thr;
gc(gc < 0) = 0;

CC = corr(gc');
diagC = tril(ones(size(CC)),-1);
mRes(i).gcCorr = nanmean(CC(diagC == 1));

target = target./max(target);
[results] = RunFilter(input,target,gc,options);
results.mse(end)
cat_copy = results.W * gc;
cat_copy = reshape(cat_copy, size(t));
figure(AA)
subplot(2,2,1)
hold on
imshow(cat_copy)
title(results.mse(end))
subplot(2,2,2)
hold on
plot(target)
plot(results.W * gc)
subplot(2,2,3)
plot(results.mse(2:end))

figure(99999)
subplot(3,3,AA)
imshow(cat_copy)
title("Threshold: " + A + ", MSE: " + results.mse(end))

[results] = RunFilter(input,target,input,options);
results.mse(end)
cat_copy = results.W * input;
cat_copy = reshape(cat_copy, size(t));
figure(AA+100)
subplot(2,2,1)
hold on
imshow(cat_copy)
title(results.mse(end))
subplot(2,2,2)
hold on
plot(target)
plot(results.W * input)
subplot(2,2,3)
plot(results.mse(2:end))

figure(999992)
subplot(3,3,AA)
imshow(cat_copy)
title("Threshold: " + A + ", MSE: " + results.mse(end))

    end
    drawnow
end
end
end
end
