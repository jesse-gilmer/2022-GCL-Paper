%% Make Figure 3:
clear all
close all

fdir = "Figure 3 Updated"
mkdir(fdir)

%% 3A:

nMF = 50;
nGC = 3000;
conv = 4;
options.mu = .5;
options.sigma = .1;
options.lrnRate = 1e-2;
options.tau = 100;
nPts = 200;

options2 = options;
options2.tau = 200;
target = MakeSignal(nPts,1,1,options2);
target = (target-min(target))./(max(target)-min(target));

input = MakeSignal(nPts,nMF,1,options);

W1 = zeros(nGC, nMF);
for k = 1:nGC
    sel = randsample(nMF, conv);
    W1(k,sel) = 1/conv;
end
mu = mean(mean(input));
sigma = mean(std(input,[],2));
thr = mu + 0 * sigma;

gc = W1 * input;
gc = gc - thr;
gc(gc < 0) = 0;

h = figure();
ii = 1;
for n = [1 5 50]
    options.nTrials = n;
    subplot(3,1,ii)
    ii = ii + 1;

    [results] = RunFilter(input,target,gc,options);
    plot(target,'k')
    hold on
    plot(sum(gc.*results.W'),'r')
    
end

prettify(gcf);
saveas(h,fdir + "/Figure_3A.emf")



