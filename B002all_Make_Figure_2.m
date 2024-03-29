%% HEADER

% B002all_Make_Figure_2.m
% From: An emergent temporal basis set robustly supports cerebellar time-series learning
% Authors: Jesse I. Gilmer, Michael A. Farries, Zachary Kilpatrick, Ioannis Delis, Abigail L. Person
% 2022, BiorXiv, https://www.biorxiv.org/content/10.1101/2022.01.06.475265v1
%
% This script was written by JG.
% Data used in this script comes from https://zenodo.org/record/5140528#.YuA103bMJD9
%
% INPUTS: N/A.
%
% OUTPUTS: 6 figure outputs to /Figure 2/


%% Make Figure 2:
close all;
clear all;

sdir = 'Figure 2/';
mkdir(sdir)

nMF = 50; %Default 50
nGC = 300; %Default 3000
nPk = 1;
conv = 4;

mu = 4;
sigma = 1;

options.mu = mu;
options.sigma = sigma;
options.lrnRate = 1e-5;
options.nTrials = 1000; %Default 1000

nPts = 1000; %Default 1000
i = 1;

% Project params:
alpha = 0;
options.tau = 100;

input = MakeSignal(nPts,nMF,1,options);
nPts = size(input,2);

%% Make baseline GC activity:

nPts = nPts;
nMF = size(input,1);

% generate connection matrix to Purkinjie cells
W2 = rand(nPk, nGC);
% scale weight matrix so that mean response = mean target
W1 = zeros(1,nMF);
for k = 1:nGC
    sel = randsample(nMF, 4);
    W1(k,sel) = 1/4;
end
mu = mean(mean(input));
sigma = mean(std(input,[],2));

alpha = 0;
thr = mu + alpha * sigma;

gc = W1 * input;
gc = gc - thr;
gc(gc < 0) = 0;

[sgc,xgc] = max(gc');
[meow,mix] = sort(xgc);
outmix  = mix;

figure();
imagesc(gc(mix,:))
prettify(gcf)
xlabel('time (ms)')
ylabel('GC unit')
colorbar(); colormap turbo
caxis([0 2])
print(gcf, '-dsvg', sdir + "/Figure_2_B_Noise0.svg")


figure()
plot(input(1,:))

prettify(gcf)
print(gcf, '-dsvg', sdir + "/Figure_2_A_Noise0example.svg")
%% Make noise:

mll = [1:100];
ii = 0;
sets = 3;
noiseAmp = [.33]
noise = randn(size(input));
noise = noise'.*var(input');
noise = noise'.*noiseAmp;
input2 = input + noise;

mu = mean(mean(input2));
sigma = mean(std(input2,[],2));
thr2 = mu + alpha * sigma;

gc2 = W1 * input2;
gc2 = gc2 - thr;
gc2(gc2 < 0) = 0;


[sgc,xgc2] = max(gc2');
[meow2,mix2] = sort(xgc2);

figure();
imagesc(gc2(outmix,:))
prettify(gcf)
xlabel('time (ms)')
ylabel('GC unit')
colorbar(); colormap turbo
% 

xlabel('time (ms)')
caxis([0 2])
disp('Painters') %https://www.mathworks.com/matlabcentral/answers/92521-why-does-matlab-not-export-eps-files-properly
set(gcf,'renderer','Painters')
print(gcf, '-dsvg', sdir + "/Figure_2_B_Noisep25.svg")

figure()
subplot(3,1,1)
plot(input(1,:))
subplot(3,1,2)
plot(noise(1,:))
subplot(3,1,3)
plot(input2(1,:))

prettify(gcf)
print(gcf, '-dsvg', sdir + "/Figure_2_A_Noisep25examples.svg")

%%
mll = [1:100];
ii = 0;
sets = 3;
noiseAmp = [1]
noise = randn(size(input));
noise = noise'.*var(input');
noise = noise'.*noiseAmp;
input2 = input + noise;

mu = mean(mean(input2));
sigma = mean(std(input2,[],2));
thr2 = mu + alpha * sigma;

gc2 = W1 * input2;
gc2 = gc2 - thr;
gc2(gc2 < 0) = 0;


[sgc,xgc2] = max(gc2');
[meow2,mix2] = sort(xgc2);

figure();
imagesc(gc2(outmix,:))
prettify(gcf)
xlabel('time (ms)')
ylabel('GC unit')
colorbar(); colormap turbo

caxis([0 2])

xlabel('time (ms)')


disp('Painters') %https://www.mathworks.com/matlabcentral/answers/92521-why-does-matlab-not-export-eps-files-properly
set(gcf,'renderer','Painters')
print(gcf, '-dsvg', sdir + "/Figure_2_B_Noisep50.svg")

figure()
subplot(3,1,1)
plot(input(1,:))
subplot(3,1,2)
plot(noise(1,:))
subplot(3,1,3)
plot(input2(1,:))

prettify(gcf)
print(gcf, '-dsvg', sdir + "/Figure_2_A_Noisep50examples.svg")
%%
mholder = [];
miholder = [];
nMF = 50; %Default 50
nGC = 3000; %Default 3000
nPk = 1;

for i = 1:100
    
    input = MakeSignal(nPts,nMF,1,options);
    % generate connection matrix to Purkinjie cells
    W2 = rand(nPk, nGC);
    % scale weight matrix so that mean response = mean target
    W1 = zeros(1,nMF);
    for k = 1:nGC
        sel = randsample(nMF, 4);
        W1(k,sel) = 1/4;
    end
    mu = mean(mean(input));
    sigma = mean(std(input,[],2));
    
    alpha = 0;
    thr = mu + alpha * sigma;
    
    gc = W1 * input;
    gc = gc - thr;
    gc(gc < 0) = 0;
    
    [sgc,xgc] = max(gc');
    [meow,mix] = sort(xgc);
    
    noiseAmp = [.33];
    noise = randn(size(input));
    noise = noise'.*var(input');
    noise = noise'.*noiseAmp;
    input2 = input + noise;
    
    mu = mean(mean(input2));
    sigma = mean(std(input2,[],2));
    thr2 = mu + alpha * sigma;
    
    gc2 = W1 * input2;
    gc2 = gc2 - thr;
    gc2(gc2 < 0) = 0;
    
    [sgc,xgc2] = max(gc2');
    [meow2,mix2] = sort(xgc2);
    
    mholder = [mholder abs(xgc-xgc2)];
    
    [sgc,xgc] = max(input');
    [sgc,xgc2] = max(input2');
    miholder = [miholder abs(xgc-xgc2)];
end

figure(); hold on; 
histogram(mholder,'binwidth',1,'displaystyle','stairs','edgecolor','r','normalization','cdf')
histogram(miholder,'binwidth',1,'displaystyle','stairs','edgecolor','k','normalization','cdf')
legend(["GC","MF"])
prettify(gcf)
xlabel('Abs. time offset (ms)')
ylabel('Probability')
xlim([0 100])
print(gcf, '-dsvg', sdir + "/Figure_2_C_LowNoise.svg")
%%
mholder = [];
miholder = [];
nMF = 50; %Default 50
nGC = 3000; %Default 3000
nPk = 1;

for i = 1:100
    
    input = MakeSignal(nPts,nMF,1,options);
    % generate connection matrix to Purkinjie cells
    W2 = rand(nPk, nGC);
    % scale weight matrix so that mean response = mean target
    W1 = zeros(1,nMF);
    for k = 1:nGC
        sel = randsample(nMF, 4);
        W1(k,sel) = 1/4;
    end
    mu = mean(mean(input));
    sigma = mean(std(input,[],2));
    
    alpha = 0;
    thr = mu + alpha * sigma;
    
    gc = W1 * input;
    gc = gc - thr;
    gc(gc < 0) = 0;
    
    [sgc,xgc] = max(gc');
    [meow,mix] = sort(xgc);
    
    noiseAmp = [1];
    noise = randn(size(input));
    noise = noise'.*var(input');
    noise = noise'.*noiseAmp;
    input2 = input + noise;
    
    mu = mean(mean(input2));
    sigma = mean(std(input2,[],2));
    thr2 = mu + alpha * sigma;
    
    gc2 = W1 * input2;
    gc2 = gc2 - thr;
    gc2(gc2 < 0) = 0;
    
    [sgc,xgc2] = max(gc2');
    [meow2,mix2] = sort(xgc2);
    
    mholder = [mholder abs(xgc-xgc2)];
    
    [sgc,xgc] = max(input');
    [sgc,xgc2] = max(input2');
    miholder = [miholder abs(xgc-xgc2)];
end

figure(); hold on; hold on;
histogram(mholder,'binwidth',1,'displaystyle','stairs','edgecolor','r','normalization','cdf')
histogram(miholder,'binwidth',1,'displaystyle','stairs','edgecolor','k','normalization','cdf')
legend(["GC","MF"])
prettify(gcf)
xlabel('Abs. time offset (ms)')
ylabel('Probability')
xlim([0 400])
print(gcf, '-dsvg', sdir + "/Figure_2_C_HighNoise.svg")

%%

mlem = struct();
ax = 0;
for a = -2:.25:2
    ax = ax + 1;
    mholder = [];
    miholder = [];
    nMF = 50; %Default 50
    nGC = 3000; %Default 3000
    nPk = 1;

for i = 1:100
    
    input = MakeSignal(nPts,nMF,1,options);
    % generate connection matrix to Purkinjie cells
    W2 = rand(nPk, nGC);
    % scale weight matrix so that mean response = mean target
    W1 = zeros(1,nMF);
    for k = 1:nGC
        sel = randsample(nMF, 4);
        W1(k,sel) = 1/4;
    end
    mu = mean(mean(input));
    sigma = mean(std(input,[],2));
    
    alpha = a;
    thr = mu + alpha * sigma;
    
    gc = W1 * input;
    gc = gc - thr;
    gc(gc < 0) = 0;
    
    [sgc,xgc] = max(gc');
    [meow,mix] = sort(xgc);
    
    noiseAmp = [.33];
    noise = randn(size(input));
    noise = noise'.*var(input');
    noise = noise'.*noiseAmp;
    input2 = input + noise;
    
    mu = mean(mean(input2));
    sigma = mean(std(input2,[],2));
    thr2 = mu + alpha * sigma;
    
    gc2 = W1 * input2;
    gc2 = gc2 - thr;
    gc2(gc2 < 0) = 0;
    
    [sgc,xgc2] = max(gc2');
    [meow2,mix2] = sort(xgc2);
    
    mholder = [mholder abs(xgc-xgc2)];
    
    [sgc,xgc] = max(input');
    [sgc,xgc2] = max(input2');
    miholder = [miholder abs(xgc-xgc2)];
end

mlem(ax).a = a
mlem(ax).mholder = mholder
mlem(ax).miholder = miholder
end

figure()
s = []
p = []
for i = 1:ax
    s = [s mlem(i).mholder];
    p = [p mlem(i).mholder*0+mlem(i).a];
end

boxplot(s,p,'Symbol','')

figure()
s2 = []
p2 = []
for i = 1:ax
    s2 = [s2 mlem(i).miholder];
    p2 = [p2 mlem(i).miholder*0+mlem(i).a];
    [~,sq(i)] = ttest2(mlem(i).mholder,mlem(i).miholder)
end

boxplot(s,p,'Symbol','')

figure()
opt2.scatter = 0;
opt2.SE = 0;
prettyWormPlot(s,p,opt2,[1 0 0])
opt2.scatter = 1;
opt2.SE = 0;
prettyWormPlot(s2,p2,opt2,[0 0 0])
legend(["GC","MF"])
prettify(gcf)
xlabel('Threshold')
ylabel('Offset')
print(gcf, '-dsvg', sdir + "/Figure_2_D_LowNoise2.svg")

figure()
opt2.scatter = 0;
opt2.SE = 1;
prettyWormPlot(s,p,opt2,[1 0 0])
opt2.scatter = 1;
opt2.SE = 1;
prettyWormPlot(s2,p2,opt2,[0 0 0])
legend(["GC","MF"])
prettify(gcf)
xlabel('Threshold')
ylabel('Offset')
print(gcf, '-dsvg', sdir + "/Figure_2_D_LowNoise2SE.svg")

%%

mlem = struct();
ax = 0;
for a = -2:.25:2
    ax = ax + 1;
    mholder = [];
    miholder = [];
    nMF = 50; %Default 50
    nGC = 3000; %Default 3000
    nPk = 1;

for i = 1:100
    
    input = MakeSignal(nPts,nMF,1,options);
    % generate connection matrix to Purkinjie cells
    W2 = rand(nPk, nGC);
    % scale weight matrix so that mean response = mean target
    W1 = zeros(1,nMF);
    for k = 1:nGC
        sel = randsample(nMF, 4);
        W1(k,sel) = 1/4;
    end
    mu = mean(mean(input));
    sigma = mean(std(input,[],2));
    
    alpha = a;
    thr = mu + alpha * sigma;
    
    gc = W1 * input;
    gc = gc - thr;
    gc(gc < 0) = 0;
    
    [sgc,xgc] = max(gc');
    [meow,mix] = sort(xgc);
    
    noiseAmp = [1];
    noise = randn(size(input));
    noise = noise'.*var(input');
    noise = noise'.*noiseAmp;
    input2 = input + noise;
    
    mu = mean(mean(input2));
    sigma = mean(std(input2,[],2));
    thr2 = mu + alpha * sigma;
    
    gc2 = W1 * input2;
    gc2 = gc2 - thr;
    gc2(gc2 < 0) = 0;
    
    [sgc,xgc2] = max(gc2');
    [meow2,mix2] = sort(xgc2);
    
    mholder = [mholder abs(xgc-xgc2)];
    
    [sgc,xgc] = max(input');
    [sgc,xgc2] = max(input2');
    miholder = [miholder abs(xgc-xgc2)];
end

mlem(ax).a = a
mlem(ax).mholder = mholder
mlem(ax).miholder = miholder
end

figure()
s = []
p = []
for i = 1:ax
    s = [s mlem(i).mholder];
    p = [p mlem(i).mholder*0+mlem(i).a];
end

boxplot(s,p,'Symbol','')

figure()
s2 = []
p2 = []
for i = 1:ax
    s2 = [s2 mlem(i).miholder];
    p2 = [p2 mlem(i).miholder*0+mlem(i).a];
    [~,sq(i)] = ttest2(mlem(i).mholder,mlem(i).miholder)
end

boxplot(s,p,'Symbol','')

figure()
opt2.scatter = 0;
opt2.SE = 0;
prettyWormPlot(s,p,opt2,[1 0 0])
opt2.scatter = 1;
opt2.SE = 0;
prettyWormPlot(s2,p2,opt2,[0 0 0])
legend(["GC","MF"])
prettify(gcf)
xlabel('Threshold')
ylabel('Offset')
print(gcf, '-dsvg', sdir + "/Figure_2_D_HiNoise2.svg")

figure()
opt2.scatter = 0;
opt2.SE = 1;
prettyWormPlot(s,p,opt2,[1 0 0])
opt2.scatter = 1;
opt2.SE = 1;
prettyWormPlot(s2,p2,opt2,[0 0 0])
legend(["GC","MF"])
prettify(gcf)
xlabel('Threshold')
ylabel('Offset')
print(gcf, '-dsvg', sdir + "/Figure_2_D_HiNoise2SE.svg")

%%


mlem = struct();
ax = 0;
for n = [0:.1:1 2:1:5]
    a = 1
    ax = ax + 1;
    mholder = [];
    miholder = [];
    nMF = 50; %Default 50
    nGC = 3000; %Default 3000
    nPk = 1;

for i = 1:100
    
    input = MakeSignal(nPts,nMF,1,options);
    % generate connection matrix to Purkinjie cells
    W2 = rand(nPk, nGC);
    % scale weight matrix so that mean response = mean target
    W1 = zeros(1,nMF);
    for k = 1:nGC
        sel = randsample(nMF, 4);
        W1(k,sel) = 1/4;
    end
    mu = mean(mean(input));
    sigma = mean(std(input,[],2));
    
    alpha = a;
    thr = mu + alpha * sigma;
    
    gc = W1 * input;
    gc = gc - thr;
    gc(gc < 0) = 0;
    
    [sgc,xgc] = max(gc');
    [meow,mix] = sort(xgc);
    
    noiseAmp = [n];
    noise = randn(size(input));
    noise = noise'.*var(input');
    noise = noise'.*noiseAmp;
    input2 = input + noise;
    
    mu = mean(mean(input2));
    sigma = mean(std(input2,[],2));
    thr2 = mu + alpha * sigma;
    
    gc2 = W1 * input2;
    gc2 = gc2 - thr;
    gc2(gc2 < 0) = 0;
    
    [sgc,xgc2] = max(gc2');
    [meow2,mix2] = sort(xgc2);
    
    mholder = [mholder abs(xgc-xgc2)];
    
    [sgc,xgc] = max(input');
    [sgc,xgc2] = max(input2');
    miholder = [miholder abs(xgc-xgc2)];
end

mlem(ax).n = n
mlem(ax).mholder = mholder
mlem(ax).miholder = miholder
end


figure()
s = []
p = []
for i = 1:ax
    s = [s mlem(i).mholder];
    p = [p mlem(i).mholder*0+mlem(i).n];
end

boxplot(s,p,'Symbol','')

figure()
s2 = []
p2 = []
for i = 1:ax
    s2 = [s2 mlem(i).miholder];
    p2 = [p2 mlem(i).miholder*0+mlem(i).n];
    [~,sq(i)] = ttest2(mlem(i).mholder,mlem(i).miholder)
end

boxplot(s2,p2,'Symbol','')

figure()
opt2.scatter = 0;
opt2.SE = 0;
prettyWormPlot(s,p,opt2,[1 0 0])
opt2.scatter = 1;
opt2.SE = 0;
prettyWormPlot(s2,p2,opt2,[0 0 0])
legend(["GC","MF"])
prettify(gcf)
xlabel('Noise Level')
ylabel('Offset')
print(gcf, '-dsvg', sdir + "/Figure_2_D_VarNoise2.svg")

figure()
opt2.scatter = 0;
opt2.SE = 1;
prettyWormPlot(s,p,opt2,[1 0 0])
opt2.scatter = 1;
opt2.SE = 1;
prettyWormPlot(s2,p2,opt2,[0 0 0])
legend(["GC","MF"])
prettify(gcf)
xlabel('Noise Level')
ylabel('Offset')
print(gcf, '-dsvg', sdir + "/Figure_2_D_VarNoise2SE.svg")

%%


mlem = struct();
ax = 0;
for n = [0:.1:1 2:1:5]
    a = 0
    ax = ax + 1;
    mholder = [];
    miholder = [];
    nMF = 50; %Default 50
    nGC = 3000; %Default 3000
    nPk = 1;

for i = 1:100
    
    input = MakeSignal(nPts,nMF,1,options);
    % generate connection matrix to Purkinjie cells
    W2 = rand(nPk, nGC);
    % scale weight matrix so that mean response = mean target
    W1 = zeros(1,nMF);
    for k = 1:nGC
        sel = randsample(nMF, 4);
        W1(k,sel) = 1/4;
    end
    mu = mean(mean(input));
    sigma = mean(std(input,[],2));
    
    alpha = a;
    thr = mu + alpha * sigma;
    
    gc = W1 * input;
    gc = gc - thr;
    gc(gc < 0) = 0;
    
    [sgc,xgc] = max(gc');
    [meow,mix] = sort(xgc);
    
    noiseAmp = [n];
    noise = randn(size(input));
    noise = noise'.*var(input');
    noise = noise'.*noiseAmp;
    input2 = input + noise;
    
    mu = mean(mean(input2));
    sigma = mean(std(input2,[],2));
    thr2 = mu + alpha * sigma;
    
    gc2 = W1 * input2;
    gc2 = gc2 - thr;
    gc2(gc2 < 0) = 0;
    
    [sgc,xgc2] = max(gc2');
    [meow2,mix2] = sort(xgc2);
    
    mholder = [mholder abs(xgc-xgc2)];
    
    [sgc,xgc] = max(input');
    [sgc,xgc2] = max(input2');
    miholder = [miholder abs(xgc-xgc2)];
end

mlem(ax).n = n
mlem(ax).mholder = mholder
mlem(ax).miholder = miholder
end


figure()
s = []
p = []
for i = 1:ax
    s = [s mlem(i).mholder];
    p = [p mlem(i).mholder*0+mlem(i).n];
end

boxplot(s,p,'Symbol','')

figure()
s2 = []
p2 = []
for i = 1:ax
    s2 = [s2 mlem(i).miholder];
    p2 = [p2 mlem(i).miholder*0+mlem(i).n];
    [~,sq(i)] = ttest2(mlem(i).mholder,mlem(i).miholder)
end

boxplot(s2,p2,'Symbol','')

figure()
opt2.scatter = 0;
opt2.SE = 0;
prettyWormPlot(s,p,opt2,[1 0 0])
opt2.scatter = 1;
opt2.SE = 0;
prettyWormPlot(s2,p2,opt2,[0 0 0])
legend(["GC","MF"])
prettify(gcf)
xlabel('Noise Level')
ylabel('Offset')
print(gcf, '-dsvg', sdir + "/Figure_2_D_t0VarNoise2.svg")

figure()
opt2.scatter = 0;
opt2.SE = 1;
prettyWormPlot(s,p,opt2,[1 0 0])
opt2.scatter = 1;
opt2.SE = 1;
prettyWormPlot(s2,p2,opt2,[0 0 0])
legend(["GC","MF"])
prettify(gcf)
xlabel('Noise Level')
ylabel('Offset')
print(gcf, '-dsvg', sdir + "/Figure_2_D_t0VarNoise2SE.svg")

