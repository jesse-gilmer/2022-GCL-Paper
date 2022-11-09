%% HEADER

% B003c_Make_Figure_3_C.m
% From: An emergent temporal basis set robustly supports cerebellar time-series learning
% Authors: Jesse I. Gilmer, Michael A. Farries, Zachary Kilpatrick, Ioannis Delis, Abigail L. Person
% 2022, BiorXiv, https://www.biorxiv.org/content/10.1101/2022.01.06.475265v1
%
% This script was written by JG.
% Data used in this script comes from https://zenodo.org/record/5140528#.YuA103bMJD9
%
% INPUTS: N/A.
%
% OUTPUTS: 3 figure panels to /Figure 3/

%%
clear all
close all

T = 25;
mRes = struct();

nMF = 100;
nGC = 6000;
conv = 4;

mu = 1;
sigma = .33;

options.mu = mu;
options.sigma = sigma;
options.lrnRate = .0002;
options.nTrials = 2000;

t = imread('Datasets/cat_g.png');
t = t(:,:,1);
tel = numel(t);

nPts = tel;

AA = 0;
for A = [-.75 0 .75]
%% Project params:
AA = AA+ 1;
alpha = A;
options.tau = T;


%% Make input signal:
input = MakeSignal(nPts,nMF,1,options);

%% Make target signal:
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

target = normalize(target,'range');
[results] = RunFilter(input,target,gc,options);
cat_copy = results.W * gc;
cat_copy = reshape(cat_copy, size(t));

figure(1); hold on;
plot(results.mse,'r')

figure(99999)
subplot(3,3,AA)
imshow(cat_copy)
title("Threshold: " + A + ", MSE: " + results.mse(end))

[results] = RunFilter(input,target,input,options);
cat_copy = results.W * input;
cat_copy = reshape(cat_copy, size(t));

figure(1); hold on;
plot(results.mse,'k')

figure(999992)
subplot(3,3,AA)
imshow(cat_copy)
title("Threshold: " + A + ", MSE: " + results.mse(end))
drawnow
end

mkdir("Figure 3");

figure(99999)
hold on;
print(gcf, '-dpng', "Figure 3/Figure_3_C_GC.png")

figure(999992)
hold on;
print(gcf, '-dpng', "Figure 3/Figure_3_C_MF.png")

figure();
imshow(t)
print(gcf, '-dpng', "Figure 3/Figure_3_C_TF.png")