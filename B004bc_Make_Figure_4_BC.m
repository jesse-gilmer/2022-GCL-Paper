%% HEADER

% B004bc_Make_Figure_4_BC.m
% From: An emergent temporal basis set robustly supports cerebellar time-series learning
% Authors: Jesse I. Gilmer, Michael A. Farries, Zachary Kilpatrick, Ioannis Delis, Abigail L. Person
% 2022, BiorXiv, https://www.biorxiv.org/content/10.1101/2022.01.06.475265v1
%
% This script was written by MF.
% Utilities must be on path.
%
% INPUTS: N/A.
% OUTPUTS: 3 figure files to ".../Figure 4/"


%% Set Network Parameters
clear all
close all
dt = 0.01;

fdir = 'Figure 4/';
mkdir(fdir)

clear mf
mf.N = 50;                  % number of mossy fiber inputs
mf.T = 25;                  % trial duration
mf.mean = 5;                % MF rate mean
mf.sd = 3;                  % MF rate standard deviation
mf.tau = 2;                 % MF rate autocorrelation time

clear gc
gc.N = 500;                   % number of granule cells (set to 0 for direct MF-->Pk network)
gc.conv = 4;                % number of MF inputs innervating each GC
gc.alphaZ = nan;              % GC threshold (z-scored)
gc.thr = nan;               % GC threshold (absolute)

clear target
target.N = 1;               % number of Purkinjie cells
target.mean = 5;            % target rate mean
target.sd = 3;              % target rate SD
target.tau = 2;             % target rate AC time

train.nTr = 300;             % number of training trials
train.rate = 1e-6;          % learning rate
train.batch = false;        % if true, weight changes implemented at end of each trial
train.normalize = true;     % if true, GC-Pk weight matrix adjusted to make initial Pk match target mean
train.negWeights = true;   % if true, allows learning rule to make positive weights negative

noise.sd = 0;               % noise SD (added to Pk activity)
noise.tau = 0;              % noise AC time (added to Pk activity)


%% Assess Performance as a GC thr

alphaZ = -1.4:0.25:1.4;
nReps = 50;
learnCurve = cell(1, length(alphaZ));

for j = 1:length(alphaZ)
    gc.alphaZ = alphaZ(j);
    learnCurve{j} = nan(nReps, train.nTr+1);
    for n = 1:nReps
        learnCurve{j}(n,:) = CbLearn(mf, gc, target, train, noise, dt);
    end  
end 


%% Plot Results
h = figure()
validReps = cell(size(learnCurve));
nValidReps = zeros(size(learnCurve));
for j = 1:length(learnCurve)
    validReps{j} = find(arrayfun(@(n) all(~isnan(learnCurve{j}(n,:))), 1:nReps));
    nValidReps(j) = length(validReps{j});
end

minReps = 10;
validThr = find(nValidReps >= minReps);
alphaZ(validThr);

learnTau1 = nan(size(learnCurve));
learnTau2 = nan(size(learnCurve));
x = (2:301)';
ft = fittype('a1*exp(-x/tau1) + a2*exp(-x/tau2) + c')
for j = validThr
    y = mean(learnCurve{j}(validReps{j},x))';
    f = fit(x, y, ft, 'StartPoint',  [15 1 3 0.5 40], 'Lower', [1 0.01 0 0.01 1]);
    learnTau1(j) = f.tau1;
    learnTau2(j) = f.tau2;
    figure
    scatter(x, y, 'filled')
    hold on
    fitted = f.a1*exp(-x/f.tau1) + f.a2*exp(-x/f.tau2) + f.c;
    plot(x, fitted, 'r')
    title(num2str(alphaZ(j)))
end

axis tight
saveas(h,fdir + "/Figure_4_B.emf")
saveas(h,fdir + "/Figure_4_B.png")

h = figure
scatter(alphaZ(validThr), learnTau1(validThr))
axis tight
saveas(h,fdir + "/Figure_4_C_upper.emf")
saveas(h,fdir + "/Figure_4_C_upper.png")
h = figure
scatter(alphaZ(validThr), learnTau2(validThr))
axis tight
saveas(h,fdir + "/Figure_4_C_lower.emf")
saveas(h,fdir + "/Figure_4_C_lower.png")
