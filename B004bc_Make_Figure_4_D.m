%% HEADER

% B004bc_Make_Figure_4_D.m
% From: An emergent temporal basis set robustly supports cerebellar time-series learning
% Authors: Jesse I. Gilmer, Michael A. Farries, Zachary Kilpatrick, Ioannis Delis, Abigail L. Person
% 2022, BiorXiv, https://www.biorxiv.org/content/10.1101/2022.01.06.475265v1
%
% This script was written by MF and JG.
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

train.nTr = 300;            % number of training trials
train.rate = 1e-6;          % learning rate
train.batch = false;        % if true, weight changes implemented at end of each trial
train.normalize = true;     % if true, GC-Pk weight matrix adjusted to make initial Pk match target mean
train.negWeights = true;    % if true, allows learning rule to make positive weights negative

noise.sd = 0;             % noise SD (added to Pk activity)
noise.tau = 0;              % noise AC time (added to Pk activity)


%% Assess Performance as a GC thr

saveTau1 = [];
saveTau2 = [];

ii = 1;
for tr = 1:10
for Z = 0
alphaZ = Z;
nReps = 10;
learnCurve = cell(1, length(alphaZ));
noise.sd = 0;  

for j = 1:length(alphaZ)
    gc.alphaZ = alphaZ(j);
    learnCurve{j} = nan(nReps, train.nTr+1);
    for n = 1:nReps
        learnCurve{j}(n,:) = CbLearn_mfNoise(mf, gc, target, train, noise, dt);
    end  
end 

learnCurvelowN = cell(1, length(alphaZ));
noise.sd = .25;      

for j = 1:length(alphaZ)
    gc.alphaZ = alphaZ(j);
    learnCurvelowN{j} = nan(nReps, train.nTr+1);
    for n = 1:nReps
        learnCurvelowN{j}(n,:) = CbLearn_mfNoise(mf, gc, target, train, noise, dt);
    end  
end 

learnCurvehighN = cell(1, length(alphaZ));
noise.sd = 1;      

for j = 1:length(alphaZ)
    gc.alphaZ = alphaZ(j);
    learnCurvehighN{j} = nan(nReps, train.nTr+1);
    for n = 1:nReps
        learnCurvehighN{j}(n,:) = CbLearn_mfNoise(mf, gc, target, train, noise, dt);
    end  
end 

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
    fitted = f.a1*exp(-x/f.tau1) + f.a2*exp(-x/f.tau2) + f.c;
end
fitnonoise = fitted;

validReps = cell(size(learnCurve));
nValidReps = zeros(size(learnCurve));
for j = 1:length(learnCurve)
    validReps{j} = find(arrayfun(@(n) all(~isnan(learnCurve{j}(n,:))), 1:nReps));
    nValidReps(j) = length(validReps{j});
end

% learnCurvelowN
validThr = find(nValidReps >= minReps);

learnTau1lowN = nan(size(learnCurvelowN));
learnTau2lowN = nan(size(learnCurvelowN));
x = (2:301)';
ft = fittype('a1*exp(-x/tau1) + a2*exp(-x/tau2) + c')
for j = validThr
    y = mean(learnCurvelowN{j}(validReps{j},x))';
    f = fit(x, y, ft, 'StartPoint',  [15 1 3 0.5 40], 'Lower', [1 0.01 0 0.01 1]);
    learnTau1lowN(j) = f.tau1;
    learnTau2lowN(j) = f.tau2;
    fitted = f.a1*exp(-x/f.tau1) + f.a2*exp(-x/f.tau2) + f.c;
end
fitlownoise = fitted;

validReps = cell(size(learnCurve));
nValidReps = zeros(size(learnCurve));
for j = 1:length(learnCurve)
    validReps{j} = find(arrayfun(@(n) all(~isnan(learnCurve{j}(n,:))), 1:nReps));
    nValidReps(j) = length(validReps{j});
end

% highN 
validThr = find(nValidReps >= minReps);
alphaZ(validThr);

learnTau1highN = nan(size(learnCurvehighN));
learnTau2highN = nan(size(learnCurvehighN));
x = (2:301)';
ft = fittype('a1*exp(-x/tau1) + a2*exp(-x/tau2) + c')
jc = learnCurvehighN{j}
for j = validThr
    y = mean(learnCurvehighN{j}(validReps{j},x))';
    f = fit(x, y, ft, 'StartPoint',  [15 1 3 0.5 40], 'Lower', [1 0.01 0 0.01 1]);
    learnTau1highN(j) = f.tau1;
    learnTau2highN(j) = f.tau2;
    fitted = f.a1*exp(-x/f.tau1) + f.a2*exp(-x/f.tau2) + f.c;
end
fithighnoise = fitted;

saveTau1(ii,:) = [alphaZ 0 learnTau1 .1 learnTau1lowN 1 learnTau1highN];
saveTau2(ii,:) = [alphaZ 0 learnTau2 .1 learnTau2lowN 1 learnTau2highN];
ii = ii + 1;

% figure(111)
% hold off
% scatter(saveTau1(:,2), saveTau1(:,3),10,'r')
% hold on
% scatter(saveTau1(:,4), saveTau1(:,5),10,'g')
% scatter(saveTau1(:,6), saveTau1(:,7),10,'b')
% xlabel('Threshold')
% ylabel('Tau')
% legend(["no noise", "low noise", "high noise"])
% 
% figure(112)
% hold off
% scatter(saveTau2(:,2), saveTau2(:,3),10,'r')
% hold on
% scatter(saveTau2(:,4), saveTau2(:,5),10,'g')
% 
% scatter(saveTau2(:,6), saveTau2(:,7),10,'b')
% xlabel('Threshold')
% ylabel('Tau')
% legend(["no noise", "low noise", "high noise"])
% drawnow

[tr]

end
end

h = figure
hold on
decoder = saveTau1(:,[2 4 6]);
decoder = decoder(:);
encoder = saveTau1(:,[3 5 7]);
encoder = encoder(:);
boxplot(encoder,decoder)

[a,b,c,d] = ttest2(saveTau1(:,[3]),saveTau1(:,[5]));
if (b <= .05)
    plot([1.1  1.9],[max(encoder) max(encoder)], 'r')
    scatter([1.5],[max(encoder)],200, 'w','filled')
    scatter([1.5],[max(encoder)], 'r*')
end
disp("no noise vs low noise t-test, tau2 p = " + b)
[a,b,c,d] = ttest2(saveTau1(:,[3]),saveTau1(:,[7]));
if (b <= .05)
    plot([2.1  2.9],[max(encoder) max(encoder)], 'r')
    scatter([2.5],[max(encoder)],200, 'w','filled')
    scatter([2.5],[max(encoder)], 'r*')
end
disp("no noise vs high noise t-test, tau2 p = " + b)
xlabel('noise S.D.')
ylabel('tau 1 (fast)')
saveas(h,fdir + "/Figure_4_D_upper.emf")
saveas(h,fdir + "/Figure_4_D_upper.png")

h = figure
decoder = saveTau2(:,[2 4 6]);
decoder = decoder(:);
encoder = saveTau2(:,[3 5 7]);
encoder = encoder(:);
boxplot(encoder,decoder)

[a,b,c,d] = ttest2(saveTau2(:,[3]),saveTau2(:,[5]));
disp("no noise vs low noise t-test, tau2 p = " + b)
hold on
if (b <= .05/sqrt(4))
    plot([1.1  1.9],[max(encoder) max(encoder)], 'r')
    scatter([1.5],[max(encoder)],200, 'w','filled')
    scatter([1.5],[max(encoder)], 'r*')
end

[a,b,c,d] = ttest2(saveTau2(:,[3]),saveTau2(:,[7]));
if (b <= .05/sqrt(4))
    plot([2.1  2.9],[max(encoder) max(encoder)], 'r')
    scatter([2.5],[max(encoder)],200, 'w','filled')
    scatter([2.5],[max(encoder)], 'r*')
end
disp("no noise vs high noise t-test, tau2 p = " + b)
xlabel('noise S.D.')
ylabel('tau 2 (slow)')
saveas(h,fdir + "/Figure_4_D_lower.emf")
saveas(h,fdir + "/Figure_4_D_lower.png")


mean(saveTau1(:,[3 5 7]))
std(saveTau1(:,[3 5 7]))

mean(saveTau2(:,[3 5 7]))
std(saveTau2(:,[3 5 7]))
