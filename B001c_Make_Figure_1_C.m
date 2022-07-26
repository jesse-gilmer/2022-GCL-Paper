%% HEADER

% B001c_Make_Figure_1_C.m
% From: An emergent temporal basis set robustly supports cerebellar time-series learning
% Authors: Jesse I. Gilmer, Michael A. Farries, Zachary Kilpatrick, Ioannis Delis, Abigail L. Person
% 2022, BiorXiv, https://www.biorxiv.org/content/10.1101/2022.01.06.475265v1
%
% This script was written by JG.
% Data used in this script comes from https://zenodo.org/record/5140528#.YuA103bMJD9
%
% INPUTS: This requires "Datasets/D002_EMG_Data.mat" 
%
% OUTPUTS: 2 figure files to ".../Figure 1/"

%% 
clear all
fdir = 'Figure 1';
mkdir(fdir)

floc = "Datasets/D002_EMG_Data.mat";
load(floc)
kin = data.KIN;
emg = data.EMG;

nPts = 1000; %Default 1000
i = 1;

% Project params:
alpha = 0;
options.tau = 100;

w = gausswin(100);
w = w/sum(w);
a = emg;
a =  padarray(a,[1000 0],'replicate','both');
a = filtfilt(w,1,abs(a));
a = a(1001:end-1000,:);
a = normalize(a,'range');
emg = [a];

mu = 4;
sigma = 1;

options.mu = mu;
options.sigma = sigma;
options.lrnRate = 1e-5;
options.nTrials = 1000; %Default 1000

nPts = 300; %Default 1000
i = 1;

% Project params:
options.tau = 100;

input = [emg'];

nMF = size(input,1); %Default 50
nGC = 300; %Default 3000
nPk = 1;
conv = 4;

W1 = zeros(1,nMF);
for k = 1:nGC
    sel = randsample(nMF, 4);
    W1(k,sel) = 1/4;
    %     W1(k,sel) = rand(1,4);
end

mu = mean(mean(input));
sigma = mean(std(input,[],2));


alpha = 0
thr = mu + alpha * sigma;

gc = W1 * input;
gc = gc - thr;
gc(gc < 0) = 0;


[sgc,xi] = max(input');
[mval,mix] = sort(xi);

h = figure(); hold on
colormap turbo
imagesc(flipud(input(mix,:)))
set(gcf,'color','w')
xlabel('time (ms)')
ylabel('MF unit')
prettify(gcf)
colorbar()
axis tight
saveas(h,fdir + "/Figure_1_C_upper.emf")
saveas(h,fdir + "/Figure_1_C_upper.png")

[sgc,xgc] = max(gc');
[mval,mix] = sort(xgc);

h = figure(); hold on
colormap turbo
imagesc(flipud(gc(mix,:)))
set(gcf,'color','w')
xlabel('time (ms)')
ylabel('GC unit')
prettify(gcf)
colorbar()
axis tight
saveas(h,fdir + "/Figure_1_C_lower.emf")
saveas(h,fdir + "/Figure_1_C_lower.png")
