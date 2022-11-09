%% HEADER

% B001b_Make_Figure_1_B.m
% From: An emergent temporal basis set robustly supports cerebellar time-series learning
% Authors: Jesse I. Gilmer, Michael A. Farries, Zachary Kilpatrick, Ioannis Delis, Abigail L. Person
% 2022, BiorXiv, https://www.biorxiv.org/content/10.1101/2022.01.06.475265v1
%
% This script was written by JG.
% Data used in this script comes from https://zenodo.org/record/5140528#.YuA103bMJD9
%
% INPUTS: This requires "Datasets/D001_JG_PN_Edited_Data.mat" 
%
% OUTPUTS: 2 figure files to ".../Figure 1/"

%% Clear workspace:
close all;
clear all;

%% Load data:
load('Datasets/D001_JG_PN_Edited_Data.mat')

%% Make output directory:
fdir = "Figure 1/";
mkdir(fdir);

%% Process data:
cells4raster = [1:128]
lims_raster = [0 500];
lims_pulse = [-50 50];
w_fr_pulse = 5;

rh = [];
ix = 1;
for i = 1:50
    ii = cells4raster(i);
    [r_tmp t_tmp] = firing_rate({cell2mat(t_spk_peri_lift{ii}')},10,1,lims_raster+[-200 200],'gaussian');
    rh(i,:) = r_tmp;
    ix = ix + 1;
end

% normalize inputs by range:
input = normalize(rh','range')';

% GCL properties:
nMF = size(input,1); %Default 50
nGC = 300; %Default 300
nPk = 1; % number of PKJ cells.
conv = 4;

% MF to GC projection
W1 = zeros(1,nMF);
for k = 1:nGC
    sel = randsample(nMF, 4);
    W1(k,sel) = 1/4;
    %     W1(k,sel) = rand(1,4);
end

mu = mean(mean(input));
sigma = mean(std(input,[],2));
alpha = 0;
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
saveas(h,fdir + "/Figure_1_B_upper.emf")
saveas(h,fdir + "/Figure_1_B_upper.png")

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
saveas(h,fdir + "/Figure_1_B_lower.emf")
saveas(h,fdir + "/Figure_1_B_lower.png")

%%
ft = fittype('a*exp(-x/tau)')
for i = 1:nMF
    c = xcorr(input(i,:));  
    c = c(end/2:end);
    cl = [1:length(c)];
    f = fit(cl', c', ft, 'StartPoint',  [c(1) 100]);
    tau(i) = f.tau;
end

h= figure()
boxplot(tau)
text(1.1,mean(tau),"The mean decay tau of PN is: " + mean(tau))
text(1.1,mean(tau)-std(tau),"The st. dev of the decay tau of PN is: " + std(tau))
xlim([0 4])
xlabel('EMG')
ylabel('Decay Taus')
prettify(gcf)
saveas(h,fdir + "/Figure_1_Bxx_lower.emf")
saveas(h,fdir + "/Figure_1_Bxx_lower.png")

disp("The mean decay tau of PN is: " + mean(tau))
disp("The st. dev of the decay tau of PN is: " + std(tau))

save(fdir + "/PN_tau.mat",'tau')
