%% HEADER

% B001d_Make_Figure_1_D.m
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


%% Raster / PSTH plots centered on cue for example neurons.
clear all
close all

fdir = 'Figure 1';
mkdir(fdir)

nMF = 50; %Default 50
nGC = 300; %Default 3000

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
options.tau = 100;

input = MakeSignal(nPts,nMF,1,options);
input = normalize(input','range')';

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
saveas(h,fdir + "/Figure_1_D_upper.emf")
saveas(h,fdir + "/Figure_1_D_upper.png")

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
saveas(h,fdir + "/Figure_1_D_lower.emf")
saveas(h,fdir + "/Figure_1_D_lower.png")
 %%

nMF = 5000; %Default 50
% Project params:
options.tau = 100;
input = MakeSignal(nPts,nMF,1,options);
input = normalize(input','range')';

ft = fittype('a*exp(-x/tau)')
for i = 1:nMF
    c = xcorr(input(i,:));  
    c = c(ceil(end/2):end);
    cl = [1:length(c)];
    f = fit(cl', c', ft, 'StartPoint',  [c(1) 100]);
    fitted = f.a*exp(-cl/f.tau);
    tau(i) = f.tau;
end

h = figure()
boxplot(tau)
text(1.1,mean(tau),"The mean decay tau of OU is: " + mean(tau))
text(1.1,mean(tau)-std(tau),"The st. dev of the decay tau of OU is: " + std(tau))
xlim([0 4])
xlabel('PN')
ylabel('Decay Taus')
prettify(gcf)
saveas(h,fdir + "/Figure_1_Dxx_lower.emf")
saveas(h,fdir + "/Figure_1_Dxx_lower.png")

disp("The mean decay tau of OU is: " + mean(tau))
disp("The st. dev of the decay tau of OU is: " + std(tau))

save(fdir + "/OU_tau.mat",'tau')

