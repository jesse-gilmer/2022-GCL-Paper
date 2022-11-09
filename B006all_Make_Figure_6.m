% EDIT THIS!


clear all;
close all;

%% Make Figure 1

fdir = "Figure 6"
mkdir(fdir)

%% Panel A

% This is an illustration.

%% Panel B

close all;

% i) MF Signals examples

% Load the creation options.
options.mu = 1;
options.sigma = .25;
options.lrnRate = 1e-5;
options.nTrials = 150;
options.tau = 10;

nPts = 1000;
nMF = 50;

input = MakeSignal(nPts,nMF,1,options);

h = figure;
prettify(gcf)
hold on;
for i = 1:5
    plot(input(i,:)+i)
end

saveas(h,fdir + "/Figure_1_i_MFexample.emf")
saveas(h,fdir + "/Figure_1_i_MFexample.png")

% ii) MF Signal Heatmap

h = figure;
prettify(gcf)
imagesc(input)
colormap(turbo)
colorbar
xlabel('Time (ms)')
ylabel('Unit #')
a = colorbar;

set(get(a,'label'),'string','Rate (AU)');
saveas(h,fdir + "/Figure_1_ii_MFrate.emf")
saveas(h,fdir + "/Figure_1_ii_MFrate.png")

h = figure();
imagesc(input)

% iii) GC Signal examples
nGC = 500;
conv = 4;

W1 = zeros(nGC, nMF);
for k = 1:nGC
    sel = randsample(nMF, conv);
    W1(k,sel) = 1/conv;
end
mu = mean(mean(input));
sigma = mean(std(input,[],2));
thr = mu + .8 * sigma;

gc = W1 * input;
gc = gc - thr;
gc(gc < 0) = 0;

for i = 1:nGC
    min_s = min(find(gc(i,:) == max(gc(i,:))));
    if length(min_s) <= 0
        min_s = 1000;
    end
    y(i) = min_s;
end

[~,c] = sort(y);
gc = gc(flipud(c),:);

gc2 = gc;
gc2(sum(gc2'>0)>nPts/20,:) = [];
gc2(sum(gc2'>0)<=5,:) = [];

h = figure;
prettify(gcf)
hold on;
for i = 1:7
    plot(gc2(round(i*size(gc2,1)/8),:)-i/12)
end

saveas(h,fdir + "/Figure_1_iii_GCexample.emf")
saveas(h,fdir + "/Figure_1_iii_GCexample.png")

% iv) GC Signal Heatmap

gc(sum(gc'>0)==0,:) = [];

h = figure;
prettify(gcf)
imagesc(gc)
colormap(turbo)
colorbar
xlabel('Time (ms)')
ylabel('Unit #')
a = colorbar;

set(get(a,'label'),'string','Rate (AU)');
saveas(h,fdir + "/Figure_1_iv_GCrate.emf")
saveas(h,fdir + "/Figure_1_iv_GCrate.png")

% v) Sparce GC Signal Heatmap

h = figure;
prettify(gcf)
imagesc(gc2)
colormap(turbo)
colorbar
xlabel('Time (ms)')
ylabel('Unit #')
a = colorbar;

set(get(a,'label'),'string','Rate (AU)');
saveas(h,fdir + "/Figure_1_v_GCrate_sparse.emf")
saveas(h,fdir + "/Figure_1_v_GCrate_sparse.png")

%% Panel D Three Signals....

nPts = 1000;
nMF = 30;

options.mu = 1;
options.sigma = .25;
options.lrnRate = 1e-5;
options.nTrials = 150;
options.tau = 150;

input = MakeSignal(nPts,nMF,1,options);
nGC = 300;
conv = 4;

W1 = zeros(nGC, nMF);
for k = 1:nGC
    sel = randsample(nMF, conv);
    W1(k,sel) = 1/conv;
end

mu = mean(mean(input));
sigma = mean(std(input,[],2));
for THR = [-1 0 1.0]
    thr = mu + (THR * sigma);
    
    gc = W1 * input;
    gc = gc - thr;
    gc(gc < 0) = 0;
    
    y = [];
    for i = 1:nGC
        min_s = min(find(gc(i,:) == max(gc(i,:))));
        if length(min_s) <= 0
            min_s = 1000;
        end
        y(i) = min_s;
    end
    
    [~,c] = sort(y);
    gc = gc(flipud(c),:);
    

gc2 = gc;
% gc2(sum(gc2'>0)>nPts/20,:) = [];
gc2(sum(gc2'>0)<=5,:) = [];
    
%     gc2 = normalize(gc2,'range')/2;
    h = figure;
    prettify(gcf)
    imagesc(gc2)
    colormap(turbo)
    colorbar
    xlabel('Time (ms)')
    ylabel('Unit #')
    a = colorbar;
    
    saveas(h,fdir + "/Figure_1C_MapGCexample"+THR+".emf")
    saveas(h,fdir + "/Figure_1C_MapGCexample"+THR+".png")
end

% end
%
% saveas(h,fdir + "/Figure_1_iii_GCexample.emf")
% saveas(h,fdir + "/Figure_1_iii_GCexample.png")


%% More... concrete version
w = gausswin(20);
w = w/sum(w);

a = zeros(6,150);
a(1,40:50) = 1;
a(2,50:60) = 1;
a(3,60:70) = 1;
a(4,70:80) = 1;
a(5,80:90) = 1;
a(6,90:100) = 1;
b = filtfilt(w,1,a')';

h = figure;
prettify(gcf)
hold on;
plot(b'-[1:6],'k')
plot([60 60],[0 -6],'r')
plot([90 90],[0 -6],'r')

saveas(h,fdir + "/Figure_1C_GCexample_GCmedium.emf")
saveas(h,fdir + "/Figure_1C_GCexample_GCmedium.png")

h = figure;
prettify(gcf)
hold on;

c = zeros(1,150);
c(60:90) = 1;
plot(c,'k')

saveas(h,fdir + "/Figure_1C_GCexample_TFmedium.emf")
saveas(h,fdir + "/Figure_1C_GCexample_TFmedium.png")


h = figure;
prettify(gcf)
hold on;


options.lrnRate = 1e-2;
options.nTrials = 2000;
[results] = RunFilter(b,c,b,options);

plot(c,'k')
plot([b'.*results.W],'color',[0.5 0.5 0.5])
plot(sum([b'.*results.W]'),'r--')


saveas(h,fdir + "/Figure_1C_GCexample_Allmedium.emf")
saveas(h,fdir + "/Figure_1C_GCexample_Allmedium.png")

%% too much
w = gausswin(20);
w = w/sum(w);

a = zeros(6,150);
a(1,40:50) = 1;
a(2,50:60) = 1;
a(3,60:70) = 1;
a(4,70:80) = 1;
a(5,80:90) = 1;
a(6,90:100) = 1;

a(1,20:30) = 1;
a(2,90:100) = 1;
a(3,120:145) = 1;
a(4,120:130) = 1;
a(5,40:50) = 1;
a(6,140:145) = 1;
b = filtfilt(w,1,a')';

h = figure;
prettify(gcf)
hold on;
plot(b'-[1:6],'k')
plot([60 60],[0 -6],'r')
plot([90 90],[0 -6],'r')

saveas(h,fdir + "/Figure_1C_GCexample_GChigh.emf")
saveas(h,fdir + "/Figure_1C_GCexample_GChigh.png")


h = figure;
prettify(gcf)
hold on;

c = zeros(1,150);
c(60:90) = 1;
plot(c,'k')

saveas(h,fdir + "/Figure_1C_GCexample_TFhigh.emf")
saveas(h,fdir + "/Figure_1C_GCexample_TFhigh.png")


h = figure;
prettify(gcf)
hold on;

options.lrnRate = 1e-2;
options.nTrials = 2000;
        
[results] = RunFilter(b,c,b,options);

plot(c,'k')
plot([b'.*results.W],'color',[0.5 0.5 0.5])
plot(sum([b'.*results.W]'),'r--')

saveas(h,fdir + "/Figure_1C_GCexample_Allhigh.emf")
saveas(h,fdir + "/Figure_1C_GCexample_Allhigh.png")

%% too high
w = gausswin(15);
w = w/sum(w);

a = zeros(6,150);
a(1,40:42) = 1;
a(2,50:60) = 1;
a(3,60:65) = 1;
a(4,72:72) = 1;
a(6,90:100) = 1;
b = filtfilt(w,1,a')';

h = figure;
prettify(gcf)
hold on;
plot(b'-[1:6],'k')
plot([60 60],[0 -6],'r')
plot([90 90],[0 -6],'r')

saveas(h,fdir + "/Figure_1C_GCexample_GClow.emf")
saveas(h,fdir + "/Figure_1C_GCexample_GClow.png")


h = figure;
prettify(gcf)
hold on;

c = zeros(1,150);
c(60:90) = 1;
plot(c,'k')

saveas(h,fdir + "/Figure_1C_GCexample_TFlow.emf")
saveas(h,fdir + "/Figure_1C_GCexample_TFlow.png")


h = figure;
prettify(gcf)
hold on;

options.lrnRate = 1e-2;
options.nTrials = 2000;
        
[results] = RunFilter(b,c,b,options);


plot(c,'k')
plot([b'.*results.W],'color',[0.5 0.5 0.5])
plot(sum([b'.*results.W]'),'r--')


saveas(h,fdir + "/Figure_1C_GCexample_Alllow.emf")
saveas(h,fdir + "/Figure_1C_GCexample_Alllow.png")




%% Panel E preprocessing
clearvars -except fdir
load("PCA_test_3000F_retest_compiled_Aug.mat")
FilterData

%% Panel E: Dimensionality Enhancement

h = figure();
hold on
opt.SE = 0;
opt.scatter = 1;
prettyWormPlot(dim,a,opt,[1 0 0])
prettyWormPlot([mfdim mfdim],[ones(1,length(mfdim))*min(a) ones(1,length(mfdim))*max(a)],opt,[0 0 1])
xlabel('Threshold')
ylabel('Dimensionality')
prettify(gcf)
axis square
hold on
legend({'','GCL','','',"MFs"})

saveas(h,fdir + "/Figure_1_C_Dim.emf")
saveas(h,fdir + "/Figure_1_C_Dim.png")

%% Panel D: Decorrelation

h = figure();
hold on
opt.SE = 0;
opt.scatter = 1;
prettyWormPlot(cc,a,opt,[1 0 0])
prettyWormPlot([MFcc MFcc],[ones(1,length(MFcc))*min(a) ones(1,length(MFcc))*max(a)],opt,[0 0 1])
xlabel('Threshold')
ylabel('Correlation')
prettify(gcf)
axis square
hold on
legend({'','GCL','','',"MFs"})

saveas(h,fdir + "/Figure_1_D_Corr.emf")
saveas(h,fdir + "/Figure_1_D_Corr.png")

%% Panel E: Temporal Properties

h = figure();
hold on
opt.SE = 0;
opt.scatter = 1;
prettyWormPlot(tau,a,opt,[1 0 0])
prettyWormPlot([mtau mtau],[ones(1,length(mtau))*min(a) ones(1,length(mtau))*max(a)],opt,[0 0 1])
xlabel('Threshold')
ylabel('Temporal Decay')
prettify(gcf)
axis square
hold on
legend({'','GCL','','',"MFs"})

saveas(h,fdir + "/Figure_1_E_Tau.emf")
saveas(h,fdir + "/Figure_1_E_Tau.png")

%% Panel F: STS

PA = parula(3);

h = figure();
hold on
opt.SE = 0;
opt.scatter = 1;
prettyWormPlot(ts,a,opt,PA(1,:))
ylabel('STS')
yyaxis right
ylabel('PCA Contribution')
prettyWormPlot(PCAind,a,opt,PA(2,:))
xlabel('Threshold')

prettify(gcf)
axis square

saveas(h,fdir + "/Figure_1_F_STS.emf")
saveas(h,fdir + "/Figure_1_F_STS.png")


%% Panel G: Lossiness

PA = prism(6);

h = figure();
hold on
opt.SE = 0;
opt.scatter = 1;
prettyWormPlot(1-lo,a,opt,PA(1,:))
prettyWormPlot(1-ploss,a,opt,PA(2,:))
prettyWormPlot(binSparse,a,opt,PA(5,:))
xlabel('Threshold')
ylabel('Lossiness')
prettify(gcf)
axis square
hold on
legend({'','Temporal Loss','','',"Population Loss",'','',"Mean Temporal Cover"})

saveas(h,fdir + "/Figure_1_H_Loss.emf")
saveas(h,fdir + "/Figure_1_H_Loss.png")

%% Panel H: PCA Contribution

h = figure();
hold on
opt.SE = 0;
opt.scatter = 1;
prettyWormPlot(normalize(exvar,'range'),a,opt,[0 0 1])
xlabel('Threshold')
ylabel('Input Variance Retained')
prettify(gcf)
axis square

saveas(h,fdir + "/Figure_1_H_VarianceEx.emf")
saveas(h,fdir + "/Figure_1_H_VarianceEx.png")
%% Panel I: Variability

h = figure();
hold on
opt.SE = 0;
opt.scatter = 0;
prettyWormPlot(v,a,opt,[0 0 0])
xlabel('Threshold')
ylabel('Variance per GC')
prettify(gcf)
axis square

saveas(h,fdir + "/Figure_1_I_Variance.emf")
saveas(h,fdir + "/Figure_1_I_Variance.png")

%% Panel Sup 1: Variable Correlations

Exps = normalize([dim; ts;  PCAind; v; cc; tau; 1-lo; 1-ploss; binSparse; exvar]','range')';

Names = ["Dim", "STS", "PCA", "Var", "Corr", "Tau",  ...
    "T. Loss", "P. Loss",  "T. Cover", "Var Ret.", ]


h = figure('units','normalized','outerposition',[0 0 1 1])
q = 1;
for i = 1:size(Exps,1)
    for j = 1:i
        subplot(size(Exps,1),size(Exps,1),q)
        scatter(Exps(i,:),Exps(j,:),1,'k')
        ylabel(Names(i))
        xlabel(Names(j))
        xlim([-.1 1.1])
        ylim([-.1 1.1])
        q = q + 1;
    end
    for j = i+1:size(Exps,1)
        q = q + 1;
    end
end
set(gcf,'color','w')
saveas(h,fdir + "/Figure_S1_I_Cov.emf")
saveas(h,fdir + "/Figure_S1_I_Cov.png")