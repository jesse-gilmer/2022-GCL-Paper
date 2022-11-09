clear all;
close all;

%% Make Figure 3

fdir = "Figure 3_bup"
mkdir(fdir)

%% Preprocess:

clearvars -except fdir
load("OU_Datasets_compiled.mat")

% Filter data
for i = 1:length(mRes)-1
    nGC(i) = mRes(i).nGC;
    a(i) = mRes(i).alpha;
    m(i) = (mRes(i).mse/mRes(i).MF_mse);
    mgc(i) = mRes(i).mse;
    mmf(i) = mRes(i).MF_mse;
    ts(i) = mRes(i).gcTS_rel;
    tss(i) = mRes(i).gcTS_weak;
    dim(i) = mRes(i).gcDim;
    mfdim(i) = mRes(i).mfDim;
    cc(i) = mRes(i).gcCorr;
    MFcc(i) =  mRes(i).mfCorr;
    tau(i) = mRes(i).gcTau;
    mtau(i) = mRes(i).inTau;
    v(i) = mRes(i).gcVar;
    lo(i) = mRes(i).gcLoss;
    gTau(i) = mRes(i).gcLearntau;
    mTau(i) = mRes(i).mfLearntau;
    try
        ploss(i) = mRes(i).gcPopSparse;
        exvar(i) = mRes(i).explainedVar;
        PCAind(i) = mRes(i).PCAind1;
        PCAind100(i) = mRes(i).PCAind;
        pSparse(i) = mRes(i).gcPopSparse;
        binSparse(i) = mRes(i).gcBinarySparse;
        marr(i,:) = mRes(i).mse_all;
        mfarr(i,:) = mRes(i).MF_mse_all;
    end
end

%% Panel A

h = figure();
hold on
opt.SE = 0;
opt.scatter = 1;
prettyWormPlot(mgc,a,opt,[1 0 0])
prettyWormPlot([mmf mmf],[ones(1,length(mmf))*min(a) ones(1,length(mmf))*max(a)],opt,[0 0 1])
ylabel('MSE')
yyaxis right
prettyWormPlot(m,a,opt,[0 0 0])
xlabel('Threshold')
ylabel('Relative MSE')
prettify(gcf)
axis square
hold on 
legend({'','GCL','','',"MFs"})

saveas(h,fdir + "/Figure_3_A_MSE.emf")

h = figure();
hold on
opt.SE = 0;
opt.scatter = 1;
prettyWormPlot(mgc,a,opt,[1 0 0])
prettyWormPlot([mmf mmf],[ones(1,length(mmf))*min(a) ones(1,length(mmf))*max(a)],opt,[0 0 1])
ylabel('MSE')
ylim([0 .048])
yyaxis right
prettyWormPlot(m,a,opt,[0 0 0])
ylim([0 2])
xlabel('Threshold')
ylabel('Relative MSE')
prettify(gcf)
axis square
hold on 
legend({'','GCL','','',"MFs"})

saveas(h,fdir + "/Figure_3_A_MSE_zoom_in.emf")

%% Panel B

% This is the cat panel

%% Panel C


clearvars -except fdir
load("nMF_Test_compiled_Aug.mat")

% Filter data
for i = 1:length(mRes)-1
    nGC(i) = mRes(i).nGC;
    nMF(i) = mRes(i).nMF;
    a(i) = mRes(i).alpha;
    m(i) = (mRes(i).mse./mRes(i).MF_mse);
    mgc(i) = mRes(i).mse;
    mmf(i) = mRes(i).MF_mse;
    ts(i) = mRes(i).gcTS_rel;
    tss(i) = mRes(i).gcTS_weak;
    dim(i) = mRes(i).gcDim;
    mfdim(i) = mRes(i).mfDim;
    cc(i) = mRes(i).gcCorr;
    MFcc(i) =  mRes(i).mfCorr;
    tau(i) = mRes(i).gcTau;
    mtau(i) = mRes(i).inTau;
    v(i) = mRes(i).gcVar;
    lo(i) = mRes(i).gcLoss;
    gTau(i) = mRes(i).gcLearntau;
    mTau(i) = mRes(i).mfLearntau;
    try
        ploss(i) = mRes(i).gcPopSparse;
        exvar(i) = mRes(i).explainedVar;
        PCAind(i) = mRes(i).PCAind1;
        PCAind100(i) = mRes(i).PCAind;
        pSparse(i) = mRes(i).gcPopSparse;
        binSparse(i) = mRes(i).gcBinarySparse;
        marr(i,:) = mRes(i).mse_all;
        mfarr(i,:) = mRes(i).MF_mse_all;
    end
end



h = figure();
hold on
opt.SE = 0;
opt.scatter = 1;
prettyWormPlot(mgc,nMF,opt,[1 0 0])
prettyWormPlot(mmf,nMF,opt,[0 0 1])
scatter(50,mean(mgc(nMF==50)),100,'r')
ylabel('MSE')
xlabel('# MFs')
prettify(gcf)
axis square
hold on 
legend({'','GCL','','',"MFs"})

saveas(h,fdir + "/Figure_3_C_nmf_MSE.emf")

h = figure();
hold on
opt.SE = 0;
opt.scatter = 1;
prettyWormPlot(mgc,nMF,opt,[1 0 0])
prettyWormPlot(mmf,nMF,opt,[0 0 1])
scatter(50,mean(mgc(nMF==50)),100,'r')
ylabel('MSE')
xlabel('# MFs')
ylim([0 0.035])
prettify(gcf)
axis square
hold on 
legend({'','GCL','','',"MFs"})

saveas(h,fdir + "/Figure_3_C_zoom_nmf_MSE.emf")

h = figure();
hold on
opt.SE = 0;
opt.scatter = 1;
prettyWormPlot(mgc,nMF,opt,[1 0 0])
prettyWormPlot(mmf,nMF,opt,[0 0 1])
scatter(50,mean(mgc(nMF==50)),100,'r')
scatter(50,mean(mmf(nMF==50)),100,'b')
plot([0 100],[mean(mgc(nMF==50))-std(mgc(nMF==50)) mean(mgc(nMF==50))-std(mgc(nMF==50))],'r--')
plot([0 100],[mean(mgc(nMF==50))+std(mgc(nMF==50)) mean(mgc(nMF==50))+std(mgc(nMF==50))],'r--')
plot([0 100],[mean(mmf(nMF==50))-std(mmf(nMF==50)) mean(mmf(nMF==50))-std(mmf(nMF==50))],'b--')
plot([0 100],[mean(mmf(nMF==50))+std(mmf(nMF==50)) mean(mmf(nMF==50))+std(mmf(nMF==50))],'b--')
ylabel('MSE')
xlabel('# MFs')
ylim([0 0.035])
prettify(gcf)
axis square
hold on 
legend({'','GCL','','',"MFs"})

saveas(h,fdir + "/Figure_3_C_zoom_nmf_MSE_w_std.emf")

%%


clearvars -except fdir
load("nGC_Test_compiled_Aug.mat")

% Filter data
for i = 1:length(mRes)-1
    nGC(i) = mRes(i).nGC;
    nMF(i) = mRes(i).nMF;
    a(i) = mRes(i).alpha;
    m(i) = (mRes(i).mse./mRes(i).MF_mse);
    mgc(i) = mRes(i).mse;
    mmf(i) = mRes(i).MF_mse;
    ts(i) = mRes(i).gcTS_rel;
    tss(i) = mRes(i).gcTS_weak;
    dim(i) = mRes(i).gcDim;
    mfdim(i) = mRes(i).mfDim;
    cc(i) = mRes(i).gcCorr;
    MFcc(i) =  mRes(i).mfCorr;
    tau(i) = mRes(i).gcTau;
    mtau(i) = mRes(i).inTau;
    v(i) = mRes(i).gcVar;
    lo(i) = mRes(i).gcLoss;
    gTau(i) = mRes(i).gcLearntau;
    mTau(i) = mRes(i).mfLearntau;
    try
        ploss(i) = mRes(i).gcPopSparse;
        exvar(i) = mRes(i).explainedVar;
        PCAind(i) = mRes(i).PCAind1;
        PCAind100(i) = mRes(i).PCAind;
        pSparse(i) = mRes(i).gcPopSparse;
        binSparse(i) = mRes(i).gcBinarySparse;
        marr(i,:) = mRes(i).mse_all;
        mfarr(i,:) = mRes(i).MF_mse_all;
    end
end



h = figure();
hold on
opt.SE = 0;
opt.scatter = 1;
prettyWormPlot(mgc,nGC,opt,[1 0 0])
prettyWormPlot([mmf mmf],[ones(1,length(mmf))*min(nGC) ones(1,length(mmf))*max(nGC)],opt,[0 0 1])
scatter(3000,mean(mgc(nGC==3000)),100,'r')
ylabel('MSE')
xlabel('# GCs')
prettify(gcf)
axis square
hold on 
legend({'','GCL','','',"MFs"})


saveas(h,fdir + "/Figure_3_D_nGC_MSE.emf")

h = figure();
hold on
opt.SE = 0;
opt.scatter = 1;
prettyWormPlot(mgc,nGC,opt,[1 0 0])
prettyWormPlot([mmf mmf],[ones(1,length(mmf))*min(nGC) ones(1,length(mmf))*max(nGC)],opt,[0 0 1])
scatter(3000,mean(mgc(nGC==3000)),100,'r')
ylabel('MSE')
xlabel('# GCs')
xlim([-10 5100])
ylim([0 0.025])
prettify(gcf)
axis square
hold on 
legend({'','GCL','','',"MFs"})


saveas(h,fdir + "/Figure_3_D_nGC_MSE_zoomed.emf")



h = figure();
hold on
opt.SE = 0;
opt.scatter = 1;
prettyWormPlot(mgc,log10(nGC),opt,[1 0 0])
prettyWormPlot([mmf mmf],[ones(1,length(mmf))*min(log10(nGC)) ones(1,length(mmf))*max(log10(nGC))],opt,[0 0 1])
scatter(log10(3000),mean(mgc(nGC==3000)),100,'r')
ylabel('MSE')
xlabel('log_1_0 # GCs')
prettify(gcf)
axis square
hold on 
legend({'','GCL','','',"MFs"})


saveas(h,fdir + "/Figure_3_D_log10nGC_MSE.emf")

h = figure();
hold on
opt.SE = 0;
opt.scatter = 1;
prettyWormPlot(mgc,log(nGC),opt,[1 0 0])
prettyWormPlot([mmf mmf],[ones(1,length(mmf))*min(log(nGC)) ones(1,length(mmf))*max(log(nGC))],opt,[0 0 1])
scatter(log(3000),mean(mgc(nGC==3000)),100,'r')
ylabel('MSE')
xlabel('ln # GCs')
prettify(gcf)
axis square
hold on 
legend({'','GCL','','',"MFs"})


saveas(h,fdir + "/Figure_3_D_ln_nGC_MSE.emf")