%% Make Figure 7 with Ridge methods.
close all;
clear all;

fdir = "Figure 7 Panels"
mkdir(fdir)

%% VOR-like:
clearvars -except fdir
load("Cosine_Data4_compiled_Aug.mat")

%% Process:
for i = 1:length(mRes)
    nGC(i) = mRes(i).nGC;
    a(i) = mRes(i).alpha;
    m(i) = (mRes(i).mse/(mRes(i).MF_mse+.0001));
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
        O(i) = mRes(i).offset;
    end
end

%% A: Examples of Cosines
nPts = 250;

for Oi = [.1 pi/2 pi]
ls = linspace(0,2*pi,nPts);
input = cos(ls+Oi)+1;
            
target = cos(ls)+1;
h = figure();
hold off
plot(target,'k','linewidth',2)
hold on
plot(input,'r')
prettify(gcf)
xlabel('Time (ms)')
ylabel('Rate (AU)')

saveas(h,fdir + "/Figure_7_A_upper" + Oi + ".emf")
end

%% A: Learning via MFs

h = figure
hold on
opt.SE = 0;
opt.scatter = 1;
prettyWormPlot(mmf,O,opt,[0 0 1])
prettyWormPlot(mgc(abs(a)<.05),O(abs(a)<.05),opt,[1 0 0])
xlabel('Phase Offset')
ylabel('MSE')

saveas(h,fdir + "/Figure_7_A.emf")


