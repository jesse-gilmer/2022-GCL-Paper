%% Make Figure 6, Updated, Again.
close all;
clear all;

fdir = "Figure 6 Panels"
mkdir(fdir)

%% Preprocessing
clearvars -except fdir
load("PCA_test_3000F_retest_compiled_Aug.mat")

%% Sort data
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

m = normalize(m,'range');

Exps = normalize([dim; ts; PCAind; v; cc; tau; lo;ploss; exvar]','zscore')';
Names = ["Dim.", "STS", "PCs", "Pop. Var.", "Corr.", "T. Sparse.",  ...
    "T. Loss", "P. Loss",  "Var. Ret."]

%% Make Figure 6 Panel A: LASSO Description

Lambda = logspace(-4,1,100);
cond = m>=0;
X = Exps(:,cond)';
y = m(cond)';
ba = a(cond);

X = X'; 

CVMdl = fitrlinear(X,y,'ObservationsIn','columns','KFold',10,'Lambda',Lambda,...
    'Learner','leastsquares','Solver','sparsa');
[Mdl,FitInfo] = fitrlinear(X,y,'ObservationsIn','columns','Lambda',Lambda,...
    'Learner','leastsquares','Solver','sparsa');

numCLModels = numel(CVMdl.Trained);
mse = kfoldLoss(CVMdl);

mse_mdl = kfoldLoss(CVMdl,'mode','individual');
mse_cat = repmat([1:length(Lambda)],numCLModels,1);

[fix,fax] = min(abs(mean(mse_mdl)-(mean(mse_mdl(mse_cat==1))+std(mse_mdl(mse_cat==1)))));

numNZCoeff = sum(Mdl.Beta~=0);

h = figure();
hold on

opt.SE = 0;
opt.scatter = 1;
prettyWormPlot(mse_mdl(:)',log10(Lambda(mse_cat(:))),opt,[.5 .25 1])
plot(log10([min(Lambda) max(Lambda)]),[mean(mse_mdl(mse_cat==1))+std(mse_mdl(mse_cat==1)) mean(mse_mdl(mse_cat==1))+std(mse_mdl(mse_cat==1))],'k--')
scatter(log10(Lambda(fax-1)),mean(mse_mdl(mse_cat==fax-1)),'k','filled')
prettify(gcf);
ylabel('MSE of LASSO')

yyaxis right
plot(log10(Lambda),numNZCoeff,'r-o'); 
scatter(log10(Lambda(fax-1)),numNZCoeff(fax-1),'r','filled')
xlabel('log_{10} LASSO Lambda');
ylabel('Number of components')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
saveas(h,fdir + "/Figure_6A.emf")
saveas(h,fdir + "/Figure_6A.png")

%% Make Figure 6 Panel B: Model unity plot.


idxFinal = fax-1;
MdlFinal = selectModels(Mdl,idxFinal);

xst = regstats(y,X'*MdlFinal.Beta);

h = figure;
sgtitle("Number of coefs: " + numNZCoeff(idxFinal))
mdl = fitlm(X'*MdlFinal.Beta, y)
plot(mdl)
prettify(gcf)
title("R2: " + mdl.Rsquared.Adjusted)
xlabel('MSE')
ylabel('Model')

saveas(h,fdir + "/Figure_6B.emf")
saveas(h,fdir + "/Figure_6B.png")

%% Make Figure 6 Panel C: Beta Components

h = figure();
bar(MdlFinal.Beta)
xticks(1:length(Names))
xticklabels([Names])
prettify(gcf)
grid on

saveas(h,fdir + "/Figure_6B.emf")
saveas(h,fdir + "/Figure_6B.png")

%% Make Figure 6 Panel D: Model comparison

h = figure();
C = mdl.Coefficients.Estimate;
opt.SE = 0;
opt.scatter = 1;
prettyWormPlot(mdl.Fitted,ba,opt,[1 .25 .5])
prettyWormPlot(y',ba,opt,[.5 .25 1])
xlabel('Threshold')
ylabel('MSE')
title("MSE: " + mdl.RMSE)

saveas(h,fdir + "/Figure_6D.emf")
saveas(h,fdir + "/Figure_6D.png")
