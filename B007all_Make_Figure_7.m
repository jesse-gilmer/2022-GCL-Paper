%% HEADER

% B007all_Make_Figure_7.m
% From: An emergent temporal basis set robustly supports cerebellar time-series learning
% Authors: Jesse I. Gilmer, Michael A. Farries, Zachary Kilpatrick, Ioannis Delis, Abigail L. Person
% 2022, BiorXiv, https://www.biorxiv.org/content/10.1101/2022.01.06.475265v1
%
% This script was written by JG.
% Data used in this script comes from https://zenodo.org/record/5140528#.YuA103bMJD9
%
% INPUTS: N/A.
%
% OUTPUTS: 3 figure panels to /Figure 7/

%% Make Figure 7
close all;
clear all;

fdir = "Figure 7"
mkdir(fdir)

%% Preprocessing
clearvars -except fdir
load("OU_Datasets_compiled.mat")

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

% m = normalize(m,'range');

Exps = normalize([dim; ts; PCAind; v; cc; tau; lo; ploss; exvar]','zscore')';
Names = ["Dim.", "STS", "PCs", "Pop. Var.", "Corr.", "T. Sparse.",  ...
    "T. Loss", "P. Loss",  "Var. Ret."]

%% Make Figure 6 Panel A: LASSO Description

Lambda = logspace(-4,1,100);
cond = a < 4;
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
% prettify(gcf);
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

% Make Figure 6 Panel B: Model unity plot.


idxFinal = fax-1;
MdlFinal = selectModels(Mdl,idxFinal);

xst = regstats(y,X'*MdlFinal.Beta);

h = figure;
sgtitle("Number of coefs: " + numNZCoeff(idxFinal))
mdl = fitlm(X'*MdlFinal.Beta, y)
plot(mdl)
% prettify(gcf)
title("R2: " + mdl.Rsquared.Adjusted)
xlabel('MSE')
ylabel('Model')

saveas(h,fdir + "/Figure_6B.emf")
saveas(h,fdir + "/Figure_6B.png")

% Make Figure 6 Panel C: Beta Components

h = figure();
bar(MdlFinal.Beta)
xticks(1:length(Names))
xticklabels([Names])
% prettify(gcf)
grid on

saveas(h,fdir + "/Figure_6B.emf")
saveas(h,fdir + "/Figure_6B.png")

% Make Figure 6 Panel D: Model comparison

h = figure();
C = mdl.Coefficients.Estimate;
opt.SE = 0;
opt.scatter = 1;
prettyWormPlot(mdl.Fitted,ba,opt,[1 .25 .5])
prettyWormPlot(y',ba,opt,[.5 .25 1])
ylim([-.2 5])
xlabel('Threshold')
ylabel('MSE')
title("MSE: " + mdl.RMSE)

saveas(h,fdir + "/Figure_6D.emf")
saveas(h,fdir + "/Figure_6D.png")


%%

Exps = normalize([dim; ts; PCAind; v; cc; tau; lo;ploss; exvar]')';
% [zz,ii] = find(a > -1 & a < 1)
% Exps = Exps(:,ii);

Names = ["Dim.", "STS", "PCs", "Pop. Var.", "Corr.", "T. Sparse.",  ...
    "T. Loss", "P. Loss",  "Var. Ret."] ;

cond = a < 4;
sel = 1:9;
sel = [2 3 4 7 8 9];
X = [Exps(sel,cond)']
y = m(cond)';
mdl = stepwiselm(X,y, 'Intercept',false,...
   'Criterion','BIC','PredictorVars',Names(sel))

ba = a(cond);

%%
Names = (strrep(mdl2.CoefficientNames,":" ," x "));
Names = Names(2:end);

Lambda = logspace(-4,1,100);
cond = a < 4;
Lambda = logspace(-4,1,100);
% cond = (a > -2)&(a < 1);
X = Exps(:,cond)';
y = m(cond)';
ba = a(cond);
sel = [2 3 4 8 9];
X = [(X(:,sel)) (X(:,3).*X(:,4)) (X(:,3).*X(:,8)) (X(:,4).*X(:,8)) (X(:,4).*X(:,9))];
y = m(cond)';
ba = a(cond);

X = X'; 

CVMdl = fitrlinear(X,y,'ObservationsIn','columns','KFold',10,'Lambda',Lambda,...
    'Learner','leastsquares','Solver','lbfgs');
[Mdl,FitInfo] = fitrlinear(X,y,'ObservationsIn','columns','Lambda',Lambda,...
    'Learner','leastsquares','Solver','lbfgs');

numCLModels = numel(CVMdl.Trained);
mse = kfoldLoss(CVMdl);


mse_mdl = kfoldLoss(CVMdl,'mode','individual');
mse_cat = repmat([1:length(Lambda)],numCLModels,1);

[fix,fax] = min(abs(mean(mse_mdl)-(mean(mse_mdl(mse_cat==1))+std(mse_mdl(mse_cat==1)))));

numNZCoeff = sum(Mdl.Beta~=0);

idxFinal = fax-1;
idxFinal = 1;
MdlFinal = selectModels(Mdl,idxFinal);

xst = regstats(y,X'*MdlFinal.Beta);

h = figure();
C = mdl.Coefficients.Estimate;
opt.SE = 0;
opt.scatter = 1;
prettyWormPlot(mdl.Fitted,ba,opt,[1 .25 .5])
prettyWormPlot(y',ba,opt,[.5 .25 1])
xlabel('Threshold')
ylabel('MSE')
title("MSE: " + mdl.RMSE)
ylim([-.2 5])
saveas(h,fdir + "/Figure_6E.emf")
saveas(h,fdir + "/Figure_6E.png")

h = figure();
bar(MdlFinal.Beta)
xticks(1:length(Names))
xticklabels([Names])
% prettify(gcf)
grid on

saveas(h,fdir + "/Figure_6F.emf")
saveas(h,fdir + "/Figure_6F.png")

