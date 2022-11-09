%% HEADER

% B008all_Make_Figure_8.m
% From: An emergent temporal basis set robustly supports cerebellar time-series learning
% Authors: Jesse I. Gilmer, Michael A. Farries, Zachary Kilpatrick, Ioannis Delis, Abigail L. Person
% 2022, BiorXiv, https://www.biorxiv.org/content/10.1101/2022.01.06.475265v1
%
% This script was written by JG.
% Data used in this script comes from https://zenodo.org/record/5140528#.YuA103bMJD9
%
% INPUTS: N/A.
%
% OUTPUTS: 3 figure panels to /Figure 8/


%%
clear all;
close all;

%% Make Figure 8

fdir = 'Figure 8';
mkdir(fdir);

clearvars -except fdir
load("EMG_Tests_compiled_Aug.mat")

for i = 1:length(mRes)
    nGC(i) = mRes(i).nGC;
    a(i) = mRes(i).alpha;
    thresh(i) = mRes(i).alpha;
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
        %         O(i) = mRes(i).step_prop;
        id(i) = mRes(i).TrialNumber;
    end
end

xout = [];
nGC(xout) = [];
a(xout) = [];
m(xout) = [];
mgc(xout) = [];
mmf(xout) = [];
ts(xout) = [];
tss(xout) = [];
dim(xout) = [];
mfdim(xout) = [];
cc(xout) = [];
MFcc(xout) = [];
tau(xout) = [];
mtau(xout) = [];
v(xout) = [];
lo(xout) = [];
gTau(xout) = [];
mTau(xout) = [];
ploss(xout) = [];

%%

h = figure
hold on
opt.SE = 0;
opt.scatter = 1;
prettyWormPlot(mmf,a,opt,[0 0 1])
prettyWormPlot(mgc,a,opt,[1 0 0])
xlabel('Threshold')
ylabel('MSE')
prettify(gcf)

saveas(h,fdir + "/Figure_8_A.emf")

%%
Exps = normalize([dim; ts; PCAind; v; cc; tau; lo;ploss; exvar]')';
Names = ["Dim.", "STS", "PCs", "Pop. Var.", "Corr.", "T. Sparse.",  ...
    "T. Loss", "P. Loss",  "Var. Ret."]
Names2 = [Names];

nBins = 8;
L = floor(length(mmf)/nBins);
O = [];
[S_mmf,qk] = sort(mmf,'descend');
qk(qk>nBins*L) = [];
M = [];
for i = 1:nBins
    Oi = S_mmf((i-1)*L+1:i*L);
    M(i) = mean(Oi);
end
O = [];
for i = 1:length(mmf)
    [~,mn] = min(abs(mmf(i) - M));
    O(i) = M(mn);
end

u_o = unique(O);
P = [];
P_o = [];

hold_B = [];
R_holder = [];
model_beta = [];
for i = 1:length(u_o)
    F1 = find(O == u_o(i));
    u_a = unique(a);

    Lambda = logspace(-4,2,100);
    cond = F1;
    X = Exps(:,cond)';
    % y = normalize(mgc(cond)','range');
    y = mgc(cond)';
    ba = a(cond);

    X = X';
    X(isnan(X)) = 0;
    CVMdl = fitrlinear(X,y,'ObservationsIn','columns','KFold',10,'Lambda',Lambda,...
        'Learner','leastsquares','Solver','lbfgs');
    Mdl = fitrlinear(X,y,'ObservationsIn','columns','Lambda',Lambda,...
        'Learner','leastsquares','Solver','lbfgs');

    numCLModels = numel(CVMdl.Trained)
    mse = kfoldLoss(CVMdl)
    mse_mdl = kfoldLoss(CVMdl,'mode','individual');
    mse_cat = repmat([1:length(Lambda)],numCLModels,1);


    [~,mloc] = min(mean(mse_mdl));

    [fix,fax] = min(abs(mean(mse_mdl)-(mean(mse_mdl(mse_cat==mloc))+(std(mse_mdl(mse_cat==mloc))/sqrt(10)))))


    idxFinal = fax-1
    MdlFinal = selectModels(Mdl,idxFinal)

    numNZCoeff = sum(Mdl.Beta~=0);

    Kval = 10;
    [CVMdl] = fitrlinear(X,y,'ObservationsIn','columns','KFold',Kval,'Lambda',Lambda(fax-1),...
        'Learner','leastsquares','Solver','lbfgs');
    B = [];
    for Q = 1:Kval
        B(Q,:) = CVMdl.Trained{Q}.Beta;
    end

    model_beta(i).Betas = B;
    mdl = fitlm(X'*MdlFinal.Beta, y);
    hold_B(i,:) = MdlFinal.Beta;

    C = mdl.Coefficients.Estimate;
    R_holder(i) = mdl.Rsquared.Adjusted;
end


h = figure();
hold on
p = turbo(9);
for checker = [3 9]


    Bs = [];
    As = [];
    Bx = [];

    for i = 1:length(u_o)
        B = model_beta(i).Betas;
        %     scatter(B(:,checker).*0+u_o(i),B(:,checker))
        Bs = [Bs B(:,checker)'];
        As = [As B(:,checker)'.*0+u_o(i)];
        Bx = [Bx [abs(B(:,checker))./sum(abs(B)')']'];
    end


    opt.SE = 0;
    opt.scatter = 1;
    prettyWormPlot(Bx,As,opt,p(checker,:))
end
L = [];
for i = [3 9]
    L = [L  Names(i)  "" "" ];
end
legend(L)
xlabel('MF MSE Alone')
ylabel("Coefficient Percentages")
saveas(h,fdir + "/Figure_8_B.emf")


h = figure();
hold on
p = turbo(9);
for checker = [3 9]

    Bs = [];
    As = [];
    Bx = [];

    for i = 1:length(u_o)
        B = model_beta(i).Betas;
        %     scatter(B(:,checker).*0+u_o(i),B(:,checker))
        Bs = [Bs B(:,checker)'];
        As = [As B(:,checker)'.*0+u_o(i)];
        Bx = [Bx [abs(B(:,checker))./sum(abs(B)')']'];
    end


    opt.SE = 0;
    opt.scatter = 1;
    prettyWormPlot(Bx,As,opt,p(checker,:))
end
L = [];
for i =  [3 9]
    L = [L  Names(i)  "" "" ];
end
legend(L)
xlabel('MF MSE Alone')
ylabel("Coefficient Percentages")
% saveas(h,fdir + "/Figure_8_B.emf")

h = figure();
hold on
plot(u_o,R_holder,'k-o')
ylim([0 1])
xlabel('MF MSE Alone')
ylabel("R^2 of fit")
prettify(h)
saveas(h,fdir + "/Figure_8_MFr2.emf")
saveas(h,fdir + "/Figure_8_MFr2.png")
%%

clear all
close all

fdir = 'Figure 8';
mkdir(fdir);

floc = "Datasets/D002_EMG_Data.mat";

load(floc)
kin = data.KIN;
emg = data.EMG;
kin_t = normalize(kin(:,1)','range');

w = gausswin(100);
w = w/sum(w);
a = emg;
a =  padarray(a,[1000 0],'replicate','both');
a = filtfilt(w,1,abs(a));
a = a(1001:end-1000,:);
a = normalize(a,'range');
emg = [a];


% Network params:
nMF = size(emg,2);   % number of mossy fibers. Default = 50.
nGC = 3000; % number of granule cells. Default = 3000.
conv = 4;   % convergence ratio, mf:gc. Default = 4.

% MF params:
mu = 4;    % Mean rate. Default = 4, but this value is normalized.
sigma = 1; % Rate standard deviation. Default = 1, but this value is normalized.
T = 100;   % OU Tau, the decay rate. Default = 100.

% Integrate into the options storage dictionary.
options.mu = mu;
options.sigma = sigma;

% Set the learn rate and trials.
options.lrnRate = 1e-5; %Learning rate. This will change below.
options.nTrials = 1000; %Number of trials to learn over. Default = 1000.


% Set the number of independent simulations to run:
ntr = 50; % Default = 50.

alpha = 0;

input = emg';
nPts = size(input,2);

target = kin_t;

W1 = zeros(nGC, nMF);
for k = 1:nGC
    sel = randsample(nMF, conv);
    W1(k,sel) = 1/conv;
end

mu = mean(mean(input));
sigma = mean(std(input,[],2));
thr = mu + alpha * sigma;

gc = W1 * input;
gc = gc - thr;
gc(gc < 0) = 0;

options.lrnRate = 1e-4;
[results] = RunFilter(input,target,gc,options);
WGC = results.W;


options.lrnRate = 1e-6;
[results] = RunFilter(input,target,input,options);
WMF = results.W;

h = figure();
hold on
plot(kin_t,'linewidth',5,'color', [0 0 0 .5])
plot(WGC*gc,'r','linewidth',2)
plot(WMF*input,'b','linewidth',2)
set(gcf,'color','w')
xlabel('time (ms)')
ylabel('Kinematic (with learned fits)')
prettify(gcf)
legend({'Target Kinematic', 'GCL','MF'})
axis tight
saveas(h,fdir + "/Figure_8_A_Example.png")
saveas(h,fdir + "/Figure_8_A_Example.emf")
%%
Exps = normalize([PCAind; exvar; v]')';
Names = ["PCs", "Var. Ret.",'variance']
Names2 = [Names];

nBins = 10;
L = floor(length(mmf)/nBins);
O = [];
[S_mmf,qk] = sort(mmf,'descend');
qk(qk>nBins*L) = [];
M = [];
for i = 1:nBins
    Oi = S_mmf((i-1)*L+1:i*L);
    M(i) = mean(Oi);
end
O = [];
for i = 1:length(mmf)
    [~,mn] = min(abs(mmf(i) - M));
    O(i) = M(mn);
end

u_o = unique(O);
P = [];
P_o = [];

hold_B = [];
R_holder = [];
model_beta = [];
for i = 1:length(u_o)
    F1 = find(O == u_o(i));
    u_a = unique(a);

    Lambda = logspace(-4,2,100);
    cond = F1;
    X = Exps(:,cond)';
    % y = normalize(mgc(cond)','range');
    y = mgc(cond)';
    ba = a(cond);

    X = X';
    X(isnan(X)) = 0;
    CVMdl = fitrlinear(X,y,'ObservationsIn','columns','KFold',10,'Lambda',Lambda,...
        'Learner','leastsquares','Solver','lbfgs');
    Mdl = fitrlinear(X,y,'ObservationsIn','columns','Lambda',Lambda,...
        'Learner','leastsquares','Solver','lbfgs');

    numCLModels = numel(CVMdl.Trained)
    mse = kfoldLoss(CVMdl)
    mse_mdl = kfoldLoss(CVMdl,'mode','individual');
    mse_cat = repmat([1:length(Lambda)],numCLModels,1);


    [~,mloc] = min(mean(mse_mdl));

    [fix,fax] = min(abs(mean(mse_mdl)-(mean(mse_mdl(mse_cat==mloc))+(std(mse_mdl(mse_cat==mloc))/sqrt(10)))))


    idxFinal = fax-1
    MdlFinal = selectModels(Mdl,idxFinal)

    numNZCoeff = sum(Mdl.Beta~=0);

    Kval = 10;
    [CVMdl] = fitrlinear(X,y,'ObservationsIn','columns','KFold',Kval,'Lambda',Lambda(fax-1),...
        'Learner','leastsquares','Solver','lbfgs');
    B = [];
    for Q = 1:Kval
        B(Q,:) = CVMdl.Trained{Q}.Beta;
    end

    model_beta(i).Betas = B;
    mdl = fitlm(X'*MdlFinal.Beta, y);
    hold_B(i,:) = MdlFinal.Beta;

    C = mdl.Coefficients.Estimate;
    R_holder(i) = mdl.Rsquared.Adjusted;
end


h = figure();
hold on
p = turbo(3);
for checker = [1 2 3]


    Bs = [];
    As = [];
    Bx = [];

    for i = 1:length(u_o)
        B = model_beta(i).Betas;
        %     scatter(B(:,checker).*0+u_o(i),B(:,checker))
        Bs = [Bs B(:,checker)'];
        As = [As B(:,checker)'.*0+u_o(i)];
        Bx = [Bx [abs(B(:,checker))./sum(abs(B)')']'];
    end


    opt.SE = 0;
    opt.scatter = 1;
    prettyWormPlot(Bx,As,opt,p(checker,:))
end
L = [];
for i = [1 2 3 ]
    L = [L  Names(i)  "" "" ];
end
legend(L)
xlabel('MF MSE Alone')
ylabel("Coefficient Percentages")
% saveas(h,fdir + "/Figure_8_B.emf")


h = figure();
hold on
plot(u_o,R_holder,'k-o')
ylim([0 1])
xlabel('MF MSE Alone')
ylabel("R^2 of fit")
% prettify(h)
% saveas(h,fdir + "/Figure_8_MFr2.emf")
% saveas(h,fdir + "/Figure_8_MFr2.png")

%%

Exps = normalize([PCAind]')';
Names = ["PCs"]
Names2 = [Names];

nBins = 8;
L = floor(length(mmf)/nBins);
O = [];
[S_mmf,qk] = sort(mmf,'descend');
qk(qk>nBins*L) = [];
M = [];
for i = 1:nBins
    Oi = S_mmf((i-1)*L+1:i*L);
    M(i) = mean(Oi);
end
O = [];
for i = 1:length(mmf)
    [~,mn] = min(abs(mmf(i) - M));
    O(i) = M(mn);
end

u_o = unique(O);
P = [];
P_o = [];

hold_B = [];
R_holder = [];
model_beta = [];
for i = 1:length(u_o)
    F1 = find(O == u_o(i));
    u_a = unique(a);

    Lambda = logspace(-4,2,100);
    cond = F1;
    X = Exps(:,cond)';
    % y = normalize(mgc(cond)','range');
    y = mgc(cond)';
    ba = a(cond);

    X = X';
    X(isnan(X)) = 0;
    CVMdl = fitrlinear(X,y,'ObservationsIn','columns','KFold',10,'Lambda',Lambda,...
        'Learner','leastsquares','Solver','lbfgs');
    Mdl = fitrlinear(X,y,'ObservationsIn','columns','Lambda',Lambda,...
        'Learner','leastsquares','Solver','lbfgs');

    numCLModels = numel(CVMdl.Trained)
    mse = kfoldLoss(CVMdl)
    mse_mdl = kfoldLoss(CVMdl,'mode','individual');
    mse_cat = repmat([1:length(Lambda)],numCLModels,1);


    [~,mloc] = min(mean(mse_mdl));

    [fix,fax] = min(abs(mean(mse_mdl)-(mean(mse_mdl(mse_cat==mloc))+(std(mse_mdl(mse_cat==mloc))/sqrt(10)))))


    idxFinal = fax-1
    MdlFinal = selectModels(Mdl,idxFinal)

    numNZCoeff = sum(Mdl.Beta~=0);

    Kval = 10;
    [CVMdl] = fitrlinear(X,y,'ObservationsIn','columns','KFold',Kval,'Lambda',Lambda(fax-1),...
        'Learner','leastsquares','Solver','lbfgs');
    B = [];
    for Q = 1:Kval
        B(Q,:) = CVMdl.Trained{Q}.Beta;
    end

    model_beta(i).Betas = B;
    mdl = fitlm(X'*MdlFinal.Beta, y);
    hold_B(i,:) = MdlFinal.Beta;

    C = mdl.Coefficients.Estimate;
    R_holder(i) = mdl.Rsquared.Adjusted;
end

h = figure();
hold on
plot(u_o,R_holder,'k-o')
ylim([0 1])
xlabel('MF MSE Alone')
ylabel("R^2 of fit")
prettify(h)
saveas(h,fdir + "/Figure_8_PCALONE.emf")

%%

Exps = normalize([exvar]')';
Names = ["Var. Ret."]
Names2 = [Names];

nBins = 8;
L = floor(length(mmf)/nBins);
O = [];
[S_mmf,qk] = sort(mmf,'descend');
qk(qk>nBins*L) = [];
M = [];
for i = 1:nBins
    Oi = S_mmf((i-1)*L+1:i*L);
    M(i) = mean(Oi);
end
O = [];
for i = 1:length(mmf)
    [~,mn] = min(abs(mmf(i) - M));
    O(i) = M(mn);
end

u_o = unique(O);
P = [];
P_o = [];

hold_B = [];
R_holder = [];
model_beta = [];
for i = 1:length(u_o)
    F1 = find(O == u_o(i));
    u_a = unique(a);

    Lambda = logspace(-4,2,100);
    cond = F1;
    X = Exps(:,cond)';
    % y = normalize(mgc(cond)','range');
    y = mgc(cond)';
    ba = a(cond);

    X = X';
    X(isnan(X)) = 0;
    CVMdl = fitrlinear(X,y,'ObservationsIn','columns','KFold',10,'Lambda',Lambda,...
        'Learner','leastsquares','Solver','lbfgs');
    Mdl = fitrlinear(X,y,'ObservationsIn','columns','Lambda',Lambda,...
        'Learner','leastsquares','Solver','lbfgs');

    numCLModels = numel(CVMdl.Trained)
    mse = kfoldLoss(CVMdl)
    mse_mdl = kfoldLoss(CVMdl,'mode','individual');
    mse_cat = repmat([1:length(Lambda)],numCLModels,1);


    [~,mloc] = min(mean(mse_mdl));

    [fix,fax] = min(abs(mean(mse_mdl)-(mean(mse_mdl(mse_cat==mloc))+(std(mse_mdl(mse_cat==mloc))/sqrt(10)))))


    idxFinal = fax-1
    MdlFinal = selectModels(Mdl,idxFinal)

    numNZCoeff = sum(Mdl.Beta~=0);

    Kval = 10;
    [CVMdl] = fitrlinear(X,y,'ObservationsIn','columns','KFold',Kval,'Lambda',Lambda(fax-1),...
        'Learner','leastsquares','Solver','lbfgs');
    B = [];
    for Q = 1:Kval
        B(Q,:) = CVMdl.Trained{Q}.Beta;
    end

    model_beta(i).Betas = B;
    mdl = fitlm(X'*MdlFinal.Beta, y);
    hold_B(i,:) = MdlFinal.Beta;

    C = mdl.Coefficients.Estimate;
    R_holder(i) = mdl.Rsquared.Adjusted;
end

h = figure();
hold on
plot(u_o,R_holder,'k-o')
ylim([0 1])
xlabel('MF MSE Alone')
ylabel("R^2 of fit")
prettify(h)
saveas(h,fdir + "/Figure_8_RETLONE.emf")

%%

Exps = normalize([PCAind; exvar]')';
Names = ["PCs","Var. Ret."]
Names2 = [Names];

nBins = 8;
L = floor(length(mmf)/nBins);
O = [];
[S_mmf,qk] = sort(mmf,'descend');
qk(qk>nBins*L) = [];
M = [];
for i = 1:nBins
    Oi = S_mmf((i-1)*L+1:i*L);
    M(i) = mean(Oi);
end
O = [];
for i = 1:length(mmf)
    [~,mn] = min(abs(mmf(i) - M));
    O(i) = M(mn);
end

u_o = unique(O);
P = [];
P_o = [];

hold_B = [];
R_holder = [];
model_beta = [];
for i = 1:length(u_o)
    F1 = find(O == u_o(i));
    u_a = unique(a);

    Lambda = logspace(-4,2,100);
    cond = F1;
    X = Exps(:,cond)';
    % y = normalize(mgc(cond)','range');
    y = mgc(cond)';
    ba = a(cond);

    X = X';
    X(isnan(X)) = 0;
    CVMdl = fitrlinear(X,y,'ObservationsIn','columns','KFold',10,'Lambda',Lambda,...
        'Learner','leastsquares','Solver','lbfgs');
    Mdl = fitrlinear(X,y,'ObservationsIn','columns','Lambda',Lambda,...
        'Learner','leastsquares','Solver','lbfgs');

    numCLModels = numel(CVMdl.Trained)
    mse = kfoldLoss(CVMdl)
    mse_mdl = kfoldLoss(CVMdl,'mode','individual');
    mse_cat = repmat([1:length(Lambda)],numCLModels,1);


    [~,mloc] = min(mean(mse_mdl));

    [fix,fax] = min(abs(mean(mse_mdl)-(mean(mse_mdl(mse_cat==mloc))+(std(mse_mdl(mse_cat==mloc))/sqrt(10)))))


    idxFinal = fax-1
    MdlFinal = selectModels(Mdl,idxFinal)

    numNZCoeff = sum(Mdl.Beta~=0);

    Kval = 10;
    [CVMdl] = fitrlinear(X,y,'ObservationsIn','columns','KFold',Kval,'Lambda',Lambda(fax-1),...
        'Learner','leastsquares','Solver','lbfgs');
    B = [];
    for Q = 1:Kval
        B(Q,:) = CVMdl.Trained{Q}.Beta;
    end

    model_beta(i).Betas = B;
    mdl = fitlm(X'*MdlFinal.Beta, y);
    hold_B(i,:) = MdlFinal.Beta;

    C = mdl.Coefficients.Estimate;
    R_holder(i) = mdl.Rsquared.Adjusted;
end

h = figure();
hold on
plot(u_o,R_holder,'k-o')
ylim([0 1])
xlabel('MF MSE Alone')
ylabel("R^2 of fit")
prettify(h)
saveas(h,fdir + "/Figure_8_both.emf")

%%

Exps = normalize([PCAind; v; exvar]')';
Names = ["PCs", "Pop. Var.","Var. Ret."]
Names2 = [Names];

nBins = 8;
L = floor(length(mmf)/nBins);
O = [];
[S_mmf,qk] = sort(mmf,'descend');
qk(qk>nBins*L) = [];
M = [];
for i = 1:nBins
    Oi = S_mmf((i-1)*L+1:i*L);
    M(i) = mean(Oi);
end
O = [];
for i = 1:length(mmf)
    [~,mn] = min(abs(mmf(i) - M));
    O(i) = M(mn);
end

u_o = unique(O);
P = [];
P_o = [];

hold_B = [];
R_holder = [];
model_beta = [];
for i = 1:length(u_o)
    F1 = find(O == u_o(i));
    u_a = unique(a);

    Lambda = logspace(-4,2,100);
    cond = F1;
    X = Exps(:,cond)';
    % y = normalize(mgc(cond)','range');
    y = mgc(cond)';
    ba = a(cond);

    X = X';
    X(isnan(X)) = 0;
    CVMdl = fitrlinear(X,y,'ObservationsIn','columns','KFold',10,'Lambda',Lambda,...
        'Learner','leastsquares','Solver','lbfgs');
    Mdl = fitrlinear(X,y,'ObservationsIn','columns','Lambda',Lambda,...
        'Learner','leastsquares','Solver','lbfgs');

    numCLModels = numel(CVMdl.Trained)
    mse = kfoldLoss(CVMdl)
    mse_mdl = kfoldLoss(CVMdl,'mode','individual');
    mse_cat = repmat([1:length(Lambda)],numCLModels,1);


    [~,mloc] = min(mean(mse_mdl));

    [fix,fax] = min(abs(mean(mse_mdl)-(mean(mse_mdl(mse_cat==mloc))+(std(mse_mdl(mse_cat==mloc))/sqrt(10)))))


    idxFinal = fax-1
    MdlFinal = selectModels(Mdl,idxFinal)

    numNZCoeff = sum(Mdl.Beta~=0);

    Kval = 10;
    [CVMdl] = fitrlinear(X,y,'ObservationsIn','columns','KFold',Kval,'Lambda',Lambda(fax-1),...
        'Learner','leastsquares','Solver','lbfgs');
    B = [];
    for Q = 1:Kval
        B(Q,:) = CVMdl.Trained{Q}.Beta;
    end

    model_beta(i).Betas = B;
    mdl = fitlm(X'*MdlFinal.Beta, y);
    hold_B(i,:) = MdlFinal.Beta;

    C = mdl.Coefficients.Estimate;
    R_holder(i) = mdl.Rsquared.Adjusted;
end

h = figure();
hold on
plot(u_o,R_holder,'k-o')
ylim([0 1])
xlabel('MF MSE Alone')
ylabel("R^2 of fit")
prettify(h)
saveas(h,fdir + "/Figure_8_bothwVar.emf")


h = figure();
hold on
p = turbo(3);
for checker = [1 2 3]

    Bs = [];
    As = [];
    Bx = [];

    for i = 1:length(u_o)
        B = model_beta(i).Betas;
        %     scatter(B(:,checker).*0+u_o(i),B(:,checker))
        Bs = [Bs B(:,checker)'];
        As = [As B(:,checker)'.*0+u_o(i)];
        Bx = [Bx [abs(B(:,checker))./sum(abs(B)')']'];
    end


    opt.SE = 0;
    opt.scatter = 1;
    prettyWormPlot(Bx,As,opt,p(checker,:))
end
L = [];
for i = [1 2 3 ]
    L = [L  Names(i)  "" "" ];
end
legend(L)
xlabel('MF MSE Alone')
ylabel("Coefficient Percentages")
% saveas(h,fdir + "/Figure_8_B.emf")
