
clearvars -except mRes
close all

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
        O(i) = mRes(i).offset;
    end
end

%%
close all


figure
hold on
opt.SE = 0;
opt.scatter = 1;
prettyWormPlot(mmf,O./pi,opt,[0 0 0])
xlabel('Phase offset')
ylabel('MF MSE')


figure
hold on
u_o = unique(O)
for i = 1:length(u_o);
        subplot(5,4,i)
F = find(O == u_o(i));


opt.SE = 0;
opt.scatter = 1;
prettyWormPlot(mgc(F),a(F),opt,[1/i log(i)/5 1-(1/i)])
prettyWormPlot([mmf(F) mmf(F)],[ones(1,length(mmf(F)))*min(a(F)) ones(1,length(mmf(F)))*max(a(F))],opt,[0 0 0])
xlabel('Offset')
ylabel('MSE')
prettify(gcf)
axis square
ylim([0 .25])
title(u_o(i)/(2*pi))
end


figure
hold on
u_o = unique(O);
P = [];
P_o = [];
for i = 1:length(u_o);
F1 = find(O == u_o(i));
u_a = unique(a);
p = [];
for j = 1:length(u_a)
    F2 = find(a == u_a(j));
    S = intersect(F1,F2);
    p(j) = mean(m(S));
    pv(j) = mean(exvar(S));
end
[P(i),qq] = min(p);
P_o(i) = u_a(qq);
P_v(i) = pv(qq);

end

figure();
hold on
plot(u_o/pi,P_o,'k-o')
prettify(gcf);
ylabel('Best Threshold')
xlabel('Phase offset')
axis square

figure
plot(u_o/pi,P_v,'r-o')
prettify(gcf);
ylabel('Variance Retained')
xlabel('Phase offset')
axis square

figure()
plot(u_o/pi,P,'r-o')
xlabel('Phase Offset')
ylabel('MSE Ratio')
prettify(gcf)
axis square

ylim([0 .25])
title(u_o(i)/(2*pi))


figure
hold on
u_o = unique(a)
for i = 1:length(u_o);
        subplot(5,4,i)
F = find(a == u_o(i));


opt.SE = 0;
opt.scatter = 1;
prettyWormPlot(mgc(F),O(F),opt,[1/i log(i)/5 1-(1/i)])
xlabel('Offset')
ylabel('MSE')
prettify(gcf)
axis square
ylim([0 .5])
end


figure
hold on
u_o = unique(O)
for i = 1:length(u_o);
    subplot(5,4,i)
    hold on
F = find(O == u_o(i));

opt.SE = 0;
opt.scatter = 1;
% scatter(exvar(F),log10(m(F)),5,[1/i log(i)/5 1-(1/i)])
mdl = fitlm(exvar(F),m(F))
plot(mdl)
xlabel('Ex Var')
ylabel('MSE')
prettify(gcf)
axis square
title(u_o(i)+ ", R2:" + mdl.Rsquared.Adjusted)
legend off
R2(i) = mdl.Rsquared.Adjusted;
end

figure()
plot(u_o/pi,R2,'r-o')
xlabel('Phase Offset')
ylabel('R2(Variance Retained, MSE)')
prettify(gcf)
axis square
