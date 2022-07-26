%% Make Figure 2, updated:
clear all
close all

fdir = "Figure 2 Updated"
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

%% Figure 2A:
h = figure();
hold on
opt.SE = 0;
opt.scatter = 1;
prettyWormPlot(mgc,a,opt,[1 0 0])
prettyWormPlot([mmf mmf],[ones(1,length(mmf))*min(a) ones(1,length(mmf))*max(a)],opt,[0 0 1])
ylabel('MSE')
xlabel('Threshold')
prettify(gcf)
axis square
hold on 
legend({'','GCL','','',"MFs"})

saveas(h,fdir + "/Figure_2A.emf")

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

saveas(h,fdir + "/Figure_2A_zoom_in.emf")


%% Prep for step size

clearvars -except fdir
load("Step_Test_compiled.mat")

% Filter data
for i = 1:length(mRes)-1
    if length(mRes(i).step) > 0
        held_step = mRes(i).step;
    else
        mRes(i).step = held_step;
    end
    step(i) = mRes(i).step;
    m(i) = (mRes(i).mse/mRes(i).MF_mse);
    mgc(i) = mRes(i).mse;
    mmf(i) = mRes(i).MF_mse;
    a(i) = mRes(i).alpha;
end
%%

h=figure();
u_s = unique(step);
u_a = unique(a);
par = parula(length(u_s));
for i = 1:length(u_s)
    hold on
    ca = find(step==u_s(i));
    z = [];
    az = [];
    for j = 1:length(u_a)
        ca = intersect(find(step==u_s(i)),find(a==u_a(j)));
        z(j) = mean(mgc(ca));
        zm(j) = mean(mmf(ca));
        az(j) = mean(a(ca));
    end
    scatter(log10(u_s(i)),log10(min(z)),15,'r','filled')
    scatter(log10(u_s(i)),log10(min(zm)),15,'b','filled')
    [aa,bb] = min(zm);
    z_min(i) = u_a(bb);
    s_min(i) = u_s(i);
end
xticks(log10(u_s))
xticklabels(u_s)
yticklabels(10.^yticks)
xlabel('Update magnitude (log_1_0 scale)')
ylabel('minimized MSE (log_1_0 scale)')

prettify(gcf)
saveas(h,fdir + "/Figure_2C.emf")
saveas(h,fdir + "/Figure_2C.png")

%%

h=figure();
u_s = unique(step);
u_a = unique(a);
par = parula(length(u_s))
for i = 1:length(u_s)
    hold on
    ca = find(step==u_s(i));
    z = [];
    az = [];
    for j = 1:length(u_a)
        ca = intersect(find(step==u_s(i)),find(a==u_a(j)));
        z(j) = mean(mgc(ca));
        zm(j) = mean(mmf(ca));
        az(j) = mean(a(ca));
    end
    z(z>1) = NaN;
    zm(zm>1) = NaN;
    [aa,bb] = min(z);
    
    scatter(log10(u_s(i)),az(bb),15,'r','filled')
end
xticks(log10(u_s))
xticklabels(u_s)
xlabel('Update magnitude (log_1_0 scale)')
ylabel('Threshold z at minimized MSE')
axis padded

prettify(gcf)
saveas(h,fdir + "/Figure_2D.emf")
saveas(h,fdir + "/Figure_2D.png")


%%


m(m>10) = 10; 
mmf(mmf>1) = 1; 
mgc(mmf>1) = 1; 


figure();
u_s = unique(step);
u_a = unique(a);
par = parula(length(u_s))
for i = 1:length(u_s)
    hold on
    ca = find(step==u_s(i));
    z = [];
    az = [];
    for j = 1:length(u_a)
        ca = intersect(find(step==u_s(i)),find(a==u_a(j)));
        z(j) = mean(mgc(ca));
        zm(j) = mean(mmf(ca));
        az(j) = mean(a(ca));
    end
    z(z>1) = NaN;
    zm(zm>1) = NaN;
    plot3(az,z.*0+i,z,'color',par(i,:))
    plot3(az,z.*0+i,zm,'color','k')
end
zlim([0 .2])


figure();
u_s = unique(step);
u_a = unique(a);
par = parula(length(u_s))
for i = 1:length(u_s)
    hold on
    ca = find(step==u_s(i));
    z = [];
    az = [];
    for j = 1:length(u_a)
        ca = intersect(find(step==u_s(i)),find(a==u_a(j)));
        z(j) = mean(mmf(ca));
        az(j) = mean(a(ca));
    end
    z(z>1) = NaN;
    az(az>1) = NaN;
    plot3(az,z.*0+i,z,'color',par(i,:))
end


figure();
u_s = unique(step);
u_a = unique(a);
par = parula(length(u_s));
for i = 1:length(u_s)
    hold on
    ca = find(step==u_s(i));
    z = [];
    az = [];
    for j = 1:length(u_a)
        ca = intersect(find(step==u_s(i)),find(a==u_a(j)));
        z(j) = mean(m(ca));
        az(j) = mean(a(ca));
    end
    [a,b] = nanmin(az);
    scatter(log(u_s(i)),a,5,par(i,:))
end


figure();
hold on
u_s = unique(step);
u_a = unique(a);
for i = 1:length(u_s)
   for j = 1:length(u_a)
       ca = intersect(find(step==u_s(i)),find(a==u_a(j)));
        scatter(mean(step(ca)),mean(m(ca)),5,i/length(u_s))
   end
end




    nGC(i) = mRes(i).nGC;
    a(i) = mRes(i).alpha;
    
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

