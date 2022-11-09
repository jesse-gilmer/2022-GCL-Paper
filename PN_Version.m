%% Raster / PSTH plots centered on cue for example neurons.
close all;
clear all;
load('D001_JG_PN_Edited_Data.mat')
cells4raster = [1:128]
lims_raster = [0 500];
lims_pulse = [-50 50];
w_fr_pulse = 5;
figure();
hold on
rh = [];
ix = 1
for i = 1:50
    ii = cells4raster(i);
    [r_tmp t_tmp] = firing_rate({cell2mat(t_spk_peri_lift{ii}')},10,1,lims_raster+[-200 200],'gaussian');
    plot(t_tmp,normalize(r_tmp,'range')+ix)
    rh(i,:) = r_tmp;
    ix = ix + 1;
end

input = normalize(rh','range')';

nMF = size(input,1); %Default 50
nGC = 300; %Default 3000
nPk = 1;

conv = 4;



[sgc,xgc] = max(input');
[meow,mix] = sort(xgc);
figure(); hold on
for i = 1:nMF
plot(input(mix(i),:) - i,'k')
end
axis tight
%%
W1 = zeros(1,nMF);
for k = 1:nGC
    sel = randsample(nMF, 4);
    W1(k,sel) = 1/4;
    %     W1(k,sel) = rand(1,4);
end

mu = mean(mean(input));
sigma = mean(std(input,[],2));


alpha = .5
thr = mu + alpha * sigma;

gc = W1 * input;
gc = gc - thr;
gc(gc < 0) = 0;

[sgc,xgc] = max(gc');
[meow,mix] = sort(xgc);
zGC= (max(gc')-mean(gc'))./std(gc');

figure(); hold on
colormap turbo
imagesc(flipud(gc(mix,:)))
set(gcf,'color','w')
xlabel('time (ms)')
ylabel('GC unit')
% scatter(xgc(mix),nGC-[1:nGC],2,[1 .1 .2],'filled','MarkerFaceAlpha',1)
prettify(gcf)
colorbar()
% figure(); hold on
% imagesc((gc(mix,:)-mean(gc(mix,:)))./std(gc(mix,:)))
% caxis([-4 4])
axis tight

[sgc,xgc] = max(input');
[meow,mix] = sort(xgc);
zMF = (max(input')-mean(input'))./std(input');

figure(); hold on
colormap turbo
imagesc(flipud(input(mix,:)))
set(gcf,'color','w')
xlabel('time (ms)')
ylabel('MF unit')
% scatter(xgc(mix),nMF-[1:nMF],4,[1 .1 .2],'filled','MarkerFaceAlpha',1)
prettify(gcf)
% caxis([-3 3])
colorbar()
axis tight

% figure(); hold on
% imagesc((input(mix,:)-mean(input(mix,:)')')./std(input(mix,:)')')
% caxis([-4 4])

figure();
hold on;
boxplot([zGC zMF],[zGC.*0 zMF.*0+1])
xticks([1 2])
xticklabels(["GCL" "MFs"])
xlabel('Input type')
ylabel('Peak firing z-score')
prettify(gcf)
xlim([.5 2.5])

[sgc,xgc] = max(gc');
[meow,mix] = sort(xgc);
mlem = diff(meow);
mlem(mlem == 0) = [];

[sgc,xgc] = max(input');
[meow,mix] = sort(xgc);
mlemi = diff(meow);
mlemi(mlemi == 0) = [];

figure();
hold on;
boxplot([mlem mlemi],[mlem.*0 mlemi.*0+1],'notch','on')
xticks([1 2])
xticklabels(["GCL" "MFs"])
xlabel('Input type')
ylabel('\DeltaPeak firing (U -> U+1) (ms)')
prettify(gcf)
xlim([.5 2.5])

%%
gc(sum(gc') <.001,:) = [];

gc = normalize(gc','range')'
[sgc,xgc] = max(gc');
[meow,mix] = sort(xgc);
figure(); hold on
for i = 1:1:length(gc')
plot(gc(mix(i),:) - i,'k')
end
axis padded