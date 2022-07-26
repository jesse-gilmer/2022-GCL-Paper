clear all;
close all;

%% Make Figure 5

fdir = "Figure 5"
mkdir(fdir)
%%

h = figure
T = zeros(7,7);
for i = 1:7
    T(i,i) = .8;
end

imagesc(T)
prettify(gcf)
caxis([0 1])
colormap('bone')
saveas(h,fdir + "/Figure_5_icon_highSTS.emf")
saveas(h,fdir + "/Figure_5_icon_highSTS.png")


h = figure
T = zeros(7,7);
for i = 1:7

    T(i,max([i-2 1]):min([i+2 7])) = .5;
        T(i,i) = .8;
end

imagesc(T)
prettify(gcf)
caxis([0 1])
colormap('bone')
saveas(h,fdir + "/Figure_5_icon_low_t_STS.emf")
saveas(h,fdir + "/Figure_5_icon_low_t_STS.png")

h = figure
T = zeros(7,7);
for i = 1:7

    T(i,max([i-1 1]):min([i+1 7])) = .5;
        T(i,i) = .8;
end

imagesc(T)
prettify(gcf)
caxis([0 1])
colormap('bone')
saveas(h,fdir + "/Figure_5_icon_low_t_STS_light.emf")
saveas(h,fdir + "/Figure_5_icon_low_t_STS_light.png")


h = figure
T = zeros(7,7);
for i = 1:7

    T(i,max([i-3 1]):min([i+3 7])) = .5;
        T(i,i) = .8;
end

imagesc(T)
prettify(gcf)
caxis([0 1])
colormap('bone')
saveas(h,fdir + "/Figure_5_icon_low_t_STS_heavy.emf")
saveas(h,fdir + "/Figure_5_icon_low_t_STS_heavy.png")




h = figure
T = zeros(7,7);
for i = 1:7

    T(i,randi([1 7],1,2)) = .5;
    T(i,i) = .8;
end

imagesc(T)
prettify(gcf)
caxis([0 1])
colormap('bone')
saveas(h,fdir + "/Figure_5_icon_low_r_STS.emf")
saveas(h,fdir + "/Figure_5_icon_low_r_STS.png")



h = figure
T = zeros(7,7);
for i = 1:7

    T(i,randi([1 7],1,1)) = .5;
    T(i,i) = .8;
end

imagesc(T)
prettify(gcf)
caxis([0 1])
colormap('bone')
saveas(h,fdir + "/Figure_5_icon_low_r_STS_light.emf")
saveas(h,fdir + "/Figure_5_icon_low_r_STS_light.png")

h = figure
T = zeros(7,7);
for i = 1:7

    T(i,randi([1 7],1,4)) = .5;
    T(i,i) = .8;
end

imagesc(T)
prettify(gcf)
caxis([0 1])
colormap('bone')
saveas(h,fdir + "/Figure_5_icon_low_r_STS_heavy.emf")
saveas(h,fdir + "/Figure_5_icon_low_r_STS_heavy.png")