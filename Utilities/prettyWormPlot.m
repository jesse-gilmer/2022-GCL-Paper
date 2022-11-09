function [] = prettyWormPlot(D,Cat,opt,c)
opt2.scatter = 0;
opt2.SE = 0;

f = fieldnames(opt);
for i = 1:length(f)
        if length(opt.(f{i})) < 1
        opt.(f{i}) = opt2.(f{i})
        end
end

darkB = c./1.2;
lightB = c;
ZSCORE = 2.576;

CATS = unique(Cat);

hold on;
set(gcf,'color','w')
set(gca,'fontsize',14)
ii = 1;
for i = 1:length(CATS)
      [lx,ix] = find(Cat == CATS(i));
        MeX(ii) = nanmean(D(ix));
        SXE(ii) = nanstd(D(ix));
        if opt.SE == 1
            SXE(ii) = (std(D(ix))/sqrt(numel(D(ix))));
        end
        ii = ii + 1;
end

P = patch([CATS fliplr(CATS)],[(MeX-SXE) fliplr(MeX+SXE)], lightB);
P.EdgeAlpha = 0;
P.FaceColor = lightB;
P.FaceAlpha = .5;

plot(CATS,MeX,'LineWidth',1.5,'Color',darkB)
scatter(CATS,MeX,25,darkB,'filled');

end