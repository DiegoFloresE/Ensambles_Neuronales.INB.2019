close all
clc

load ('Rasterbin.mat')
X = Rasterbin';
Y = Rasterbin;
Z = X*Y;
normX = sqrt(sum(X.^2,2));
normY = sqrt(sum(Y.^2,1));
C = bsxfun(@rdivide, bsxfun(@rdivide, Z, normX), normY);
x = size (X,1);
y = size (Y,2);

imagesc (x,y,C);
colormap jet;
title('Similarity map','fontsize',16');
xticks(1000:200:2000);xticklabels({'0','200','400','600','800','1000'});
yticks(1000:200:2000);yticklabels({'0','200','400','600','800','1000'});
c=colorbar;
c.Label.String = 'similarity index';
c.FontSize = 14;
c.Ticks =[0 1];
xlabel('vectors','fontsize',14);
ylabel('vectors','fontsize',14);
