%% Analysis of mutational scanning in Matlab
%{
    This script is part of a chapter for Methods in MolecularBiology. Please cite xxx.
    This file and its documentation is available at https://github.com/matteoferla/mutational_scanning
    This is here just for now. I will add a better documented script in the coming week
        %}

%% Import data
% data is in xlsx, so open Excel and copy paste in Mat with right click paste.
predata=[];  openvar('predata');  % each table in the stack is pasted to the right of the previous.
colname={};  openvar('colname');
rowname={};  openvar('rowname');
stackname={};openvar('stackname');

% data saved in mat file.
load('demo_data.mat')

% predata is 20x(4*4) while data is 20x4x4

data=reshape(predata,20,4,4);
%check
display([predata(1,5),data(1,1,2)]);

%% build colors and size
% black for 1 conversion
% red for 2-selectivity
% green for 4-selectivity
% blue for 3-selectivity
% Ehm that is four colours. I cannot do four only three, unless you see UV.
% Black square round each would not work.
% let's not make a image(). Let's make a scatter where circle size is *1*.

fdata=data/100;
[x,y] = meshgrid(1:numel(colname),1:numel(rowname));
dotsize=fdata(:,:,1).*200;
% RGB
Nmuts=numel(colname)*numel(rowname);
% there must be a more elegant way. schmeh.
colore=[data(1+Nmuts*1:Nmuts*2)', data(1+Nmuts*2:Nmuts*3)', data(1+Nmuts*3:Nmuts*4)'];

%% plot!
scatter(x(:),y(:),dotsize(:),colore,'filled');
xlim([0,5])
ylim([0,21])
ax=gca;
ax.XTick = 1:numel(colname);
ax.XTickLabel = colname;
ax.YTick = 1:numel(rowname);
ax.YTickLabel = rowname;

%% plot transposed...
figure;
scatter(y(:),x(:),dotsize(:),colore,'filled');
ylim([0,5])
xlim([0,21])
ax=gca;
ax.YTick = 1:numel(colname);
ax.YTickLabel = colname;
ax.XTick = 1:numel(rowname);
ax.XTickLabel = rowname;


%% plot3
% why is it yellow?
% Y=R+G. W=R+G+B.
% maybe adding a border?
figure;
scatter(y(:),x(:),dotsize(:),colore,'filled');
hold on
scatter(y(:),x(:),dotsize(:),'k');
ylim([0,5])
xlim([0,21])
ax=gca;
ax.YTick = 1:numel(colname);
ax.YTickLabel = colname;
ax.XTick = 1:numel(rowname);
ax.XTickLabel = rowname;
%% 4x heatmaps...
%I think I failed to copy paste...
figure;
% black for 1 conversion
% red for 2-selectivity
% blue for 3-selectivity
% green for 4-selectivity
schemi={'Greys','Reds','Blues','Greens'};
for i=1:4
subplot(2,2,i)
imagesc(data(:,:,i))
colorbar;
ax=gca;
ax.XTick = 1:numel(colname);
ax.XTickLabel = colname;
ax.YTick = 1:numel(rowname);
ax.YTickLabel = rowname;
title(stackname{i});
colormap(ax,brewermap([],schemi{i}))
end

%okay, this odd.
sum(data(:,:,2:4),3)
data(:,:,5)=100-sum(data(:,:,2:4),3);

%% Ternary plot

% ah. The sum of the three selectvities are 100%.
figure
% Plot the ternary axis system
terplot;
% Plot the data
% First set the colormap (can't be done afterwards)
colormap(brewermap([],'YlOrRd'))
[hd,hcb]=ternaryc(fdata(1+Nmuts*1:Nmuts*2)',...
    fdata(1+Nmuts*2:Nmuts*3)',...
    fdata(1+Nmuts*3:Nmuts*4)',...
    reshape(data(:,:,1),80,1),'o');
% Add the labels
terlabel(stackname(2),stackname(3),stackname(4));
muts=strcat(repmat(colname',numel(rowname),1),repmat(rowname,numel(colname),1));
intensita=reshape(fdata(:,:,1),80,1);
for hdii=1:numel(hd)
    hdi=hd(hdii);
    text(hdi.XData+0.02,hdi.YData,muts(hdii));
    hdi.MarkerEdgeColor='k';
    hdi.MarkerSize=10*intensita(hdii);
end


figure
% Plot the ternary axis system
[h,hg,htick]=terplot;
% Plot the data
% First set the colormap (can't be done afterwards)
colormap(brewermap([],'YlOrRd'))
[hd,hcb]=ternaryc(fdata(1+Nmuts*1:Nmuts*2)',...
    fdata(1+Nmuts*2:Nmuts*3)',...
    fdata(1+Nmuts*3:Nmuts*4)',...
    reshape(data(:,:,1),80,1),'o');
% Add the labels
hlabels=terlabel(stackname(2),stackname(3),stackname(4));
muts=strcat(repmat(colname',numel(rowname),1),repmat(rowname,numel(colname),1));
intensita=reshape(fdata(:,:,1),80,1);
for hdii=1:numel(hd)
    hdi=hd(hdii);
    text(hdi.XData+0.02,hdi.YData,muts(hdii));
    hdi.MarkerEdgeColor='k';
    %hdi.MarkerSize=10*intensita(hdii);
end

%% triangles?
figure;
%glyphplot(colore,'grid',[numel(rowname),numel(colname)],'obslabels',strcat(repmat(colname',numel(rowname),1),repmat(rowname,numel(colname),1)))
% scaled by the conversion.
glyphplot(colore.*repmat(reshape(data(:,:,1),80,1),1,3),'grid',[numel(rowname),numel(colname)],'obslabels',strcat(repmat(colname',numel(rowname),1),repmat(rowname,numel(colname),1)))

%% plots with scaling
figure;
scolore=(colore-repmat(min(colore),80,1))./(repmat(max(colore),80,1)-repmat(min(colore),80,1));
scatter(y(:),x(:),dotsize(:),scolore,'filled');
hold on
scatter(y(:),x(:),dotsize(:),'k');
ylim([0,5])
xlim([0,21])
ax=gca;
ax.YTick = 1:numel(colname);
ax.YTickLabel = colname;
ax.XTick = 1:numel(rowname);
ax.XTickLabel = rowname;
title({'Relative activities of the mutants','Colors = selectivities scales $\frac{x-min}{max-min}$','Size= substrate usage'});

%% AA order
% the mutations can be clustered by distance.
figure;
dendrogram(linkage(pdist(reshape(data(:,:,1:4),80,4))),'labels',reshape(strcat(repmat(colname',numel(rowname),1),repmat(rowname,numel(colname),1)),80,1));
title('mutations linked based on average distance of all four variables')

% Even better, the AA can be clustered. The order useful for the heatmaps.
figure;
dendrogram(linkage(pdist(reshape(data(:,:,1:4),20,16))),'labels',rowname);
title('With first pane (substrate use)')
figure;
dendrogram(linkage(pdist(reshape(data(:,:,2:4),20,12))),'labels',rowname);
title('Without first pane (substrate use)')
ax=gca;
niceOrder=ax.XTickLabel;
niceOrderIndex=zeros(20,1);
for ni=1:numel(niceOrder)
    niceOrderIndex(ni)=find(strcmp(rowname,niceOrder{ni}));
end

%% combo

figure;
subplot(2,1,1)
dendrogram(linkage(pdist(reshape(data(:,:,2:4),20,12))),'labels',rowname);
subplot(2,1,2)
niceOrder80=reshape(repmat(niceOrderIndex,1,4)+repmat(0:3,20,1)*20,80,1);
niceOrder320=reshape(repmat(niceOrderIndex,1,16)+repmat(0:15,20,1)*20,320,1);
scolore=(colore-repmat(min(colore),80,1))./(repmat(max(colore),80,1)-repmat(min(colore),80,1));
scatter(y(:),x(:),dotsize(niceOrder80),scolore(niceOrder80,:),'filled');
hold on
scatter(y(:),x(:),dotsize(niceOrder80),'k');
ylim([0,5])
xlim([0,21])
ax=gca;
ax.YTick = 1:numel(colname);
ax.YTickLabel = colname;
ax.XTick = 1:numel(rowname);
ax.XTickLabel = niceOrder;
%title({'Relative activities of the mutants','Colors = selectivities scales $\frac{x-min}{max-min}$','Size= substrate usage'});


%% requested heatmap

combo=[data(:,:,1)'; data(:,:,2)'; data(:,:,3)'; data(:,:,4)'];
K=(data(:,:,1)'-min(data(1:80)))/(max(data(1:80))-min(data(1:80)));
R=(data(:,:,2)'-min(data(81:160)))/(max(data(81:160))-min(data(81:160)));
G=(data(:,:,3)'-min(data(161:240)))/(max(data(81:160))-min(data(81:160)));
B=(data(:,:,4)'-min(data(241:320)))/(max(data(241:320))-min(data(241:320)));
combocol=[cat(3,K,K,K);...
    cat(3,R,zeros(4,20),zeros(4,20));...
    cat(3,zeros(4,20),G,zeros(4,20));...
    cat(3,zeros(4,20),zeros(4,20),B);...
    ];
image(combocol);
ax=gca;
ax.YTick = 1:16;
ax.YTickLabel = reshape(repmat(colname,1,4),16,1);
ax.XTick = 1:numel(rowname);
ax.XTickLabel = rowname;


