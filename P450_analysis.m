%% Analysis of mutational scanning in Matlab
% This script is part of a chapter for Methods in MolecularBiology. Please cite it.
% This file and its documentation is available at https: / / github.com / matteoferla / mutational_scanning
% This script was tested under R2016a
% It requires http://uk.mathworks.com/matlabcentral/fileexchange/7210-ternary-plots


%% Import data
% data is normally stored in Excel spreadsheets and the quickest way is to copy-paste the data of interest in the variable editor.
% First, open Excel and copy. In Mat make a blank variable (array or cell array) open it, 
% and paste or right click and 'paste from Excel'.
colname = {}; openvar('colname');
rowname = {}; openvar('rowname');
stackname = {}; openvar('stackname');
% the copy-paste trick does not work with >2 order arrays, so 2nd order
% matrix first.
% In this case multiple competing activities are being checked.
% In this case, each table in the stack is pasted to the right of the previous.
% predata is 20x(4*4) while data is 20x4x4
predata = []; openvar('predata');
% convert predata to a 3rd order tensor of dimensions 20, 4, 4.
data = reshape(predata, 20, 4, 4);
% Alway check when you reshape data
display([predata(1, 5), data(1, 1, 2)]);

% if you are using the demo values, load the mat file.
load('demo_data.mat')
% if it does not work you might be in the wrong directory ?chance current folder then? 
% or have unicode ("fancy letters") in your path.
%% derived data
% The data is in %. Let's make it fractional.
fdata = data / 100;
% number of mutants
Nmuts = numel(colname) * numel(rowname);

%% build colors and size
% black for 1 conversion
% red for 2-selectivity
% green for 4-selectivity
% blue for 3-selectivity
% let's make a scatter plot with dot size.
% RGB... this ought to be done with a reshape, but this uncoooth way is easier. 
colore=reshape(data(81:end),80,3);
% scale min = 0 max = 1
scolore = (colore - repmat(min(colore), 80, 1)) ./ (repmat(max(colore), 80, 1) - repmat(min(colore), 80, 1));
% gamma correction mu=0.5    x_gcorr=x^gamma
%{
<latex>
\begin{equation} 
\boldsymbol{\xi}=\left( \frac{\boldsymbol{x}-min(\boldsymbol{x})}{max(\boldsymbol{x})-min(\boldsymbol{x})}\right)^\gamma
\textrm{ where }\gamma=\frac{log(0.5)}{log(median(\boldsymbol{x}))}
\end{equation}
</latex>
%}
gammas=log(0.5)./log(median(scolore));
gcolore=scolore.^repmat(gammas, 80,1);

%% Scatter plot
% This is a fake scatter plot. That is I pick the (x,y) to draw what I
% like.
% grid
[x, y] = meshgrid(1:numel(colname), 1:numel(rowname));
% size based on activity
dotsize = fdata(:, :, 1) .* 200;
figure;
% raw?
% scatter(y(:), x(:), dotsize(:), colore, 'filled');
% scaled?
% scatter(y(:), x(:), dotsize(:), scolore, 'filled');
% scaled with gamma correction
scatter(y(:), x(:), dotsize(:), gcolore, 'filled');
hold on
% border.
scatter(y(:), x(:), dotsize(:), 'k');
ylim([0, 5.5])
xlim([0, 21])
ax = gca;
ax.YTick = 1:numel(colname);
ax.YTickLabel = colname;
ax.XTick = 1:numel(rowname);
ax.XTickLabel = rowname;
%title({'Relative activities of the mutants', 'Colors = selectivities scales $\frac{x-min}{max-min}$', 'Size= substrate usage'}, 'interpreter','latex');
title({'Relative activities of the mutants', 'Colors = selectivities scaled by (\frac{x-min}{max-min})^gamma where mean^gamma=0.5', 'Size= substrate usage'});

ax2=axes('Position',[.6 .8 .2 .1]);
% box on
% dotsize = fdata(:, :, 1) .* 200;
scatter(1:4, ones(4,1), 50:50:200, 'k');
axis off
text(1:4, ones(4,1)/3,{'25%','50%','75%','100%'})
text(-2, 1, 'activity rel. to wt')



%% Heatmap

figure
K = 1 - reshape(fdata(1:80),20,4)';  % blackness, not whiteness is activity
R = reshape(gcolore(1:80),20,4)';
G = reshape(gcolore(81:160),20,4)';
B = reshape(gcolore(161:240),20,4)';
colore2 = [cat(3, K, K, K); cat(3, R, G, B)];
image(colore2);
ax = gca;
ax.YTick = 1:8;
ax.YTickLabel = reshape(repmat(colname, 1, 2), 8, 1);
ax.XTick = 1:numel(rowname);
ax.XTickLabel = rowname;
title({'Heatmap with black for conversion, red for 2-selectivity','green for 4-selectivity and blue for 3-selectivity'})

%% Rank
% this is a bad idea. Here for illustration why.
figure
[Ks, Ki]=sort(data(1:80));
K = reshape(Ki,4,20);
[Rs, Ri]=sort(data(81:160));
R = reshape(Ri,4,20);
[Gs, Gi]=sort(data(161:240));
G = reshape(Gi,4,20);
[Bs, Bi]=sort(data(241:320));
B = reshape(Bi,4,20);
colore2 = [cat(3, K, K, K)./80; cat(3, R, G, B)./80];
image(colore2);
ax = gca;
ax.YTick = 1:8;
ax.YTickLabel = reshape(repmat(colname, 1, 2), 8, 1);
ax.XTick = 1:numel(rowname);
ax.XTickLabel = rowname;

%% Ternary plot
% Version 1. With dots scaled to activity
figure
% Plot the ternary axis system
terplot;
% Plot the data
colormap(brewermap([], 'YlOrRd'))
hd = ternaryc(fdata(1 + Nmuts * 1:Nmuts * 2)',...
  fdata(1 + Nmuts * 2:Nmuts * 3)',...
  fdata(1 + Nmuts * 3:Nmuts * 4)',...
  reshape(data(:, :, 1), 80, 1), 'o');
% Add the labels
terlabel(stackname(2), stackname(3), stackname(4));
muts = strcat(repmat(colname',numel(rowname),1),repmat(rowname,numel(colname),1));
intensita = reshape(fdata(:, :, 1), 80, 1);
for hdii = 1:numel(hd)
    hdi = hd(hdii);
    text(hdi.XData + 0.02, hdi.YData, muts(hdii));
    hdi.MarkerEdgeColor = 'k';
    hdi.MarkerSize = 10 * intensita(hdii);
end
c=colorbar;
c.Label.String = '% Activity relative to wild type';

% Version 2. With dots not scaled
figure;
% Plot the ternary axis system
[h, hg, htick] = terplot;
% Plot the data
% First set the colormap (can't be done afterwards)
colormap(brewermap([], 'YlOrRd'))
[hd, hcb] = ternaryc(fdata(1 + Nmuts * 1:Nmuts * 2)',...
  fdata(1 + Nmuts * 2:Nmuts * 3)',...
  fdata(1 + Nmuts * 3:Nmuts * 4)',...
  reshape(data(:, :, 1), 80, 1), 'o');
% Add the labels
hlabels = terlabel(stackname(2), stackname(3), stackname(4));
muts = strcat(reshape(repmat(colname,numel(rowname),1),80,1),repmat(rowname,numel(colname),1));
intensita = reshape(fdata(:, :, 1), 80, 1);
txt=cell(numel(hd),1);
txtx=zeros(numel(hd),1);
txty=zeros(numel(hd),1);
for hdii = 1:numel(hd)
    hdi = hd(hdii);
    txtx(hdii)=hdi.XData + 0.02;
    txty(hdii)=hdi.YData;
    hdi.MarkerEdgeColor = 'k';
    %hdi.MarkerSize=10*intensita(hdii);
end
% place the text next to the dots
text(txtx,txty,muts);
% or try and be fancy and fit them in nicely (fun from SO not Exchange):
% textfit(txtx,txty,muts);
c=colorbar;
c.Label.String = '% Activity relative to wild type';

