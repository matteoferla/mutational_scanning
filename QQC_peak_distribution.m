%% Peak height
% One assumption is that peak height is constant which it isn't.
% So the first question is: do different channels behave differently?

%% load
[Sample, Probability] = scfread('temp_via_matlab.scf');
chromamap=[0.4660 0.6740 0.1880;
    0.8500    0.3250    0.0980;
    0.3250    0.3250    0.3250;
    0         0.4470    0.7410];
set(groot,'defaultAxesColorOrder',chromamap)

%% Raw data
figure
plot(Probability.peak_index(Probability.prob_A>0),Sample.A(Probability.peak_index(Probability.prob_A>0)),'.')
hold on
plot(Probability.peak_index(Probability.prob_T>0),Sample.T(Probability.peak_index(Probability.prob_T>0)),'.')
plot(Probability.peak_index(Probability.prob_G>0),Sample.G(Probability.peak_index(Probability.prob_G>0)),'.')
plot(Probability.peak_index(Probability.prob_C>0),Sample.C(Probability.peak_index(Probability.prob_C>0)),'.')
xlabel('sequential number of appearance ?not interpolated!')
ylabel('RFU')
legend({'A','T','G','C'})
title('Raw data')

%% Refined
figure
A=interp1(Probability.peak_index(Probability.prob_A>0),smooth(Probability.peak_index(Probability.prob_A>0),Sample.A(Probability.peak_index(Probability.prob_A>0)),50,'rlowess'),100:200:2e4);
T=interp1(Probability.peak_index(Probability.prob_T>0),smooth(Probability.peak_index(Probability.prob_T>0),Sample.T(Probability.peak_index(Probability.prob_T>0)),50,'rlowess'),100:200:2e4);
G=interp1(Probability.peak_index(Probability.prob_G>0),smooth(Probability.peak_index(Probability.prob_G>0),Sample.G(Probability.peak_index(Probability.prob_G>0)),50,'rlowess'),100:200:2e4);
C=interp1(Probability.peak_index(Probability.prob_C>0),smooth(Probability.peak_index(Probability.prob_C>0),Sample.C(Probability.peak_index(Probability.prob_C>0)),50,'rlowess'),100:200:2e4);
plot([A; T; G; C]')
xlabel('Length of Sequence')
ylabel('RFU')
legend({'A','T','G','C'})
title('Smoothed')
%% fit
% Is it best to fit to a lognormal or gamma?
% I'll go polynomial for now.
abscissa=Probability.peak_index(Probability.prob_A>0);
ordinate=Sample.A(Probability.peak_index(Probability.prob_A>0));
p=polyfit(abscissa,ordinate, 3);
f1 = polyval(p,abscissa);
figure
plot([abscissa, abscissa] ,[ordinate, f1])

% all together
abscissa=[Probability.peak_index(Probability.prob_A>0);Probability.peak_index(Probability.prob_T>0); Probability.peak_index(Probability.prob_G>0);Probability.peak_index(Probability.prob_C>0)];
[a,ai]=sort(abscissa);
ordinate=[Sample.A(Probability.peak_index(Probability.prob_A>0)); Sample.T(Probability.peak_index(Probability.prob_T>0));Sample.G(Probability.peak_index(Probability.prob_G>0));Sample.C(Probability.peak_index(Probability.prob_C>0))];
p=polyfit(abscissa(ai),ordinate(ai), 3);
f1 = polyval(p,abscissa(ai));
figure
plot(abscissa(ai) ,ordinate(ai),'.b')
hold on
plot(abscissa(ai), f1,'-r')
% the data is too noisy...




% Proportions
figure
pie(sum([Probability.prob_A>0, Probability.prob_T>0, Probability.prob_G>0, Probability.prob_C>0])/1821)
legend({'A','T','G','C'})

figure
A=Sample.A(Probability.peak_index(Probability.prob_A>0));
T=Sample.T(Probability.peak_index(Probability.prob_T>0));
G=Sample.G(Probability.peak_index(Probability.prob_G>0));
C=Sample.C(Probability.peak_index(Probability.prob_C>0));
bar(cat(1,[mean(A), mean(T), mean(G), mean(C)]/mean([A; T; G; C])/4,[sum(Probability.prob_A>0),sum(Probability.prob_T>0),sum(Probability.prob_G>0),sum(Probability.prob_C>0)]/1821)')
ax=gca;
ax.XTickLabel={'A','T','G','C'};
legend({'Mean peak intensity','Number of calls'})

figure
violin({A,T,G,C})

%% Distance
% one thing that seems true is that peak height vs. distance from previous
% might correlate.
gap=(numel(Sample.A)/numel(Probability.prob_A));
figure
x=Probability.peak_index(Probability.prob_A>0);
plot((x(2:end)-x(1:end-1))/gap,Sample.A(x(2:end)),'.')
hold on
x=Probability.peak_index(Probability.prob_T>0);
plot((x(2:end)-x(1:end-1))/gap,Sample.T(x(2:end)),'.')
x=Probability.peak_index(Probability.prob_G>0);
plot((x(2:end)-x(1:end-1))/gap,Sample.G(x(2:end)),'.')
x=Probability.peak_index(Probability.prob_C>0);
plot((x(2:end)-x(1:end-1))/gap,Sample.C(x(2:end)),'.')
legend({'A','T','G','C'});
xlabel('Distance from upstream base of same type')
ylabel('Intensity')
ylim([0 1000])
xlim([0 20])


% but I know there is a trend with distance too.
figure
x=Probability.peak_index(Probability.prob_A>0);
plot((x(2:end)-x(1:end-1))/gap,Sample.A(x(2:end))./polyval(p,x(2:end)),'.')
hold on
x=Probability.peak_index(Probability.prob_T>0);
plot((x(2:end)-x(1:end-1))/gap,Sample.T(x(2:end))./polyval(p,x(2:end)),'.')
x=Probability.peak_index(Probability.prob_G>0);
plot((x(2:end)-x(1:end-1))/gap,Sample.G(x(2:end))./polyval(p,x(2:end)),'.')
x=Probability.peak_index(Probability.prob_C>0);
plot((x(2:end)-x(1:end-1))/gap,Sample.C(x(2:end))./polyval(p,x(2:end)),'.')
legend({'A','T','G','C'});
xlabel('Distance from upstream base of same type')
ylabel('Intensity')
ylim([0 2])
xlim([0 20])

