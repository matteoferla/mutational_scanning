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

%% Refined
figure
A=interp1(Probability.peak_index(Probability.prob_A>0),Sample.A(Probability.peak_index(Probability.prob_A>0)),1:200:2e4);
T=interp1(Probability.peak_index(Probability.prob_T>0),Sample.T(Probability.peak_index(Probability.prob_T>0)),1:200:2e4);
G=interp1(Probability.peak_index(Probability.prob_G>0),Sample.G(Probability.peak_index(Probability.prob_G>0)),1:200:2e4);
C=interp1(Probability.peak_index(Probability.prob_C>0),Sample.C(Probability.peak_index(Probability.prob_C>0)),1:200:2e4);
plot([A; T; G; C]')
xlabel('Length of Sequence')

figure
A=spline(Probability.peak_index(Probability.prob_A>0),Sample.A(Probability.peak_index(Probability.prob_A>0)),100:200:2e4);
T=spline(Probability.peak_index(Probability.prob_T>0),Sample.T(Probability.peak_index(Probability.prob_T>0)),100:200:2e4);
G=spline(Probability.peak_index(Probability.prob_G>0),Sample.G(Probability.peak_index(Probability.prob_G>0)),100:200:2e4);
C=spline(Probability.peak_index(Probability.prob_C>0),Sample.C(Probability.peak_index(Probability.prob_C>0)),100:200:2e4);
plot([A; T; G; C]')
xlabel('Length of Sequence')
ylabel('RFU')
legend({'A','T','G','C'})

figure
A=interp1(Probability.peak_index(Probability.prob_A>0),smooth(Probability.peak_index(Probability.prob_A>0),Sample.A(Probability.peak_index(Probability.prob_A>0)),50,'rlowess'),100:200:2e4);
T=interp1(Probability.peak_index(Probability.prob_T>0),smooth(Probability.peak_index(Probability.prob_T>0),Sample.T(Probability.peak_index(Probability.prob_T>0)),50,'rlowess'),100:200:2e4);
G=interp1(Probability.peak_index(Probability.prob_G>0),smooth(Probability.peak_index(Probability.prob_G>0),Sample.G(Probability.peak_index(Probability.prob_G>0)),50,'rlowess'),100:200:2e4);
C=interp1(Probability.peak_index(Probability.prob_C>0),smooth(Probability.peak_index(Probability.prob_C>0),Sample.C(Probability.peak_index(Probability.prob_C>0)),50,'rlowess'),100:200:2e4);
plot([A; T; G; C]')
xlabel('Length of Sequence')
ylabel('RFU')
legend({'A','T','G','C'})

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
