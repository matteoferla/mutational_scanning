
%% Prerequisite
% This needs Staden tools to convert ab1 to scf.
% Download from: https://sourceforge.net/projects/staden/files/
% unzip
% On a mac OS X open terminal
% navigate to the folder, probably this:
% $ cd ~/Downloads/io_lib-1.14.6
% read the instructuctions
% $ cat README
% do as they say
% $ ./configure
% $ make
% $ make install

%% Get scf via bash
% (If you want matlab to do this section see next one.)
% $ convert_trace -out_format scf < ACE-AA-088-01-55Â°C-BM3-A82_19C-T7-T7minus1.ab1 > test.scf
% the $ sign means in your terminal.
% Note the file name had the classic Unicode misinterpreted as ISO-8859-1
% when it was made ?Matlab is blameless this time.
[file,pathname]=uigetfile('scf');
[Sample, Probability] = scfread(fullfile(pathname,file));

%% Get scf via Matlab
[file,pathname]=uigetfile('ab1');
% read the ab1 file to be converted
% In Mac El captain and greater stuff in not added to bin but usr/local/bin, so the
% following fails...
% system(sprintf('convert_trace -out_format scf < %s > temp.scf',f{1}));
% this works
system(sprintf('/usr/local/bin/convert_trace -out_format scf < %s > temp.scf', fullfile(pathname,file)));
% although do change your Matlab's list of path if you haven't.
[Sample, Probability] = scfread(fullfile(pathname,'temp.scf'));

%% Scheme
% what scheme is used? % uncomment the appropriate...
% scheme='NNN';
% scheme_pred=ones(3, 4)./4;
% scheme='NNK'; % K is keto, G or T. 
% scheme_pred=[ones(2, 4)./4; [0 0.5 0 0.5]];
% scheme='NNS'; % S is strong, G or C. 
% scheme_pred=[ones(2, 4)./4; [0 0 0.5 0.5]];
% scheme='NDT'; % D is not C. 
% scheme_pred=[ones(1, 4)./4; [1/3 1/3 1/3 0]; [0 1 0 0]];
% scheme='22C'; % see dx.doi.org/10.1038/srep10654
% scheme_pred=[[0.27 0.19 0.27 0.27]; [0.32 0.32 0.18 0.18]; [0 0.55 0.45 0]];
% scheme='20C'; % see dx.doi.org/10.1038/srep10654 aka. Tang
% scheme_pred=[[0.30 0.20 0.25 0.25]; [0.35 0.25 0.20 0.20]; [0.30 0.60 0.10 0]];
% scheme='19C'; % Ask Carlos for reference.
% scheme_pred=[[0.263 0.263 0.211 0.263]; [0.368 0.263 0.211 0.158]; [0.368 0.053 0.105 0.474]];


%% Find mutation
% Roughly were the mutation is:
mutation = 378; % you give this! This example is specific to the example file!
zone = mutation-4:mutation+4;
% The range is because the initial few bases may or may not be called and
% never come out consistent.
% It may be better giving a sequence to find in front or something, as everyone has a
% different need, drop matteo.ferla@gmail.com an email and he will change
% it for your need.


%% Plot sequence zome.
% Abi traces are normally with this scheme.
% A in g
% T in r
% G in k
% C in b
% harsh primaries:
chromamap=[0 1 0; 1 0 0; 1 1 1; 0 0 1];
% Matlab default ColorOrder copied softer colors:
chromamap=[0.4660 0.6740 0.1880;
    0.8500    0.3250    0.0980;
    0.3250    0.3250    0.3250;
    0         0.4470    0.7410];
set(groot,'defaultAxesColorOrder',chromamap)
figure;
pa=Probability.peak_index(zone(1));
pb=Probability.peak_index(zone(end));
plot([Sample.A(pa:pb), Sample.T(pa:pb), Sample.G(pa:pb), Sample.C(pa:pb)]);
text(Probability.peak_index(zone)-Probability.peak_index(zone(1)),ones(numel(zone),1)*20,Probability.base(zone))
xlabel('Arbitrary chromatogram time')
ylabel('RFU')
legend({'A','T','G','C'})

%% get range
% find all the bases in the given range that are less than 0.9 pure.
interpeak = mean(Probability.peak_index(2:end) - Probability.peak_index(1:end-1));
noisy=max([Sample.A(Probability.peak_index(zone)),Sample.T(Probability.peak_index(zone)), Sample.C(Probability.peak_index(zone)), Sample.G(Probability.peak_index(zone))]')./sum([Sample.A(Probability.peak_index(zone)),Sample.T(Probability.peak_index(zone)), Sample.C(Probability.peak_index(zone)), Sample.G(Probability.peak_index(zone))]');
ni=find(noisy<0.9);
% crash if they are not 3 consecutive bases --an issue if NDT: so comment
% it out if needed.
if (numel(ni) ~= 3) || (max(ni)-min(ni) ~= 2)
    error('Noisy values not sequential')
end
m=zeros(3,4);
%% Measure
% loop per base, but convert to the chromatogram time axis (pa:pb).
for i=1:3
    % base rough boundary. Fitting to gaussian may be better, but a
    % whole can of worms, due to bell curves merging.
    pa=Probability.peak_index(zone(1)+ni(i)-1)-floor(interpeak/2);
    pb=Probability.peak_index(zone(1)+ni(i)-1)+ceil(interpeak/2);
    % RFU of the peak top
    m(i,1:4)=max([Sample.A(pa:pb), Sample.T(pa:pb), Sample.G(pa:pb), Sample.C(pa:pb)]);
end
%% Analyse
% make fraction of 1.
m2=m./repmat(sum(m'),4,1)';
% http://dx.doi.org/10.1016/j.enzmictec.2013.02.012
deviation=sum(scheme_pred-abs(scheme_pred-m2),2);
weights=sum(scheme_pred>0,2)./sum(scheme_pred(:)>0);
wsum=sum(deviation.*weights);
worse=sum(scheme_pred-abs(scheme_pred-[ones(3,1), zeros(3,3)]),2);
wmin=sum(worse.*weights);
Qpool=(wsum+abs(wmin))/(1+abs(wmin));

display(array2table([m2, deviation,repmat(Qpool,3,1)],...
    'VariableNames', {'A','T','G','C','sumDev','Qpool'},...
    'RowNames',cellstr(strcat(file,' N',int2str(zone(1)+ni')))));
%{
>> ans = 

                                                                 A           T           G          C        sumDev      Qpool 
                                                             _________    ________    _______    _______    ________    _______

    ACE-AA-088-01-55A?°C-BM3-A82_19C-T7-T7minus1.ab1 N377     0.010811    0.075676    0.38559    0.52793     0.17297    0.38074
    ACE-AA-088-01-55A?°C-BM3-A82_19C-T7-T7minus1.ab1 N378    0.0049669    0.071192    0.38411    0.53974     0.15232    0.38074
    ACE-AA-088-01-55A?°C-BM3-A82_19C-T7-T7minus1.ab1 N379    0.0042644           0    0.80597    0.18977    -0.11194    0.38074

%}
    

%% pie time
clean_m2=m2;
clean_m2(clean_m2<=0)=0.0001;
figure;
set(groot,'defaultAxesColorOrder',chromamap)
for i=1:3
subplot(2,3,i)
pa=Probability.peak_index(zone(1)+ni(i)-1)-floor(interpeak/2);
pb=Probability.peak_index(zone(1)+ni(i)-1)+ceil(interpeak/2);
chroma=plot([Sample.A(pa:pb), Sample.T(pa:pb), Sample.G(pa:pb), Sample.C(pa:pb)]);
subplot(2,3,i+3)
% The codons were per row to make humans unaccostomed to matlab happy...
slices=pie(clean_m2(i,:),{'A','T','G','C'});
colormap(chromamap)
end



