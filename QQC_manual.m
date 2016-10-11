
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
system(sprintf('/usr/local/bin/convert_trace -out_format scf < "%s" > temp_via_matlab.scf', fullfile(pathname,file)));
% although do change your Matlab's list of path if you haven't.
[Sample, Probability] = scfread('temp_via_matlab.scf');

%% Unicode in path or something as bad?
% Change Current folder...
[Sample, Probability] = scfread(fullfile('temp.scf'));

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
    pa=Probability.peak_index(zone(1)+ni(i)-1)-floor(interpeak/2*0.70);
    pb=Probability.peak_index(zone(1)+ni(i)-1)+ceil(interpeak/2*0.70);
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
codon=m2;
codon(codon<=0)=0.0001;
figure;
set(groot,'defaultAxesColorOrder',chromamap)
for i=1:3
subplot(2,3,i)
pa=Probability.peak_index(zone(1)+ni(i)-1)-floor(interpeak/2);
pb=Probability.peak_index(zone(1)+ni(i)-1)+ceil(interpeak/2);
chroma=plot([Sample.A(pa:pb), Sample.T(pa:pb), Sample.G(pa:pb), Sample.C(pa:pb)]);
g=gca;
w=floor(interpeak/2)+ceil(interpeak/2);
line([w/2*0.3+1 w/2*0.3+1],g.YLim,'Color','k');
line([w/2*1.7+1 w/2*1.7+1],g.YLim,'Color','k');
subplot(2,3,i+3)
% The codons were per row to make humans unaccostomed to matlab happy...
slices=pie(codon(i,:),{'A','T','G','C'});
colormap(chromamap)
end
suptitle(sprintf('Qpool= %0.2f for scheme %s',Qpool,scheme))


%% AA composition
bases={'A','T','G','C'};
% cellfun(@strcat,bases',bases) gives an error.
codnames=cell(4,4,4);
for i=1:4
    for j=1:4
     codnames(i,j,:)=strcat(bases{i},bases{j},bases);
    end
end

% tensor product
codprob=bsxfun(@mtimes, codon(1,:)'*codon(2,:), reshape(codon(3,:),1,1,4));
% proof of no errors.
% sum(codprob(:)) is one.
% codnames(find(codprob==max(codprob(:)))) is the best according to
% figure 2.
codtable=[array2table(codprob(:),'RowNames',codnames(:),'VariableNames', {'Prob'}), cell2table(nt2aa(codnames(:)),'RowNames',codnames(:),'VariableNames', {'AA'})];
aaprob=grpstats(codtable,'AA',@sum);
saaprob=sortrows(aaprob,'AA');
display(aaprob)

%ref
codpred=bsxfun(@mtimes, scheme_pred(1,:)'*scheme_pred(2,:), reshape(scheme_pred(3,:),1,1,4));
predtable=[array2table(codpred(:),'RowNames',codnames(:),'VariableNames', {'Prob'}), cell2table(nt2aa(codnames(:)),'RowNames',codnames(:),'VariableNames', {'AA'})];
aapred=grpstats(predtable,'AA',@sum);
saapred=sortrows(aapred,'AA');

figure;
bar([saaprob.sum_Prob, saapred.sum_Prob])
ax=gca;
ax.XTick=1:21;
ax.XTickLabel=saaprob.AA;
ax.XTickLabelRotation=45;
title(file)
legend({'Experimental','Expected'})


%% Disentangle mixed codons
% This section is experimental.

% Some schemes, use mixes of codons.
% 22c is NDT?+?9 eq. VHG?+?1 eq. TGG
% Tang is 12 eq. NDT?+?6 eq. VHA?+?1 eq. TGG?+?1 eq. ATG
% There are two ways. Solve the equation or just Suduku it.
% Fitting the equation would be better but requires writing the matrices.

% Tang
NDT = [ones(1, 4)./4; [1/3 1/3 1/3 0]; [0 1 0 0]];
VHA = [1/3 0 1/3 1/3; 1/3 1/3 0 1/3; 1 0 0 0];
TGG = [0 1 0 0; 0 0 1 0; 0 0 1 0];
ATG  = [1 0 0 0; 0 1 0 0; 0 0 1 0];
con = cat(3, NDT, VHA, TGG, ATG);
pro = [12 6 1 1];

% 22c
VHG = [1/3 0 1/3 1/3; 1/3 1/3 0 1/3; 0 0 1 0];
con = cat(3, NDT, VHG, TGG);
pro = [1 9 1];

ppro = squeeze(pro ./sum(pro));

objfunc = @(offness) sum(abs(reshape(abs(offness(:,:,1) .* ppro(1) .* con(:,:,1)) + abs(offness(:,:,2) .* ppro(2) .* con(:,:,2)) + abs(offness(:,:,3) .* ppro(3) .* con(:,:,3)) - codon, 12,1)));
% ideally a weight would be good to balance the two.
[offness,fval] = fminsearch(objfunc,ones(3,4,numel(ppro)),optimset('MaxFunEvals',5e4));
offness=abs(offness);
pcodon=abs(offness(:,:,1) .* ppro(1) .* con(:,:,1)) + abs(offness(:,:,2) .* ppro(2) .* con(:,:,2)) + abs(offness(:,:,3) .* ppro(3) .* con(:,:,3));

figure;
image(cat(3,pcodon,codon,zeros(3,4)))
ax=gca;
ax.XTick = 1:4;
ax.XTickLabel = bases;
ax.YTick = 1:3;
ax.YTickLabel = {'Pos1','Pos2','Pos3'};
title({'redness=predicted by model', 'greenness=actual'})

% tensor product
codprob=zeros(4,4,4);
codpred=zeros(4,4,4);


for i=1:numel(ppro)
    dcodon=con(:,:,i);
    tcodon=offness(:,:,i) .* dcodon;
    codprob=codprob+ppro(i)*bsxfun(@mtimes, tcodon(1,:)'*tcodon(2,:), reshape(tcodon(3,:),1,1,4));
    codpred=codpred+ppro(i)*bsxfun(@mtimes, dcodon(1,:)'*dcodon(2,:), reshape(dcodon(3,:),1,1,4));
end
% proof of no errors.
% sum(codprob(:)) is one. it is 0.32
% codnames(find(codprob==max(codprob(:)))) is the best according to
% figure 2.
codtable=[array2table(codprob(:),'RowNames',codnames(:),'VariableNames', {'Prob'}), cell2table(nt2aa(codnames(:)),'RowNames',codnames(:),'VariableNames', {'AA'})];
aaprob=grpstats(codtable,'AA',@sum);
saaprob=sortrows(aaprob,'AA');
display(aaprob)

%ref
predtable=[array2table(codpred(:),'RowNames',codnames(:),'VariableNames', {'Prob'}), cell2table(nt2aa(codnames(:)),'RowNames',codnames(:),'VariableNames', {'AA'})];
aapred=grpstats(predtable,'AA',@sum);
saapred=sortrows(aapred,'AA');

figure;
bar([saaprob.sum_Prob, saapred.sum_Prob])
ax=gca;
ax.XTick=1:21;
ax.XTickLabel=saaprob.AA;
ax.XTickLabelRotation=45;
title(file)
legend({'Experimental','Expected'})





