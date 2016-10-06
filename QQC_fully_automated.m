%% QQC Full automation
% This script is an example of fully automated analysis.
% It assumes the scheme is written in the file name and that all ab1 files
% in a folder are to be analysed.

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

%% Analysis presets
% Roughly were the mutation is:
zone = 370:390;
% The range is because the initial few bases may or may not be called and
% never come out consistent.
% It may be better giving a sequence to find in front or something, as everyone has a
% different need, drop matteo.ferla@gmail.com an email and he will change
% it for your need.

%% Script
%%% Parse
% Navigate to the right folder and that the folder free of unicode symbols
d = dir; %returns a structure
files = {d.name}; %names of all the files.
abfilelogi = cellfun(@sum, strfind(files, '.ab1')); %fudgelogical of ab1
fi=0;  % index.
% initiate blank table
Tmain  = cell2table(cell(0,7), 'VariableNames', {'A','T','G','C','sumDev','Qpool','scheme'});
% loop through
for f = files(abfilelogi > 0)
    display(f);
    fi=fi+1;
    
    %%% figure out the scheme.
    % in this case the user gives it in the sequencing file.
    % A T G C
    % see dx.doi.org/10.1038/srep10654
    if ~isempty(strfind(f{1},'NNN'))
        scheme='NNN';
        scheme_pred=ones(3, 4)./4;
    elseif ~isempty(strfind(f{1},'NNK'))
        scheme='NNK'; % K is keto, G or T. 
        scheme_pred=[ones(2, 4)./4; [0 0.5 0 0.5]];
    elseif ~isempty(strfind(f{1},'NNS'))
        scheme='NNS'; % S is strong, G or C. 
        scheme_pred=[ones(2, 4)./4; [0 0 0.5 0.5]];
    elseif ~isempty(strfind(f{1},'NDT'))
        scheme='NDT'; % D is not C. 
        scheme_pred=[ones(1, 4)./4; [1/3 1/3 1/3 0]; [0 1 0 0]];
    elseif ~isempty(strfind(f{1},'22C'))
        scheme='22C'; % see dx.doi.org/10.1038/srep10654
        scheme_pred=[[0.27 0.19 0.27 0.27]; [0.32 0.32 0.18 0.18]; [0 0.55 0.45 0]];
    elseif ~isempty(strfind(f{1},'20C'))
        scheme='19C'; % see dx.doi.org/10.1038/srep10654
        scheme_pred=[[0.30 0.20 0.25 0.25]; [0.35 0.25 0.20 0.20]; [0.30 0.60 0.10 0]];
    else
        error('Cannot recognise scheme from filename')
    end
    
    %%% get the scf....
    % Yes, the sudo bin strikes again!
    % In Mac El captain and greater stuff in not added to bin but usr/local/bin, so the
    % following fails...
    % system(sprintf('convert_trace -out_format scf < %s > temp.scf',f{1}));
    % this works
    system(sprintf('/usr/local/bin/convert_trace -out_format scf < %s > temp.scf', f{1}));
    % although do change your Matlab's list of $PATH
    [Sample, Probability] = scfread('temp.scf');
    %%% get range
    % find all the bases in the given range that are less than 0.9 pure.
    interpeak = mean(Probability.peak_index(2:end) - Probability.peak_index(1:end-1));
    noisy=max([Sample.A(Probability.peak_index(zone)),Sample.T(Probability.peak_index(zone)), Sample.C(Probability.peak_index(zone)), Sample.G(Probability.peak_index(zone))]')./sum([Sample.A(Probability.peak_index(zone)),Sample.T(Probability.peak_index(zone)), Sample.C(Probability.peak_index(zone)), Sample.G(Probability.peak_index(zone))]');
    ni=find(noisy<0.9);
    % crash if they are not 3 consecutive bases --an issue if NDT.
    if (numel(ni) ~= 3) || (max(ni)-min(ni) ~= 2)
        error('Noisy values not sequential')
    end
    m=zeros(3,4);
    ii=0;
    %%% Measure
    % loop per base, but in the chromatogram time axis.
    for i=(-1:1)*interpeak
        % base rough boundary. Fitting to gaussian may be better, but a
        % whole can of worms, due to bell curves merging.
        pa=Probability.peak_index(zone(1)+ni(1))+i;
        pb=Probability.peak_index(zone(1)+ni(1))+i+interpeak;
        ii=ii+1;
        % RFU of the peak top
        m(ii,1:4)=max([Sample.A(pa:pb), Sample.T(pa:pb), Sample.G(pa:pb), Sample.C(pa:pb)]);
    end
    %%% Analyse
    % make fraction of 1.
    m2=m./repmat(sum(m'),4,1)';
    % http://dx.doi.org/10.1016/j.enzmictec.2013.02.012
    deviation=sum(scheme_pred-abs(scheme_pred-m2),2);
    weights=sum(scheme_pred>0,2)./sum(scheme_pred(:)>0);
    wsum=sum(deviation.*weights);
    worse=sum(scheme_pred-abs(scheme_pred-[ones(3,1), zeros(3,3)]),2);
    wmin=sum(worse.*weights);
    Qpool=(wsum+abs(wmin))/(1+abs(wmin));
    %%% Store
    % "add" to the table Tmain.
    Tmain = [Tmain; [array2table([m2, deviation,repmat(Qpool,3,1)], 'VariableNames', {'A','T','G','C','sumDev','Qpool'},'RowNames',cellstr(strcat(f{1},' N',int2str(zone(1)+ni')))),cell2table(cellstr(repmat(scheme,3,1)),'VariableNames',{'scheme'})]];
end
%%% Output
display(Tmain)
writetable(Tmain,'summary.csv','WriteRowNames',true)