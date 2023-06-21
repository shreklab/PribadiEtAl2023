cd('E:\Dropbox\chalasanilabsync\mrieger\Manuscripts\PribadiEtAl-2022\FINALv2\Figures\FIG 8\')
d=readtable('Figure8A_dopaminoacidsequences.txt','Delimiter','\t','ReadVariableNames',true);
% Aligment performed on Clustal Omega Server at:
% https://www.ebi.ac.uk/Tools/msa/clustalo/
% with alignments downloaded in fasta format:
alignmentFasta=fastaread('dop-receptors_aligned.txt'); 

z=repmat('-',nrow(d),length(alignmentFasta(1).Sequence));
for i = 1:nrow(z)
    z(i,:) = alignmentFasta(i).Sequence;
end
wukabat = readtable('Fig8A_wukabatvariability.txt','Delimiter','\t','ReadVariableNames',true);
%wukabat variability determined on online server at: http://imed.med.ucm.es/PVS/
% download as tab delimited text file and read in.

alphabet='-GAILMVPFWYNCQSTDERHK';
alphabetFactor=nominal(cellstr(alphabet'),cellstr(alphabet'),cellstr(alphabet'));
alphabetFactorNumber=double(alphabetFactor);
alphabetcolors = hsv(length(alphabet)-1);
alphabetcolors = [0 0 0; alphabetcolors];


tmpCells = cell(nrow(z),ncol(z));
for i = 1:nrow(tmpCells)
    x = cellstr(z(i,:)')';
    tmpCells(i,:)=x;
end
tmpFactor=nominal(tmpCells,getlabels(alphabetFactor),getlabels(alphabetFactor));
tmpFactorNumber=double(tmpFactor);


figure;
subplot(4,1,1)
bar(wukabat.wukabat(1:500));
set(gca,'Xtick',[],'TickDir','out',...
    'Xtick',[1 500],'Xticklabel',sub(1:500,[1 500]),...
    'Xlim',[1 500],'Ylim',[0 15],'Ytick',[5 10 15]);
box off;
subplot(4,1,2)
imagesc(tmpFactorNumber(:,1:500));
colormap(alphabetcolors);
box off;
set(gca,'Ytick',1:nrow(z),'TickDir','out',...
    'Xtick',[1,50:50:500],'Xticklabel',[1,50:50:500],'Xlim',[1 500]);

subplot(4,1,3)
bar(wukabat.wukabat(501:ncol(z)));
set(gca,'Xtick',[],'TickDir','out',...
    'Xtick',[1 ncol(z)-501+1],'Xticklabel',sub(501:ncol(z),[1 ncol(z)-501+1]),...
    'Xlim',[1 500],'Ylim',[0 15],'Ytick',[5 10 15]);
box off;
subplot(4,1,4)
imagesc(tmpFactorNumber(:,501:ncol(z)));
colormap(alphabetcolors);
box off;
set(gca,'Ytick',1:nrow(z),'TickDir','out',...
    'Xtick',[1,50,100,ncol(z)-501+1],'Xticklabel',[1,50,100,ncol(z)-501+1]+500,'Xlim',[1 500]);

figure;imagesc(alphabetFactorNumber);
colormap(alphabetcolors);
set(gca,'TickDir','out','Xtick',[],'Ytick',1:nrow(alphabetFactorNumber),'YTickLabel',...
    cellstr(alphabetFactor));
box off;
