clc
clear

startPos = 1;
endPos = 1;
nSymbols = 4;

%Prosze odkomentowac wczytywanie wybranego pliku fasta z nukleotydami:
Data=fastaread('nukleotydy1.fasta');
%Data=fastaread('nukleotydy2.fasta');

test=cell(1,length(Data));
for N=1:length(Data)
test{N}=Data(N).Sequence;
end

%%%  a,t,c,g,

count_A=zeros(length(test{1}),1);
count_T=zeros(length(test{1}),1);
count_C=zeros(length(test{1}),1);
count_G=zeros(length(test{1}),1);

for M=1:length(test)
    
    for N=1:length(test{1})
        if strcmp(test{M}(N),'a')
            count_A(N)=count_A(N)+1;
        elseif strcmp(test{M}(N),'t')
            count_T(N)=count_T(N)+1;
        elseif strcmp(test{M}(N),'c')
            count_C(N)=count_C(N)+1;
        elseif strcmp(test{M}(N),'g')
            count_G(N)=count_G(N)+1;          
        end 
    end 
end 


wyrazy = repmat({''}, 1, length(test{1}));
for ii = 1:length(test{1})

    for jj = 1:count_A(ii)  
        wyrazy(ii) = append(wyrazy(ii),'A');
    end
    for jj = 1:count_C(ii)
        wyrazy(ii) = append(wyrazy(ii),'C');
    end
    for jj = 1:count_G(ii)
        wyrazy(ii) = append(wyrazy(ii),'G');
    end
    for jj = 1:count_T(ii)
        wyrazy(ii) = append(wyrazy(ii),'T');
    end
    
end

P = wyrazy;

P = P(:);
P = strrep(P,' ','-'); % usuniecie spacji i myslnikow 
P = char(P); % tablica znakow
P = P';
seqs = upper(P);
[numSeq, nPos] = size(seqs); 
uniqueList = unique(seqs);
m = length(uniqueList);
pcM = zeros(m, nPos);
for i = 1:nPos
    for j = 1:m
        pcM(j,i) = sum(seqs(:,i) == uniqueList(j));
    end
end
    
% Oblicz macierz wag do wyswietalnia seqlogo
freqM = [];
[symbolList, tmpIdx] = regexpi(uniqueList', '[A-Z]', 'match');
symbolList = char(symbolList');
 
if ~isempty(tmpIdx)
    for i = 1:length(tmpIdx)
        freqM(i, :) = pcM(tmpIdx(i),:);  %#ok
    end
end

% czestotliwosc symbolu w okreslonej pozycji sekwencji
freqM = freqM/numSeq;

%maxLen - maksymalna dlugosc sekwencji
maxLen = size(freqM, 2);

if endPos == 1
    endPos = maxLen;
end
freqM = freqM(:, startPos:endPos);
wtM = freqM; 
S_before = log2(nSymbols);
freqM(freqM == 0) = 1;

% niepewnosc kazdej pozycji
S_after = -sum(log2(freqM).*freqM, 1);
R = S_before - S_after;
nPos = (endPos - startPos) + 1;
for i =1:nPos
    wtM(:, i) = wtM(:, i) * R(i);
end

hFigure = seqshowlogo(wtM, symbolList, true, startPos);
   
function hFigure = seqshowlogo(varargin)

seqType = 'NT';
wtMatrix = [];
symbols = [];
startPos = 1;

if nargin == 4 
    wtMatrix = varargin{1};
    symbols = varargin{2};
end

import com.mathworks.toolbox.bioinfo.sequence.*;
import com.mathworks.mwswing.MJScrollPane;

% Viewer
logoViewer = SequenceViewer(wtMatrix, symbols, startPos, seqType);
awtinvoke(logoViewer,'addSeqLogo()');
scrollpanel = MJScrollPane(logoViewer, MJScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,MJScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);

% Figure seqlogo
hFigure = figure('WindowStyle', 'normal','Resize', 'on','Toolbar', 'none','NumberTitle','off','Name', 'Sequence Logo');
        
% Dopasowanie figure
d = awtinvoke(scrollpanel, 'getPreferredSize()');
pos = getpixelposition(hFigure);
pos(3) = d.getWidth;
pos(4) = d.getHeight;
setpixelposition(hFigure,pos);
[logoP, logoContainer] = javacomponent(scrollpanel, [0, 0, pos(3), pos(4)], hFigure);
set(logoContainer, 'units', 'normalized');
set(hFigure, 'visible', 'on')
end %seqshowlogo

%sortuj macierz wag
function [p,s] = sortWeightOrder(weight, symbollist)
[s, index] = sort(symbollist);
p=weight;
for i = 1:size(weight, 2)
    x=weight(:,i);
    p(:,i) = x(index);
end
end %sortWeightOrder

