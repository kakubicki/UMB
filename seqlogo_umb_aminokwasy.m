clc
clear

%wyswietl seqlogo dla zestawu dopasowanych sekwencji nukleodytowych
startPos = 1;
endPos = 1;
nSymbols = 23;

%Prosze odkomentowac wczytywanie wybranego pliku fasta z aminokwasami:
%Data=fastaread('aminokwasy1.fasta');
Data=fastaread('aminokwasy2.fasta');

test=cell(1,length(Data));
for N=1:length(Data)
test{N}=Data(N).Sequence;
end

%%%  a, r, n, d, c, q, e, g, h, i, l, k, m, f, p, s, t, w, y, v

count_A=zeros(length(test{1}),1);
count_R=zeros(length(test{1}),1);
count_N=zeros(length(test{1}),1);
count_D=zeros(length(test{1}),1);
count_C=zeros(length(test{1}),1);
count_Q=zeros(length(test{1}),1);
count_E=zeros(length(test{1}),1);
count_G=zeros(length(test{1}),1);
count_H=zeros(length(test{1}),1);
count_I=zeros(length(test{1}),1);
count_L=zeros(length(test{1}),1);
count_K=zeros(length(test{1}),1);
count_M=zeros(length(test{1}),1);
count_F=zeros(length(test{1}),1);
count_P=zeros(length(test{1}),1);
count_S=zeros(length(test{1}),1);
count_T=zeros(length(test{1}),1);
count_W=zeros(length(test{1}),1);
count_Y=zeros(length(test{1}),1);
count_V=zeros(length(test{1}),1);

for M=1:length(test)
    
    for N=1:length(test{1})
        if strcmp(test{M}(N),'A')
            count_A(N)=count_A(N)+1;
        elseif strcmp(test{M}(N),'R')
            count_R(N)=count_R(N)+1;
        elseif strcmp(test{M}(N),'N')
            count_N(N)=count_N(N)+1;
        elseif strcmp(test{M}(N),'D')
            count_D(N)=count_D(N)+1;   
        elseif strcmp(test{M}(N),'C')
            count_C(N)=count_C(N)+1;
        elseif strcmp(test{M}(N),'Q')
            count_Q(N)=count_Q(N)+1;
        elseif strcmp(test{M}(N),'E')
            count_E(N)=count_E(N)+1;
        elseif strcmp(test{M}(N),'G')
            count_G(N)=count_G(N)+1;
        elseif strcmp(test{M}(N),'H')
            count_H(N)=count_H(N)+1;
        elseif strcmp(test{M}(N),'I')
            count_I(N)=count_I(N)+1;
        elseif strcmp(test{M}(N),'L')
            count_L(N)=count_L(N)+1;
        elseif strcmp(test{M}(N),'K')
            count_K(N)=count_K(N)+1;
        elseif strcmp(test{M}(N),'M')
            count_M(N)=count_M(N)+1;
        elseif strcmp(test{M}(N),'F')
            count_F(N)=count_F(N)+1;
        elseif strcmp(test{M}(N),'P')
            count_P(N)=count_P(N)+1;
        elseif strcmp(test{M}(N),'S')
            count_S(N)=count_S(N)+1;
        elseif strcmp(test{M}(N),'T')
            count_T(N)=count_T(N)+1;
        elseif strcmp(test{M}(N),'W')
            count_W(N)=count_W(N)+1;
        elseif strcmp(test{M}(N),'Y')
            count_Y(N)=count_Y(N)+1;
        elseif strcmp(test{M}(N),'V')
            count_V(N)=count_V(N)+1;         
        end 
    end 
end 


wyrazy = repmat({''}, 1, length(test{1}));
for ii = 1:length(test{1})

    for jj = 1:count_A(ii)  
        wyrazy(ii) = append(wyrazy(ii),'A');
    end
    for jj = 1:count_R(ii)
        wyrazy(ii) = append(wyrazy(ii),'R');
    end
    for jj = 1:count_N(ii)
        wyrazy(ii) = append(wyrazy(ii),'N');
    end
    for jj = 1:count_D(ii)
        wyrazy(ii) = append(wyrazy(ii),'D');
    end
    for jj = 1:count_C(ii)
        wyrazy(ii) = append(wyrazy(ii),'C');
    end
    for jj = 1:count_Q(ii)
        wyrazy(ii) = append(wyrazy(ii),'Q');
    end
    for jj = 1:count_E(ii)
        wyrazy(ii) = append(wyrazy(ii),'E');
    end
    for jj = 1:count_G(ii)
        wyrazy(ii) = append(wyrazy(ii),'G');
    end
    for jj = 1:count_H(ii)
        wyrazy(ii) = append(wyrazy(ii),'H');
    end
    for jj = 1:count_I(ii)
        wyrazy(ii) = append(wyrazy(ii),'I');
    end
    for jj = 1:count_L(ii)
        wyrazy(ii) = append(wyrazy(ii),'L');
    end
    for jj = 1:count_K(ii)
        wyrazy(ii) = append(wyrazy(ii),'K');
    end
    for jj = 1:count_M(ii)
        wyrazy(ii) = append(wyrazy(ii),'M');
    end
    for jj = 1:count_F(ii)
        wyrazy(ii) = append(wyrazy(ii),'F');
    end
    for jj = 1:count_P(ii)
        wyrazy(ii) = append(wyrazy(ii),'P');
    end
    for jj = 1:count_S(ii)
        wyrazy(ii) = append(wyrazy(ii),'S');
    end
    for jj = 1:count_T(ii)
        wyrazy(ii) = append(wyrazy(ii),'T');
    end
    for jj = 1:count_W(ii)
        wyrazy(ii) = append(wyrazy(ii),'W');
    end
    for jj = 1:count_Y(ii)
        wyrazy(ii) = append(wyrazy(ii),'Y');
    end
    for jj = 1:count_V(ii)
        wyrazy(ii) = append(wyrazy(ii),'V');
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

seqType = 'AA';
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
end 

