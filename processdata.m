clc;
clear;
close all;

% PARAMETERS NEEDED TO BE ADJUSTED
totallen=1280+240;
load('path to/longdata.mat')
load('path to/lgthdata.mat')
load('path to/corrdata.mat')


datalen = length(corrdata);
dataoffset = ceil(0.3*datalen); 
totalpeaks = sum(corrdata(dataoffset:datalen));
previndex=0;
currindex=0;
initfalsecnt=0;
% Find the location of potention correct peaks and discard the wrong
for i=totallen+dataoffset:datalen
    pflag(i)=0;
    if corrdata(i)==1
        for j=i-totallen:-totallen:1
            if corrdata(j)==1
                pflag(i)=1;
                break;
            end
        end
        if pflag(i)==0
            initfalsecnt = initfalsecnt+1;
        end
    end
end
doublepeakcnt=0;
missedpeakcnt=0;
for i=dataoffset:datalen
    indexdiff(i)=0;
    flag(i)=0;
    flagval(i)=0;
    correctpeakflag(i) = 0;
    if pflag(i)==1
        currindex = i;
        if previndex>0
            indexdiff(i) = currindex-previndex;
            if rem(indexdiff(i),totallen)>1 && rem(indexdiff(i),totallen)<totallen-1
                flag(i)=1;
                flagval(i)=indexdiff(i);
            elseif i>1 && rem(indexdiff(i-1),totallen)==totallen-1 && rem(indexdiff(i),totallen)~=1
                flag(i)=1;
                flagval(i)=indexdiff(i);
            elseif i>1 && rem(indexdiff(i-1),totallen)~=totallen-1 &&...
                    rem(indexdiff(i-1),totallen)>0 && rem(indexdiff(i),totallen)==1
                flag(i)=1;
                flagval(i)=indexdiff(i);
            elseif i>1 && rem(indexdiff(i-1),totallen)>0 && rem(indexdiff(i),totallen)==1
                correctpeakflag(i) = 1;
                doublepeakcnt=doublepeakcnt+1;
                missedpeakcnt = missedpeakcnt+(indexdiff(i-1)/totallen)-1;
            elseif rem(indexdiff(i),totallen)==0
                correctpeakflag(i) = 1;
                missedpeakcnt = missedpeakcnt+(indexdiff(i)/totallen)-1;
            else
                a=0;
            end
%             if i>1 && rem(indexdiff(i),totallen)==1 &&
%                 flag(i)=1;
%                 flagval(i)=indexdiff(i);
%             end
        end
    end
    previndex=currindex;
end
falsepeakcnt=sum(flag)+initfalsecnt;
correctpeakcnt = sum(correctpeakflag);
maxcorr = 1.2*max(abs(longdata));
hold on;
plot(maxcorr*flag,'r');
plot(maxcorr*correctpeakflag,'g');
plot(abs(longdata),'b');
plot(lgthdata,'c');

str = sprintf('False alarm rate: %f',falsepeakcnt/(totalpeaks-doublepeakcnt));
disp(str);
str = sprintf('Missed rate: %f',missedpeakcnt/((datalen-dataoffset)/totallen));
disp(str);
