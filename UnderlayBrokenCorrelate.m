function [rxdata, final_short_SynCorr] = UnderlayBrokenCorrelate...
    (rxdata, orgdata, ShPNseq, LgPNseq)
VarDB = readVariables;
fft_len = VarDB.fft_len;
% pnSeqLen = VarDB.pnSeqLen; 
% short_seq_rep = VarDB.short_seq_rep; % NOTE: keep it even number
% longPnSeqLen = VarDB.longPnSeqLen;
data_len = VarDB.data_len;
%%% Threshold for short sequence
ShortPkThres = VarDB.ShortPkThres;
ShortRepThres = VarDB.ShortRepThres;
   
% Length of short sequence used for synchronization
global pnSeqLen;

% Number of repetition of short sequence
global short_seq_rep; 

% Length of long sequence
global longPnSeqLen;   
global LongCorr;
global LongThreshold;
past_peak = 0;
current_peak = 0;
i_os = 0;
avgcalcstart = 100;
ComShPNseq = [];
for index=1:short_seq_rep
    ComShPNseq = [ComShPNseq ShPNseq(index,:)];
end

pkind = 0;
pkCount = zeros(2,length(rxdata));
totalPkCount=0;
CorrectPkCount=0;
MissedPkCount=0;
IncorrPkCount=0;
qpskConstell = (1/sqrt(2))*[1+1i 1-1i -1+1i -1-1i];
    %%%%%%%%%%%% ORBIT %%%%%%%%%%%%%%
    ShortPkThres = 1;%0.035;%10
    LongSeqThres = 0;
    LongThresFlag = 0;
    LgThresUpdateTime = 1920;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for index=1:(length(rxdata) - short_seq_rep*pnSeqLen - longPnSeqLen)  
    for rep=1:short_seq_rep
        short_SynCorr(rep,index) = tCorrelate(ShPNseq(rep,:), rxdata(index:index+pnSeqLen-1));
    end
    % Add up the sync to get better peaks
    final_short_SynCorr(index) = 0;
    for rep=1:short_seq_rep
        if (index-(rep-1)*pnSeqLen <= 0)
            tempval = 0;
        else
            tempval = abs(short_SynCorr(rep,index-(rep-1)*pnSeqLen));
        end
        final_short_SynCorr(index) = final_short_SynCorr(index) + tempval;
    end
        
    % Normalize the peaks
    if index<=avgcalcstart
        norm_short_SynCorr(index) = final_short_SynCorr(index)/mean(final_short_SynCorr(1:index));
    else
        norm_short_SynCorr(index) = final_short_SynCorr(index)/mean(final_short_SynCorr(index-avgcalcstart:index));
    end
    %% 
    % Freq offset calculation/correction ; 
    % NOTE: this is not a generic method. only for same pn seq repetition
    SynCorr(index) = short_SynCorr(1,index);
    abs_SynCorr(index) = abs(SynCorr(index));
    %%%% Long threshold init
    if index<longPnSeqLen
        tempThres = mean(abs(rxdata(1:index)))*longPnSeqLen/sqrt(12);
        LongSeqThres = 0.7*LongSeqThres + 0.3*tempThres*0.8;
    end
    if index<pnSeqLen+1
        [flag, flag1, flag2, SuccInd, FosEst] = tPeakDetect(SynCorr, abs_SynCorr, rxdata, index,...
            ComShPNseq, LgPNseq, ShortPkThres, LongSeqThres);
    else
        [flag, flag1, flag2, SuccInd, FosEst] = tPeakDetect(SynCorr, abs_SynCorr, rxdata, index,...
            ComShPNseq, LgPNseq, ShortPkThres, LongSeqThres);
    end
    LongCorr(index) = flag2;
    FreqOS(index) = flag1;
    %%%% long threshold update
    if flag2>0
        LongThresFlag = flag2;
        for rep = 1:short_seq_rep
            tempSParr(rep) = abs_SynCorr(index-(rep-1)*pnSeqLen);
        end
        %disp(mean(abs(orgdata(index-100:index)))*320/sqrt(12));
        sortThres = sort(tempSParr);
        newThres = sortThres((short_seq_rep-ShortRepThres));
        ShortPkThres = 0.7*ShortPkThres +0.3*newThres*0.9;
    end
    if rem(index, LgThresUpdateTime) == 0
        if LongThresFlag>0            
            tempThres = LongThresFlag;
        else
%             newThres = mean(abs(rxdata(index-longPnSeqLen:index)))*pnSeqLen/sqrt(12);
%             ShortPkThres = 0.7*ShortPkThres +0.3*newThres*0.9;
            tempThres = mean(abs(rxdata(index-longPnSeqLen:index)))*longPnSeqLen/sqrt(12);
            ShortPkThres = 0.7*ShortPkThres +0.3*1;
        end
        LongSeqThres = 0.7*LongSeqThres + 0.3*tempThres*0.8;
        LongThresFlag = 0;
        disp(LongSeqThres);
    end
    LongThreshold(index) = LongSeqThres;
end