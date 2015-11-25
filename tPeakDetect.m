function [stg1succ, stg2succ, stg3succ, SuccInd, FosEst] = tPeakDetect(SynCorr, abs_SynCorr,...
    data, index, ComShPNseq, LgPNseq, pkthre, LongSeqThres)

%VarDB = readVariables;
% pnSeqLen = VarDB.pnSeqLen; 
% short_seq_rep = VarDB.short_seq_rep; % NOTE: keep it even number
% longPnSeqLen = VarDB.longPnSeqLen;
%data_len = VarDB.data_len;
%ShortRepThres = VarDB.ShortRepThres;
noofsamples = length(data);
%LongSeqThres = VarDB.LongSeqThres;
% Length of short sequence used for synchronization
global pnSeqLen;

% Number of repetition of short sequence
global short_seq_rep; 

% Length of long sequence
global longPnSeqLen;
global ShortRepThres;
global data_len;
global orgdata;
global datablockcount;
global DataPNseq;
global dataSeqLen;
global BERinfo;
global next_expected_peak;
global correctpkcounter;
% if short_seq_rep<4 || rem(short_seq_rep, 2)>0
%     error('short_seq_rep should be >= 4 and an even number');
% end
stg1succ = 0;
stg2succ = 0;
stg3succ = 0;
SuccInd = 0;
avg_freq_offset = -100;
avg_freq_offset2 = -100;
thre_factor = 0.9;
stg1succ_count = 0;
FosEst = 0;
indexoffset = 45433;
pilot_ind = 8;

if index > (short_seq_rep-1)*pnSeqLen + data_len
    stg1succ =1;
    for rep = 1:short_seq_rep
        if abs_SynCorr(index-(rep-1)*pnSeqLen) > pkthre
            stg1succ_count = stg1succ_count+1;            
        end
        tempSParr(rep) = abs_SynCorr(index-(rep-1)*pnSeqLen);
    end
    %{
    if mod(index-indexoffset, 1920) == 0
        disp(num2str([stg1succ_count,  sort(tempSParr), pkthre]));
    end
    %}
    if stg1succ_count < ceil(ShortRepThres) % approx frc% of the peaks are detected
        stg1succ = 0;
    %elseif sum(PkInfo) < short_seq_rep
    %    stg1succ = 0;
    end
    %%%%%%%%%%
    i_max=index;
    %%%%%%%%%%
    if stg1succ == 1   %  stage 1 successful !!!   
        %msg = ['Stage 1 successful ', num2str(index)];
        %        disp(msg); 
        stg2succ = 1;
        for rep = 1:short_seq_rep-1
            pk_phase_offset = ...
                SynCorr(index-(rep-1)*pnSeqLen) / SynCorr(index-rep*pnSeqLen);
            if SynCorr(index-rep*pnSeqLen) == 0
                pk_phase_offset = 0;
            end
            pk_freq_offset(rep) = angle(pk_phase_offset)/(2*pi*pnSeqLen);
        end
        avg_freq_offset = mean(pk_freq_offset);
        numerator_corr = sum(SynCorr(index : -pnSeqLen : index-(ceil(short_seq_rep/2)-1)*pnSeqLen));
        denominator_corr = sum(SynCorr(index-ceil(short_seq_rep/2)*pnSeqLen:...
                -pnSeqLen : index-(short_seq_rep-1)*pnSeqLen));
        pk_phase_offset = numerator_corr/denominator_corr;
        pk_freq_offset2 = angle(pk_phase_offset)/(short_seq_rep*pi*pnSeqLen);
        avg_freq_offset2 = mean(pk_freq_offset2);


        if stg2succ == 1 % Now stage 2 is successful!!! Start long sequence correlation
            stg3succ = 0; 
            
            i_os = index+pnSeqLen;
            i_os_end = index+pnSeqLen+longPnSeqLen-1;%size(rxdata,2);%i_os+2*pnSeqLen-1;

            f_os_corr = avg_freq_offset;%avg_freq_offset
            stg2succ = f_os_corr;
            %f_os_corr = 0;
            
            
            tempdata = data(i_os:i_os_end)...
                .*exp(-1j*2*pi*f_os_corr*(0:longPnSeqLen-1));
            
            %LongCorr = tCorrelate(LgPNseq', data(index+pnSeqLen:index+pnSeqLen+longPnSeqLen-1));
            LongCorr = tCorrelate(LgPNseq', tempdata);
            %%%%%%%%%%%% ORBIT %%%%%%%%%%%%%%
            %LongSeqThres = 60;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            stg3succ = abs(LongCorr);
            if stg3succ > LongSeqThres
                correctpkcounter = correctpkcounter+1;
                msg = [num2str(correctpkcounter),'### ', num2str(index), ': ',...
                num2str(index+1920), ': ', num2str(f_os_corr),...
                ': ', num2str(abs(LongCorr)), ': ', num2str(stg1succ_count),...
                ': ',num2str(LongSeqThres), ': ', num2str(mean(abs(orgdata(index-100:index)))*320/sqrt(12)),...
                ];
                disp(msg);
                if mod(index-indexoffset, 1920) ~= 0
                    disp('***False***');
                end
                next_expected_peak = index+1920;
                 %############## finding the data bits #####################
                phaseoffset = 0;
                bitinfo = [];
                bitinfobu = [];
                angbit = [];
                datablockcount = datablockcount+1;
                %hold on;
                %%%%%%%%%%% using phase average %%%%%%%%
%                 for inloop=1:data_len/dataSeqLen          
%                     i_os = index+pnSeqLen+longPnSeqLen+(inloop-1)*dataSeqLen;
%                     i_os_end = index+pnSeqLen+longPnSeqLen+inloop*dataSeqLen-1;%size(rxdata,2);%i_os+2*pnSeqLen-1;
%                     if i_os_end>noofsamples
%                         break;
%                     end
%                     tempdata = data(i_os:i_os_end)...
%                         .*exp(-1j*2*pi*f_os_corr*((inloop-1)*dataSeqLen:inloop*dataSeqLen-1));
%                     bitinfobu(inloop) = (tCorrelate(DataPNseq', tempdata));
%                     angbit(inloop) = angle(bitinfobu(inloop));
% %                     if(angbit(inloop)<0)
% %                         angbit(inloop) = angbit(inloop) + 2*pi;
% %                     end
%                     temp_offset(inloop) = 0; 
%                     if rem(inloop-1,pilot_ind)==0
%                         % following line assumes that +1 was sent as pilot.
%                         temp_offset(inloop) = bitinfobu(inloop)/(abs(bitinfobu(inloop)));%mean(angbit);
%                         % following line assumes that -1 was sent as pilot.
%                         %phaseoffset = -pi+angbit(inloop);%mean(angbit);
%                     end
%                 end
%                 tempX = sum(temp_offset);
%                 phaseoffset = angle(tempX);%/(data_len/(dataSeqLen*pilot_ind));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for inloop=1:data_len/dataSeqLen            
                    i_os = index+pnSeqLen+longPnSeqLen+(inloop-1)*dataSeqLen;
                    i_os_end = index+pnSeqLen+longPnSeqLen+inloop*dataSeqLen-1;%size(rxdata,2);%i_os+2*pnSeqLen-1;
                    if i_os_end>noofsamples
                        break;
                    end
                    tempdata = data(i_os:i_os_end)...
                        .*exp(-1j*2*pi*f_os_corr*((inloop-1)*dataSeqLen:inloop*dataSeqLen-1));
                    bitinfobu(inloop) = (tCorrelate(DataPNseq', tempdata));
                    angbit(inloop) = angle(bitinfobu(inloop));
                    if(angbit(inloop)<0)
                        angbit(inloop) = angbit(inloop) + 2*pi;
                    end
                    if rem(inloop-1,pilot_ind)==0
                        % following line assumes that +1 was sent as pilot.
                        phaseoffset = angbit(inloop);%mean(angbit);
                        % following line assumes that -1 was sent as pilot.
                        %phaseoffset = -pi+angbit(inloop);%mean(angbit);
                    end
                    bitinfo(inloop) = bitinfobu(inloop)*exp(-1j*phaseoffset);
                    %plot(15*bitinfobu./(abs(bitinfobu)),'bo');axis([-30 30 -30 30]); grid on;
                    %plot(15*bitinfo./(abs(bitinfo)),'bo');axis([-30 30 -30 30]); grid on;
                    %pause(0.1);
                    if rem(inloop-1,pilot_ind)==0 % #This 'if' is not needed, right now calc BER
                        bitval(inloop) = 0;
                    else
                        if bitinfo(inloop)>0
                            bitval(inloop) = 1;
                        else
                            bitval(inloop) = 0;%-1; % ###switch it to -1 when doing comprehansive
                        end
                    end
                end
                %hold off;
                %close all;
                if i_os_end<noofsamples
                    BitErrorRate(datablockcount) = 1-(sum(bitval)/floor((data_len/dataSeqLen)*(1-1/pilot_ind)));
                    BERinfo(1,datablockcount) = index;
                    BERinfo(2,datablockcount) = BitErrorRate(datablockcount);
                end
                %##########################################################
            elseif index == next_expected_peak
                disp('*******MISSED THE PEAK*******');
            end
        end
    end
end
if stg3succ == 0
    if rem(index+pnSeqLen+longPnSeqLen, data_len + short_seq_rep * pnSeqLen+longPnSeqLen) == 1
        msg = ['### correct peak missed ', num2str(index), ': ',...
            num2str(stg1succ), ': ', num2str(stg3succ), ': ',...
            num2str(avg_freq_offset), ': ', num2str(avg_freq_offset2)];
        %disp(msg); 
        SuccInd = 2;
    end
end