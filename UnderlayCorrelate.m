function [rxdata, sync_corr] = UnderlayCorrelate(rxdata, sync_len, gold_seq, pk_threshold)
f_os_corr = 0;
pk_freq_offset = 0;
past_peak = 0;
current_peak = 0;
ind_offset = 0;
rxdatabu = rxdata;
for i=1:(size(rxdata,2)-sync_len)
    r_sync_corr(i) = sum(transpose(real(gold_seq)).*real(rxdata(i:i+sync_len-1)));
    i_sync_corr(i) = sum(transpose(real(gold_seq)).*imag(rxdata(i:i+sync_len-1)));
    sync_corr(i) = r_sync_corr(i)^2 + i_sync_corr(i)^2;
    %%%%%%%%%%%%%% Freq offset calculation and correction %%%%%%%%%%%%%
    %%{
    if i>=sync_len
        if sync_corr(i) > pk_threshold
            current_peak = i;
            if current_peak == past_peak+sync_len % found two consecutive peaks
                pk_phase_offset = sync_corr(current_peak) / sync_corr(past_peak);
                pk_freq_offset = angle(pk_phase_offset)/(4*pi*sync_len);
                ind_offset = i+1;
                ind_offset_end = noofsamples;%ind_offset+2*sync_len-1;
                f_os_corr = pk_freq_offset;
                rxdata(ind_offset:ind_offset_end) = rxdata(ind_offset:ind_offset_end)...
                    .*exp(-1j*2*pi*f_os_corr*(0 : ind_offset_end-ind_offset));
            end
            past_peak = current_peak;
        end
    end
    %}
end