function [txdata, data] = get_txdata()

VarDB = readVariables;
fft_len = VarDB.fft_len;
cp_len = VarDB.cp_len;
noofsamples = VarDB.noofsamples;
desired_SPR = VarDB.desired_SPR; % in dB
%% Prepare data
    real_data = 2*round(rand(1,noofsamples))-1;
    img_data = 2*round(rand(1,noofsamples))-1;
    data = complex(real_data, img_data)/sqrt(2);
    % normalized power data
    norm_txdata = [];
    for i=1:(noofsamples/fft_len)
        norm_txdata = cat(2, norm_txdata, (1/sqrt(fft_len))*fft(data((i-1)*fft_len+1:i*fft_len)));
    end
    % HERE CP SHOULD BE ADDED %
    % Data power
    P_data = 10^(desired_SPR/10);
    
    txdata = sqrt(P_data)*(norm_txdata);