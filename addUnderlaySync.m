function [rxdata, gold_seq] = addUnderlaySync(fft_len, cp_len, txdata, sync_len, P_sync)
%% Prepare PN sequence
% Find primitive polynomial
p_poly_size = ceil(log2(sync_len));
p_poly_all = primpoly(p_poly_size,'all','nodisplay');
p_poly_1 = fliplr(de2bi(p_poly_all(1)));
p_poly_2 = fliplr(de2bi(p_poly_all(2)));
init_cond_1 = p_poly_2;
init_cond_1(1) = [];
init_cond_2 = p_poly_1;
init_cond_2(1) = [];

DiffCode=0;

if (DiffCode == 0)
    overall_len = sync_len;
else
    overall_len = sync_len+1;
end
% Generating gold sequence
gold_seq_poly = comm.GoldSequence('FirstPolynomial',p_poly_1,...
             'SecondPolynomial', p_poly_2,...
             'FirstInitialConditions', init_cond_1,...
             'SecondInitialConditions', init_cond_2,...
             'Index', 4, 'SamplesPerFrame', overall_len);
if (DiffCode == 0)  
    real_gs = (2*step(gold_seq_poly))-1;
    imag_gs = zeros(overall_len,1);
    %gold_seq = complex(real_gs,imag_gs);
    gold_seq = real_gs;
else    
    real_gs = (2*step(gold_seq_poly))-1;
    imag_gs = zeros(overall_len,1);
    %gold_seq = complex(real_gs,imag_gs);
    gold_seq = real_gs;
    shifted_gold_seq = [0; gold_seq(1:sync_len)];
    %diff_gold_seq = (gold_seq - shifted_gold_seq);
    %diff_gold_seq = xor(gold_seq, shifted_gold_seq);
    diff_gold_seq(1) = gold_seq(1);
    for i=2:sync_len
        diff_gold_seq(i) = gold_seq(i)/diff_gold_seq(i-1);
    end
end
%% Add the PN as underlay
for i=1:size(txdata,2)
    if (DiffCode == 0)
        rxdata(i) = txdata(i)+sqrt(P_sync)*gold_seq(rem(i-1,sync_len)+1);
    else
        rxdata(i) = txdata(i)+sqrt(P_sync)*diff_gold_seq(rem(i-1,sync_len)+1);
    end
end
