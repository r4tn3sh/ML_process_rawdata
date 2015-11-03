function [data, ShPNseq, LgPNseq] = addBrokenUnderlaySync(txdata)

VarDB = readVariables;
% pnSeqLen = VarDB.pnSeqLen; 
% short_seq_rep = VarDB.short_seq_rep; % NOTE: keep it even number
% longPnSeqLen = VarDB.longPnSeqLen;
data_len = VarDB.data_len;
P_sync = VarDB.P_sync;
                
% Length of short sequence used for synchronization
global pnSeqLen;

% Number of repetition of short sequence
global short_seq_rep; 

% Length of long sequence
global longPnSeqLen;    

%% Add short PN sequences
p_poly_size = ceil(log2(pnSeqLen));
p_poly_all = primpoly(p_poly_size,'all','nodisplay');

init_cond_1 = round(rand(1,p_poly_size));
init_cond_2 = round(rand(1,p_poly_size));
p_poly_1 = fliplr(de2bi(p_poly_all(1)));
p_poly_2 = fliplr(de2bi(p_poly_all(2)));
% using different PN sequence in each block, set diffPN to 1
diffPN = 0;
if short_seq_rep>0
    for i=1:short_seq_rep
        if diffPN == 1
            p_poly_1 = fliplr(de2bi(p_poly_all(2*i-1)));
            p_poly_2 = fliplr(de2bi(p_poly_all(2*i)));
            init_cond_1 = round(rand(1,p_poly_size));
            init_cond_2 = round(rand(1,p_poly_size));
        end
        overall_len = pnSeqLen;
    %%  Generating short Gold sequence
%         pn_seq_poly = comm.GoldSequence('FirstPolynomial',p_poly_1,...
%                      'SecondPolynomial', p_poly_2,...
%                      'FirstInitialConditions', init_cond_1,...
%                      'SecondInitialConditions', init_cond_2,...
%                      'Index', 4, 'SamplesPerFrame', overall_len);
%         real_gs = (2*step(pn_seq_poly))-1;
    %%  Generating short PN sequence
        pn_seq_poly = comm.PNSequence('Polynomial',p_poly_1,...
                     'InitialConditions', init_cond_1,...
                     'SamplesPerFrame', overall_len);
        real_gs = (2*step(pn_seq_poly))-1;
    %%  Generating short hadamard sequence
    %     pn_seq_poly = comm.HadamardCode('Length', 2^p_poly_size,'Index', randi(overall_len),'SamplesPerFrame', overall_len);
    %     real_gs = step(pn_seq_poly);


        ShPNseq (i,:) = real_gs;
    end
end

%% Add long PN sequences
if longPnSeqLen>0
    p_poly_size = ceil(log2(longPnSeqLen));
    p_poly_all = primpoly(p_poly_size,'all','nodisplay');

    init_cond_1 = round(rand(1,p_poly_size));
    init_cond_2 = round(rand(1,p_poly_size));
    p_poly_1 = fliplr(de2bi(p_poly_all(1)));
    p_poly_2 = fliplr(de2bi(p_poly_all(2)));
    overall_len = longPnSeqLen;
    % Generating long pn sequence
    gold_seq_poly = comm.GoldSequence('FirstPolynomial',p_poly_1,...
                 'SecondPolynomial', p_poly_2,...
                 'FirstInitialConditions', init_cond_1,...
                 'SecondInitialConditions', init_cond_2,...
                 'Index', 4, 'SamplesPerFrame', overall_len);
    real_gs = (2*step(gold_seq_poly))-1;

    %%  Generating long PN sequence
%         pn_seq_poly = comm.PNSequence('Polynomial',p_poly_1,...
%                      'InitialConditions', init_cond_1,...
%                      'SamplesPerFrame', overall_len);
%         real_gs = (2*step(pn_seq_poly))-1;


    LgPNseq = real_gs;
else
    LgPNseq = [];
end

%% Add the PN as underlay
for i=1:length(txdata)
    j = rem(i-1,data_len + pnSeqLen*short_seq_rep + longPnSeqLen)+1;
    if j <= data_len
        data(i) = txdata(i);
    elseif j<= data_len + pnSeqLen*short_seq_rep
        rep = ceil((j - data_len)/pnSeqLen);
        if j<= data_len + rep*pnSeqLen
            data(i) = txdata(i)+sqrt(P_sync)*ShPNseq(rep, rem(j-data_len-1, pnSeqLen)+1);
        end
    else
        data(i) = txdata(i)+sqrt(P_sync)*LgPNseq(j-data_len - pnSeqLen*short_seq_rep);
    end
end
