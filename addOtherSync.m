function [rxdata] = addOtherSync(txdata)

VarDB = readVariables;
fft_len = VarDB.fft_len;
cp_len = VarDB.cp_len;
% pnSeqLen = VarDB.pnSeqLen; 
P_sync = VarDB.P_sync;
P_noise = VarDB.P_noise;
f_os = VarDB.f_os;

% Length of short sequence used for synchronization
global pnSeqLen;

%% Prepare PN sequence
% Find primitive polynomial
p_poly_size = ceil(log2(pnSeqLen));
p_poly_all = primpoly(p_poly_size,'all','nodisplay');
p_poly_1 = fliplr(de2bi(p_poly_all(3)));
p_poly_2 = fliplr(de2bi(p_poly_all(4)));
init_cond_1 = round(rand(1,p_poly_size));
init_cond_2 = round(rand(1,p_poly_size));

overall_len = pnSeqLen;
    
%% Generating Gold sequence
pn_seq_poly = comm.GoldSequence('FirstPolynomial',p_poly_1,...
             'SecondPolynomial', p_poly_2,...
             'FirstInitialConditions', init_cond_1,...
             'SecondInitialConditions', init_cond_2,...
             'Index', 4, 'SamplesPerFrame', overall_len);
real_gs = (2*step(pn_seq_poly))-1;

%%  Generating short PN sequence
% pn_seq_poly = comm.PNSequence('Polynomial',p_poly_1,...
%              'InitialConditions', init_cond_1,...
%              'SamplesPerFrame', overall_len);
% real_gs = (2*step(pn_seq_poly))-1;
%% Generating PN sequence
% pn_seq_poly = comm.PNSequence('Polynomial',p_poly_1,...
%              'InitialConditions', init_cond_1,...
%              'SamplesPerFrame', overall_len);
%real_gs = (2*step(pn_seq_poly))-1;
%% Generating Hadamard sequence
% pn_seq_poly = comm.HadamardCode('Length', 2^p_poly_size,'Index', randi(overall_len),'SamplesPerFrame', overall_len);
% real_gs = step(pn_seq_poly);


pn_seq = real_gs;
noofsamples = length(txdata);
%% Add the PN as underlay and additive noise
temp_i = randi(noofsamples);
for i=1:noofsamples
    rxdata(rem(temp_i+i, noofsamples)+1) = txdata(rem(temp_i+i, noofsamples)+1)+0*...
        pn_seq(rem(i-1,pnSeqLen)+1);%+P_noise*(randn + 1i*randn)/sqrt(2);
end

%rxdata=rxdata.*exp(2*pi*f_os*1j.*(1:length(rxdata)));
