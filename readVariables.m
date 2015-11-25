function [data] = readVariables()

global ShortPkThres;
% Length of short sequence used for synchronization
global pnSeqLen;

% Number of repetition of short sequence
global short_seq_rep; 

% Length of long sequence
global longPnSeqLen;
global f_os;
global ShortRepThres;
global data_len;
% FFT length
fft_len=64;
% CP length
cp_len = fft_len/4;
total_len = fft_len+cp_len;
% Overall number of bits to be sent
noofsamples = 50000;
noofsamples = noofsamples - rem(noofsamples,fft_len);

% Desired signal to PN-sequence power ratio
desired_SPR = 10; % in dB

% Number of blocks of data between two sync blocks
data_len = 2*(pnSeqLen*short_seq_rep + longPnSeqLen);%6*total_len;

% Frequency offset
%f_os=0.001;

% Noise power
P_noise = 1;

% Sync power
P_sync = 1;

%% Thresholds
ShortPkThres = 0.25*pnSeqLen;   % Peak threshold for short correlation
ShortRepThres = ceil(0.5*short_seq_rep); % How many peaks to be detected
LongSeqThres = 0.7*longPnSeqLen; % Peak threshold for long correlation

data = struct('fft_len', fft_len,...
                  'cp_len', cp_len,...
                  'total_len', total_len,...
                  'noofsamples', noofsamples,...
                  'desired_SPR', desired_SPR,...
                  'pnSeqLen', pnSeqLen,...
                  'short_seq_rep', short_seq_rep,...
                  'longPnSeqLen', longPnSeqLen,...
                  'data_len', data_len,...
                  'f_os', f_os,...
                  'P_noise', P_noise,...
                  'P_sync', P_sync,...
                  'ShortPkThres', ShortPkThres,...
                  'ShortRepThres', ShortRepThres,...
                  'LongSeqThres', LongSeqThres...
                  );