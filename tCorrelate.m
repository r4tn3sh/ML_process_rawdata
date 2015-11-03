function [tCorr] = tCorrelate(pnSeq, data)

if length(pnSeq) == length(data)
    r_sh_sync_corr = sum(pnSeq.*real(data));%dot(pnSeq,real(data));
    i_sh_sync_corr = sum(pnSeq.*imag(data));%dot(pnSeq,imag(data));
    tCorr = complex(r_sh_sync_corr, i_sh_sync_corr);
else
    error('Length of PN sequence and data should be same');
end