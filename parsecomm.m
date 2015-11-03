clc, clear all, close all
f = fopen('commoutput1.txt','rt');  % Read Input File
C = textscan(f, '%s', 'Delimiter', ' \r\n,');
C = C{1}; % Store Input File in a Cell
fclose(f);

startIndx = find(~cellfun(@isempty, regexp(C, '************************', 'match')));

%assert(all(size(startIndx) == size(endIndx)))
loop = 0;
loop1 = 0;
loop2 = 0;
for rep=1:length(startIndx)   
    loop = loop+1;
    if rem(rep, 5) == 1
        loop = 1;
        loop1 = loop1+1;
    end
    if rem(rep, 25) == 1
        loop1 = 1;
        loop2 = loop2+1;
    end
    SuccRate(loop, loop1, loop2) = str2double(C(startIndx(rep)+3));
    FalseAlarm(loop, loop1, loop2) = str2double(C(startIndx(rep)+9));
    Tpks = str2double(C(startIndx(rep)-11));
    Fpks = str2double(C(startIndx(rep)-6));
    PFalseAlarm(loop, loop1, loop2) = Fpks/(Fpks+Tpks);
    MissedDetection(loop, loop1, loop2) = str2double(C(startIndx(rep)+15));
end