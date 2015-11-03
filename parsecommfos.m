clc, clear all, close all
f = fopen('commoutputnoise.txt','rt');  % Read Input File
C = textscan(f, '%s', 'Delimiter', ' \r\n,');
C = C{1}; % Store Input File in a Cell
fclose(f);

startIndx = find(~cellfun(@isempty, regexp(C, '************************', 'match')));

%assert(all(size(startIndx) == size(endIndx)))
loop = 0;
for rep=1:length(startIndx)   
    loop = loop+1;
    SuccRate(loop) = str2double(C(startIndx(rep)+3));
    FalseAlarm(loop) = str2double(C(startIndx(rep)+9));
    MissedDetection(loop) = str2double(C(startIndx(rep)+15));
    
    fa = str2double(C(startIndx(rep)-6));
    md = str2double(C(startIndx(rep)-1));
    cp = str2double(C(startIndx(rep)-11));
    tp = str2double(C(startIndx(rep)-16));
    PFalseAlarm(loop) = fa/(fa+cp);
    PMissedDetection(loop) = md/tp;
end