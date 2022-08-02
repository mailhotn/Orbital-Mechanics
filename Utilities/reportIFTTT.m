function reportIFTTT(dbPath,time)
if nargin < 2
    time = pi;
end

c = clock;
fileID = fopen([dbPath '\zzzMatlab Notification\OptDone_'...
    num2str(c(3)) '-' num2str(c(2)) '_' num2str(c(4)) ...
    '.' num2str(c(5)) '.txt'],'w'); % create error file for IFTT phone notification
fprintf(fileID,['Elapsed time: ' num2str(time) newline...
    '☜(⌒▽⌒)☞ ']);