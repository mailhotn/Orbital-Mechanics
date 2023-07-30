function reportIFTTT(GoogleDrivePath,time)
if nargin < 2
    time = pi;
end

c = clock;
fileID = fopen([GoogleDrivePath '\Doc Data\zNotifications\CalcDone_'...
    num2str(c(3)) '-' num2str(c(2)) '_' num2str(c(4)) ...
    '.' num2str(c(5)) '.txt'],'w'); % create error file for IFTT phone notification
fprintf(fileID,['Elapsed time: ' num2str(time) newline...
    '☜(⌒▽⌒)☞ ']);