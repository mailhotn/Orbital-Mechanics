function errorIFTTT(dbPath,message)
if nargin < 2
    message = ':(';
end

c = clock;
fileID = fopen([dbPath '\zzzMatlab Errors\Error_'...
    num2str(c(3)) '-' num2str(c(2)) '_' num2str(c(4)) ...
    '.' num2str(c(5)) '.txt'],'w'); % create error file for IFTT phone notification
fprintf(fileID,message);