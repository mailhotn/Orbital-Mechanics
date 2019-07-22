function [state,options,optchanged] = GaMoGenSave(options,state,~,OptParams,latEm)
save([OptParams.datafolder '\MoGaTemp_Lat_' num2str(latEm) '.mat'],'state');
optchanged = false;
end