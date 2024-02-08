clear
%% Load Data
dataFolder = 'C:\Users\User\Google Drive\Doc Data\Error Mapping'; % ASRI
% dataFolder = 'D:\Dropbox\Doc Fourier Data\Error Mapping'; % Laptop???
% load([dataFolder '\ErrMaps_15-6-2022_9-52.mat']); % First run e:0.005-0.55, i:0.4-90
% load([dataFolder '\ErrMaps_15-6-2022_11-13.mat']); % Singularity test e:0.5-0.7, i:60-70
% load([dataFolder '\ErrMaps_21-6-2022_4-33.mat']); % Big Mapping e:0.001-0.7, i:0.1-179.9

% Added Deprit
% load([dataFolder '\ErrMaps_12-8-2022_1-22.mat']); % First Deprit - total garbage due to error
% load([dataFolder '\ErrMaps_16-8-2022_0-53.mat']); % First functional Deprit
% load([dataFolder '\ErrMaps_19-8-2022_19-57.mat']); % Deprit w/ new elliptics
% load([dataFolder '\ErrMaps_22-8-2022_15-3.mat']); % Deprit w/ new elliptics - more thorough
% load([dataFolder '\ErrMaps_1-9-2022_13-9.mat']); % Deprit & SP
% load([dataFolder '\ErrMaps_5-9-2022_13-45.mat']); % SP + Deprit + k5

% Added mean sma to all of Fourier
% load([dataFolder '\ErrMaps_6-9-2022_13-30.mat']); % SP + Deprit + k4
% load([dataFolder '\ErrMaps_23-9-2022_14-27.mat']); % Cancellations, SP + Deprit + k4

% Use mean elements for all IC
% load([dataFolder '\ErrMaps_28-10-2022_12-17.mat']); % Cancellations, SP + Deprit + k4SP

% Variation of AOP
% load([dataFolder '\ErrMaps_12-11-2022_1-54.mat']); % SP + Deprit + 2Ord
% load([dataFolder '\ErrMaps_16-11-2022_1-35.mat']); % SP + Deprit + 2Ord w/o dAop

% Wrapping SP change
% load([dataFolder '\ErrMaps_18-11-2022_2-6.mat']); % SP + Deprit + 2Ord
% load([dataFolder '\ErrMaps_23-11-2022_2-45.mat']); % SP + Deprit + 2Ord w/o dAop
% load([dataFolder '\ErrMaps_3-12-2022_4-37.mat']); % SP + Deprit + 2Ord w/o dAop Low ECC
% load([dataFolder '\ErrMaps_20-12-2022_9-33.mat']); % SP + Deprit + 2Ord w/o dAop Lower left
% load([dataFolder '\ErrMaps_27-12-2022_10-19.mat']); % SP + Deprit + 2OrdK5 w/o dAop Low ecc 
% load([dataFolder '\ErrMaps_5-1-2023_7-47.mat']); % 2Ord: K5vK4 w/o dAop Low ecc 
% load([dataFolder '\ErrMaps_10-1-2023_16-35.mat']); % 2Ord: K5vK6 w/o dAop Low ecc
% load([dataFolder '\ErrMaps_15-3-2023_21-25.mat']); % 2Ord: K5 w/o dAop Low ecc, no Inc singularity in OeOsc 

% Change to Propagator - Normalized + ode78 + tolerances
% load([dataFolder '\ErrMaps_24-5-2023_20-48.mat']); % 2Ord: K5 w/o dAop Low ecc, no Inc singularity in OeOsc 
% load([dataFolder '\ErrMaps_31-5-2023_23-54.mat']); % same as last + use mean aop
% load([dataFolder '\ErrMaps_27-12-2023_11-45.mat']); % same as last + variable initial M, small sample
% load([dataFolder '\ErrMaps_27-12-2023_12-30.mat']); % same as last + variable M, small sample 
% load([dataFolder '\ErrMaps_27-12-2023_13-52.mat']); % variable initial M, removed Brouwer weird wrapping small
% load([dataFolder '\ErrMaps_28-12-2023_21-59.mat']); % same as last large sample
load([dataFolder '\ErrMaps_24-1-2024_3-57.mat']); % Low sma




incRange = MapData.incRange;
eccRange = MapData.eccRange;
nEcc = length(eccRange);
nInc = length(incRange);
errTenF = MapData.errTenF;
errTenF2 = MapData.errTenF2;
errTenB = MapData.errTenB;
errTenD = MapData.errTenD;

errTenB(errTenB==inf) = nan;

%% Cut off low ecc
% eccCut = 0.1;
% eccRange = eccRange(eccRange>eccCut);
% errTenF = errTenF((nEcc-length(eccRange)+1):end,:,:);
% errTenB = errTenB((nEcc-length(eccRange)+1):end,:,:);

%% Cut off high inc - shows that brouwer always had issues
% incCut = 50;
% incRange = incRange(incRange<incCut);
% errTenF = errTenF(:,1:length(incRange),:);
% errTenB = errTenB(:,1:length(incRange),:);
% errTenD = errTenD(:,1:length(incRange),:);
%% Plot Relative Errors
plotPropErrMap(eccRange,incRange,errTenF,errTenB,'Fourier','Brouwer SP',0)

% plotPropErrMap(eccRange,incRange,errTenF2,errTenB,'Fourier2','Brouwer SP',1)

% plotPropErrMap(eccRange,incRange,errTenF2,errTenF,'Fourier2','Fourier',2)

%% Plot Absolute Errors
% plotPropErrMap(eccRange,incRange,errTenF2,[],'Fourier2',[],10)

plotPropErrMap(eccRange,incRange,errTenF,[],'Fourier',[],20)

plotPropErrMap(eccRange,incRange,errTenB,[],'Brouwer',[],30)