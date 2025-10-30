function plotPropErrMapPaper(eccRange,incRange,errTen1,errTen2,label1,label2,saveFolder)
%plotPropErrMap plots error maps for comaprison of propagation methods
% Changes based on first review, same colorbar for K & F in abs, save
% images directly & create total error map
%
% ~~~~~~~~~~~~~~  Inputs  ~~~~~~~~~~~~~~~
% errTen1, errTen2: nEcc x nInc x 6 Tensors
% label1, label2: label strings - optional
% nFig - 10s digit for figure numbering - default 0
% logE - flag to use logarithmic e Axis - no longer input
if nargin < 7
    saveFolder = [];
end

if nargin < 5
    label1 = [];
    label2 = [];
end
if nargin == 3
    errTen2 = [];
end
nFig = 1;
logE = true; % Removed Option to select log scale for e, no reason to ever
% go back to linear, but if yes, just change here instead of input
if logE
    eccRange = log10(eccRange);
    eccLabel = '$\rm{Eccentricity}$'; % switched to labeling ecc directly
    eTicks = flip(eccRange(end):-0.2:eccRange(1)).';
    eTicks(1) = eccRange(1);
    eTickLabels = num2str(10.^eTicks,2);
elseif ~logE
    eccLabel = '$\rm{Eccentricity}$';
end


titleSet = {'r','s','w','\dot{r}','\dot{s}','\dot{w}','d','v'};
unitSet = {'km','km','km','\frac{km}{s}','\frac{km}{s}','\frac{km}{s}','km','\frac{km}{s}'};
fontSizeAll = 14;
fileNameSet = {'Pos-r','Pos-s','Pos-w','Vel-r','Vel-s','Vel-w','Pos-tot','Vel-tot'};
% Add totals to 7 & 8 pos
% Get Dist & Vel Mats
distMap1 = sqrt(errTen1(:,:,1).^2+errTen1(:,:,2).^2+errTen1(:,:,3).^2);
velMap1 = sqrt(errTen1(:,:,4).^2+errTen1(:,:,5).^2+errTen1(:,:,6).^2);
distMap2 = sqrt(errTen2(:,:,1).^2+errTen2(:,:,2).^2+errTen2(:,:,3).^2);
velMap2 = sqrt(errTen2(:,:,4).^2+errTen2(:,:,5).^2+errTen2(:,:,6).^2);
errTen1 = cat(3,cat(3,errTen1,distMap1),velMap1);
errTen2 = cat(3,cat(3,errTen2,distMap2),velMap2);

% Create Mesh
[incMesh, eccMesh] = meshgrid(incRange,eccRange);
% Level Vectors for colorbar
levelsOe1 = [linspace(-max(abs(errTen1(:,:,1)-errTen2(:,:,1)),[],'all'),0,10),...
    linspace(0,max(abs(errTen1(:,:,1)-errTen2(:,:,1)),[],'all'),10)];

levelsOe2 = [linspace(-max(abs(errTen1(:,:,2)-errTen2(:,:,2)),[],'all'),0,10),...
    linspace(0,max(abs(errTen1(:,:,2)-errTen2(:,:,2)),[],'all'),10)];

levelsOe3 = [linspace(-max(abs(errTen1(:,:,3)-errTen2(:,:,3)),[],'all'),0,10),...
    linspace(0,max(abs(errTen1(:,:,3)-errTen2(:,:,3)),[],'all'),10)];

levelsOe4 = [linspace(-max(abs(errTen1(:,:,4)-errTen2(:,:,4)),[],'all'),0,10),...
    linspace(0,max(abs(errTen1(:,:,4)-errTen2(:,:,4)),[],'all'),10)];

levelsOe5 = [linspace(-max(abs(errTen1(:,:,5)-errTen2(:,:,5)),[],'all'),0,10),...
    linspace(0,max(abs(errTen1(:,:,5)-errTen2(:,:,5)),[],'all'),10)];

levelsOe6 = [linspace(-max(abs(errTen1(:,:,6)-errTen2(:,:,6)),[],'all'),0,10),...
    linspace(0,max(abs(errTen1(:,:,6)-errTen2(:,:,6)),[],'all'),10)];

levelsOe7 = [linspace(-max(abs(errTen1(:,:,7)-errTen2(:,:,7)),[],'all'),0,10),...
    linspace(0,max(abs(errTen1(:,:,7)-errTen2(:,:,7)),[],'all'),10)];

levelsOe8 = [linspace(-max(abs(errTen1(:,:,8)-errTen2(:,:,8)),[],'all'),0,10),...
    linspace(0,max(abs(errTen1(:,:,8)-errTen2(:,:,8)),[],'all'),10)];

% setup for loop
levelsOe = [levelsOe1;levelsOe2;levelsOe3;levelsOe4;levelsOe5;levelsOe6;levelsOe7;levelsOe8];

% Virtual colorfix mesh - makes a small mesh off screen with the full
% range of colors so we consistently get 0 as the same color
[xV,yV] = meshgrid([1,2],[1,2]);

for iOe = 1:8
    % Plot all
    figure(nFig*10+iOe)
    contourf(xV,yV,levelsOe(iOe,end)*[1,0;0,-1],levelsOe(iOe,:))
    c = colorbar;
    c.Label.Interpreter = 'latex';
    c.Label.String = ['$\Delta\tilde{' titleSet{iOe} '} \left[\mathrm{' unitSet{iOe} '}\right]$'];
    c.Label.FontSize = fontSizeAll;
    colormap jet
    shading interp
    hold on
    contourf(incMesh,eccMesh,errTen1(:,:,iOe)-errTen2(:,:,iOe),levelsOe(iOe,:),'LineColor','none')
    xlim([incRange(1),incRange(end)])
    ylim([eccRange(1),eccRange(end)])
    xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',fontSizeAll)
    ylabel(eccLabel,'interpreter','latex','fontsize',fontSizeAll)
    if logE
        yticks(eTicks)
        yticklabels(eTickLabels)
    end
    if isempty(saveFolder)
        title(['$' titleSet{iOe} ': ' label1 '-' label2 '$'],Interpreter='latex')
    else
        exportgraphics(gcf,[saveFolder '\RelMap' fileNameSet{iOe} '.eps'],"Resolution",600)
    end
    hold off
end

% Create Mesh
[incMesh, eccMesh] = meshgrid(incRange,eccRange);
% Level Vectors for colorbar - take into account both 1 & 2
levelsOe1 = linspace(0,max([max(errTen1(:,:,1),[],'all'),max(errTen2(:,:,1),[],'all')]),20);
levelsOe2 = linspace(0,max([max(errTen1(:,:,2),[],'all'),max(errTen2(:,:,2),[],'all')]),20);
levelsOe3 = linspace(0,max([max(errTen1(:,:,3),[],'all'),max(errTen2(:,:,3),[],'all')]),20);
levelsOe4 = linspace(0,max([max(errTen1(:,:,4),[],'all'),max(errTen2(:,:,4),[],'all')]),20);
levelsOe5 = linspace(0,max([max(errTen1(:,:,5),[],'all'),max(errTen2(:,:,5),[],'all')]),20);
levelsOe6 = linspace(0,max([max(errTen1(:,:,6),[],'all'),max(errTen2(:,:,6),[],'all')]),20);
levelsOe7 = linspace(0,max([max(errTen1(:,:,7),[],'all'),max(errTen2(:,:,7),[],'all')]),20);
levelsOe8 = linspace(0,max([max(errTen1(:,:,8),[],'all'),max(errTen2(:,:,8),[],'all')]),20);
% setup for loop
levelsOe = [levelsOe1;levelsOe2;levelsOe3;levelsOe4;levelsOe5;levelsOe6;levelsOe7;levelsOe8];

% Virtual colorfix mesh - makes a small mesh off screen with the full
% range of colors so we consistently get 0 as the same color
[xV,yV] = meshgrid([1,2],[1,2]);

% Abs Errors 1
for iOe = 1:8
    % Plot all
    figure(nFig*20+iOe)
    contourf(xV,yV,levelsOe(iOe,end)*[1,0;0,-1],levelsOe(iOe,:))
    c = colorbar;
    c.Label.Interpreter = 'latex';
    c.Label.String = ['$\tilde{' titleSet{iOe} '} \left[\mathrm{' unitSet{iOe} '}\right]$'];
    c.Label.FontSize = fontSizeAll;
    colormap jet
    shading interp
    hold on
    contourf(incMesh,eccMesh,errTen1(:,:,iOe),levelsOe(iOe,:),'LineColor','none')
    xlim([incRange(1),incRange(end)])
    ylim([eccRange(1),eccRange(end)])
    xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',fontSizeAll)
    ylabel(eccLabel,'interpreter','latex','fontsize',fontSizeAll)
    if logE
        yticks(eTicks)
        yticklabels(eTickLabels)
    end
    if isempty(saveFolder)
        title(['$' titleSet{iOe} ': ' label1 '$'],Interpreter='latex')
    else
        exportgraphics(gcf,[saveFolder '\ForAbsMap' fileNameSet{iOe} '.eps'],"Resolution",600)
    end
    hold off
end
% Abs Errors 2
for iOe = 1:8
    % Plot all
    figure(nFig*30+iOe)
    contourf(xV,yV,levelsOe(iOe,end)*[1,0;0,-1],levelsOe(iOe,:))
    c = colorbar;
    c.Label.Interpreter = 'latex';
    c.Label.String = ['$\tilde{' titleSet{iOe} '} \left[\mathrm{' unitSet{iOe} '}\right]$'];
    c.Label.FontSize = fontSizeAll;
    colormap jet
    shading interp
    hold on
    contourf(incMesh,eccMesh,errTen2(:,:,iOe),levelsOe(iOe,:),'LineColor','none')
    xlim([incRange(1),incRange(end)])
    ylim([eccRange(1),eccRange(end)])
    xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',fontSizeAll)
    ylabel(eccLabel,'interpreter','latex','fontsize',fontSizeAll)
    if logE
        yticks(eTicks)
        yticklabels(eTickLabels)
    end
    if isempty(saveFolder)
        title(['$' titleSet{iOe} ': ' label2 '$'],Interpreter='latex')
    else
        exportgraphics(gcf,[saveFolder '\KozAbsMap' fileNameSet{iOe} '.eps'],"Resolution",600)
    end
    hold off
end
end