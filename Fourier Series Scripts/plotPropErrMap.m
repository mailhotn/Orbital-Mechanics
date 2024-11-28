function plotPropErrMap(eccRange,incRange,errTen1,errTen2,label1,label2,nFig,logE)
%plotPropErrMap plots error maps for comaprison of propagation methods
%
% ~~~~~~~~~~~~~~  Inputs  ~~~~~~~~~~~~~~~
% errTen1, errTen2: nEcc x nInc x 6 Tensors
% label1, label2: label strings - optional
% nFig - 10s digit for figure numbering - default 0
% logE - flag to use logarithmic e Axis

if nargin < 8
    logE = false;
end
if nargin < 7
    nFig = 0;
end
if nargin < 5
    label1 = [];
    label2 = [];
end
if nargin == 3
    errTen2 = [];
end
if logE
    eccRange = log10(eccRange);
    eccLabel = '$\rm{Eccentricity}$'; % switched to labeling ecc directly
    eTicks = flip(eccRange(end):-0.2:eccRange(1)).';
    eTicks(1) = eccRange(1);
    eTickLabels = num2str(10.^eTicks,2);
elseif ~logE
    eccLabel = '$\rm{Eccentricity}$';
end

if ~isempty(errTen2) % errTen2 exists - Comparison
    % Create Mesh
    [incMesh, eccMesh] = meshgrid(incRange,eccRange);
    % Level Vectors for colorbar
    levelsSma = [linspace(-max(abs(errTen1(:,:,1)-errTen2(:,:,1)),[],'all'),0,10),...
        linspace(0,max(abs(errTen1(:,:,1)-errTen2(:,:,1)),[],'all'),10)];
    
    levelsEcc = [linspace(-max(abs(errTen1(:,:,2)-errTen2(:,:,2)),[],'all'),0,10),...
        linspace(0,max(abs(errTen1(:,:,2)-errTen2(:,:,2)),[],'all'),10)];
    
    levelsInc = [linspace(-max(abs(errTen1(:,:,3)-errTen2(:,:,3)),[],'all'),0,10),...
        linspace(0,max(abs(errTen1(:,:,3)-errTen2(:,:,3)),[],'all'),10)];
    
    levelsRan = [linspace(-max(abs(errTen1(:,:,4)-errTen2(:,:,4)),[],'all'),0,10),...
        linspace(0,max(abs(errTen1(:,:,4)-errTen2(:,:,4)),[],'all'),10)];
    
    levelsAop = [linspace(-max(abs(errTen1(:,:,5)-errTen2(:,:,5)),[],'all'),0,10),...
        linspace(0,max(abs(errTen1(:,:,5)-errTen2(:,:,5)),[],'all'),10)];
    
    levelsMan = [linspace(-max(abs(errTen1(:,:,6)-errTen2(:,:,6)),[],'all'),0,10),...
        linspace(0,max(abs(errTen1(:,:,6)-errTen2(:,:,6)),[],'all'),10)];
    
    % setup for loop
    levelsOe = [levelsSma;levelsEcc;levelsInc;levelsRan;levelsAop;levelsMan];
    titleSet = {'sma','ecc','inc','raan','aop','M'};
    
    % Virtual colorfix mesh - makes a small mesh off screen with the full
    % range of colors so we consistently get 0 as the same color
    [xV,yV] = meshgrid([1,2],[1,2]);
    
    for iOe = 1:6
     % Plot all
    figure(nFig*10+iOe)
    contourf(xV,yV,levelsOe(iOe,end)*[1,0;0,-1],levelsOe(iOe,:))
    colorbar
    colormap jet
    shading interp
    hold on
    contourf(incMesh,eccMesh,errTen1(:,:,iOe)-errTen2(:,:,iOe),levelsOe(iOe,:),'LineColor','none')
    xlim([incRange(1),incRange(end)])
    ylim([eccRange(1),eccRange(end)])
    xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
    ylabel(eccLabel,'interpreter','latex','fontsize',12)
    if logE
        yticks(eTicks)
        yticklabels(eTickLabels)
    end
    if ~isempty(label1)
        title([titleSet{iOe} ': ' label1 '-' label2])
    end
    hold off
    end
else % no ErrTen2 - Absolute Error
    
    % Create Mesh
    [incMesh, eccMesh] = meshgrid(incRange,eccRange);
    % Level Vectors for colorbar
    levelsSma = linspace(0,max(errTen1(:,:,1),[],'all'),20);
    levelsEcc = linspace(0,max(errTen1(:,:,2),[],'all'),20);
    levelsInc = linspace(0,max(errTen1(:,:,3),[],'all'),20);
    levelsRan = linspace(0,max(errTen1(:,:,4),[],'all'),20);
    levelsAop = linspace(0,max(errTen1(:,:,5),[],'all'),20);
    levelsMan = linspace(0,max(errTen1(:,:,6),[],'all'),20);
    % setup for loop
    levelsOe = [levelsSma;levelsEcc;levelsInc;levelsRan;levelsAop;levelsMan];
    titleSet = {'sma','ecc','inc','raan','aop','M'};
    
    % Virtual colorfix mesh - makes a small mesh off screen with the full
    % range of colors so we consistently get 0 as the same color
    [xV,yV] = meshgrid([1,2],[1,2]);
    
    for iOe = 1:6
     % Plot all
    figure(nFig*10+iOe)
    contourf(xV,yV,levelsOe(iOe,end)*[1,0;0,-1],levelsOe(iOe,:))
    colorbar
    colormap jet
    shading interp
    hold on
    contourf(incMesh,eccMesh,errTen1(:,:,iOe),levelsOe(iOe,:),'LineColor','none')
    xlim([incRange(1),incRange(end)])
    ylim([eccRange(1),eccRange(end)])
    xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
    ylabel(eccLabel,'interpreter','latex','fontsize',12)
    if logE
        yticks(eTicks)
        yticklabels(eTickLabels)
    end
    if ~isempty(label1)
        title([titleSet{iOe} ': ' label1])
    end
    hold off
    end  
end