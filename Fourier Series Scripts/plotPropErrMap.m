function plotPropErrMap(eccRange,incRange,errTen1,errTen2,label1,label2,nFig)
%plotPropErrMap plots error maps for comaprison of propagation methods
%
% ~~~~~~~~~~~~~~  Inputs  ~~~~~~~~~~~~~~~
% errTen1, errTen2: nEcc x nInc x 6 Tensors
% label1, label2: label strings - optional
% nFig - 10s digit for figure numbering - default 0
%
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
    % Virtual colorfix mesh
    [xV,yV] = meshgrid([1,2],[1,2]);
    
    % Plot
    figure(nFig*10+1)
    contourf(xV,yV,levelsSma(end)*[1,0;0,-1],levelsSma)
    colorbar
    colormap jet
    shading interp
    hold on
    contourf(incMesh,eccMesh,errTen1(:,:,1)-errTen2(:,:,1),levelsSma,'LineColor','none')
    xlim([incRange(1),incRange(end)])
    ylim([eccRange(1),eccRange(end)])
    xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
    ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
    if ~isempty(label1)
        title(['sma: ' label1 '-' label2])
    end
    hold off
    
    figure(nFig*10+2)
    contourf(xV,yV,levelsEcc(end)*[1,0;0,-1],levelsEcc)
    colorbar
    colormap jet
    shading interp
    hold on
    contourf(incMesh,eccMesh,errTen1(:,:,2)-errTen2(:,:,2),levelsEcc,'LineColor','none')
    xlim([incRange(1),incRange(end)])
    ylim([eccRange(1),eccRange(end)])
    xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
    ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
    if ~isempty(label1)
        title(['ecc: ' label1 '-' label2])
    end
    hold off
    
    figure(nFig*10+3)
    contourf(xV,yV,levelsInc(end)*[1,0;0,-1],levelsInc)
    colorbar
    colormap jet
    shading interp
    hold on
    contourf(incMesh,eccMesh,errTen1(:,:,3)-errTen2(:,:,3),levelsInc,'LineColor','none')
    xlim([incRange(1),incRange(end)])
    ylim([eccRange(1),eccRange(end)])
    xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
    ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
    if ~isempty(label1)
        title(['inc: ' label1 '-' label2])
    end
    hold off
    
    figure(nFig*10+4)
    contourf(xV,yV,levelsRan(end)*[1,0;0,-1],levelsRan)
    colorbar
    colormap jet
    shading interp
    hold on
    contourf(incMesh,eccMesh,errTen1(:,:,4)-errTen2(:,:,4),levelsRan,'LineColor','none')
    xlim([incRange(1),incRange(end)])
    ylim([eccRange(1),eccRange(end)])
    xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
    ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
    if ~isempty(label1)
        title(['raan: ' label1 '-' label2])
    end
    hold off
    
    figure(nFig*10+5)
    contourf(xV,yV,levelsAop(end)*[1,0;0,-1],levelsAop)
    colorbar
    colormap jet
    shading interp
    hold on
    contourf(incMesh,eccMesh,errTen1(:,:,5)-errTen2(:,:,5),levelsAop,'LineColor','none')
    xlim([incRange(1),incRange(end)])
    ylim([eccRange(1),eccRange(end)])
    xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
    ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
    if ~isempty(label1)
        title(['aop: ' label1 '-' label2])
    end
    hold off
    
    figure(nFig*10+6)
    contourf(xV,yV,levelsMan(end)*[1,0;0,-1],levelsMan)
    colorbar
    colormap jet
    shading interp
    hold on
    contourf(incMesh,eccMesh,errTen1(:,:,6)-errTen2(:,:,6),levelsMan,'LineColor','none')
    xlim([incRange(1),incRange(end)])
    ylim([eccRange(1),eccRange(end)])
    xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
    ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
    if ~isempty(label1)
        title(['M: ' label1 '-' label2])
    end
    hold off
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
    % Virtual colorfix mesh
    [xV,yV] = meshgrid([1,2],[1,2]);
    
    
    % Plot
    figure(nFig*10+1)
    contourf(xV,yV,levelsSma(end)*[1,0;0,1],levelsSma)
    colorbar
    colormap jet
    shading interp
    hold on
    contourf(incMesh,eccMesh,errTen1(:,:,1),levelsSma,'LineColor','none')
    xlim([incRange(1),incRange(end)])
    ylim([eccRange(1),eccRange(end)])
    xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
    ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
    if ~isempty(label1)
        title(['sma: ' label1])
    end
    hold off
    
    figure(nFig*10+2)
    contourf(xV,yV,levelsEcc(end)*[1,0;0,1],levelsEcc)
    colorbar
    colormap jet
    shading interp
    hold on
    contourf(incMesh,eccMesh,errTen1(:,:,2),levelsEcc,'LineColor','none')
    xlim([incRange(1),incRange(end)])
    ylim([eccRange(1),eccRange(end)])
    xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
    ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
    if ~isempty(label1)
        title(['ecc: ' label1])
    end
    hold off
    
    figure(nFig*10+3)
    contourf(xV,yV,levelsInc(end)*[1,0;0,1],levelsInc)
    colorbar
    colormap jet
    shading interp
    hold on
    contourf(incMesh,eccMesh,errTen1(:,:,3),levelsInc,'LineColor','none')
    xlim([incRange(1),incRange(end)])
    ylim([eccRange(1),eccRange(end)])
    xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
    ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
    if ~isempty(label1)
        title(['inc: ' label1])
    end
    hold off
    
    figure(nFig*10+4)
    contourf(xV,yV,levelsRan(end)*[1,0;0,1],levelsRan)
    colorbar
    colormap jet
    shading interp
    hold on
    contourf(incMesh,eccMesh,errTen1(:,:,4),levelsRan,'LineColor','none')
    xlim([incRange(1),incRange(end)])
    ylim([eccRange(1),eccRange(end)])
    xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
    ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
    if ~isempty(label1)
        title(['raan: ' label1])
    end
    hold off
    
    figure(nFig*10+5)
    contourf(xV,yV,levelsAop(end)*[1,0;0,1],levelsAop)
    colorbar
    colormap jet
    shading interp
    hold on
    contourf(incMesh,eccMesh,errTen1(:,:,5),levelsAop,'LineColor','none')
    xlim([incRange(1),incRange(end)])
    ylim([eccRange(1),eccRange(end)])
    xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
    ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
    if ~isempty(label1)
        title(['aop: ' label1])
    end
    hold off
    
    figure(nFig*10+6)
    contourf(xV,yV,levelsMan(end)*[1,0;0,1],levelsMan)
    colorbar
    colormap jet
    shading interp
    hold on
    contourf(incMesh,eccMesh,errTen1(:,:,6),levelsMan,'LineColor','none')
    xlim([incRange(1),incRange(end)])
    ylim([eccRange(1),eccRange(end)])
    xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
    ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
    if ~isempty(label1)
        title(['M: ' label1])
    end
    hold off
    
end