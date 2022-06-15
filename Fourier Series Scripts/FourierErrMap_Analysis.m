%% Load Data


%% Plot Relative Errors

% Creat Mesh
[incMesh, eccMesh] = meshgrid(incRange,eccRange);
% Level Vectors for colorbar
levelsSma = [linspace(-max(abs(errTenF(:,:,1)-errTenB(:,:,1)),[],'all'),0,10),...
    linspace(0,max(abs(errTenF(:,:,1)-errTenB(:,:,1)),[],'all'),10)];

levelsEcc = [linspace(-max(abs(errTenF(:,:,2)-errTenB(:,:,2)),[],'all'),0,10),...
    linspace(0,max(abs(errTenF(:,:,2)-errTenB(:,:,2)),[],'all'),10)];

levelsInc = [linspace(-max(abs(errTenF(:,:,3)-errTenB(:,:,3)),[],'all'),0,10),...
    linspace(0,max(abs(errTenF(:,:,3)-errTenB(:,:,3)),[],'all'),10)];

levelsRan = [linspace(-max(abs(errTenF(:,:,4)-errTenB(:,:,4)),[],'all'),0,10),...
    linspace(0,max(abs(errTenF(:,:,4)-errTenB(:,:,4)),[],'all'),10)];

levelsAop = [linspace(-max(abs(errTenF(:,:,5)-errTenB(:,:,5)),[],'all'),0,10),...
    linspace(0,max(abs(errTenF(:,:,5)-errTenB(:,:,5)),[],'all'),10)];

levelsMan = [linspace(-max(abs(errTenF(:,:,6)-errTenB(:,:,6)),[],'all'),0,10),...
    linspace(0,max(abs(errTenF(:,:,6)-errTenB(:,:,6)),[],'all'),10)];
% Virtual colorfix mesh
[xV,yV] = meshgrid([1,2],[1,2]);


% Plot
figure(1)
contourf(xV,yV,levelsSma(end)*[1,0;0,-1],levelsSma)
colorbar
colormap jet
shading interp
hold on
contourf(incMesh,eccMesh,errTenF(:,:,1)-errTenB(:,:,1),levelsSma,'LineColor','none')
xlim([incRange(1),incRange(end)])
ylim([eccRange(1),eccRange(end)])
xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
hold off

figure(2)
contourf(xV,yV,levelsEcc(end)*[1,0;0,-1],levelsEcc)
colorbar
colormap jet
shading interp
hold on
contourf(incMesh,eccMesh,errTenF(:,:,2)-errTenB(:,:,2),levelsEcc,'LineColor','none')
xlim([incRange(1),incRange(end)])
ylim([eccRange(1),eccRange(end)])
xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
hold off

figure(3)
contourf(xV,yV,levelsInc(end)*[1,0;0,-1],levelsInc)
colorbar
colormap jet
shading interp
hold on
contourf(incMesh,eccMesh,errTenF(:,:,3)-errTenB(:,:,3),levelsInc,'LineColor','none')
xlim([incRange(1),incRange(end)])
ylim([eccRange(1),eccRange(end)])
xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
hold off

figure(4)
contourf(xV,yV,levelsRan(end)*[1,0;0,-1],levelsRan)
colorbar
colormap jet
shading interp
hold on
contourf(incMesh,eccMesh,errTenF(:,:,4)-errTenB(:,:,4),levelsRan,'LineColor','none')
xlim([incRange(1),incRange(end)])
ylim([eccRange(1),eccRange(end)])
xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
hold off

figure(5)
contourf(xV,yV,levelsAop(end)*[1,0;0,-1],levelsAop)
colorbar
colormap jet
shading interp
hold on
contourf(incMesh,eccMesh,errTenF(:,:,5)-errTenB(:,:,5),levelsAop,'LineColor','none')
xlim([incRange(1),incRange(end)])
ylim([eccRange(1),eccRange(end)])
xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
hold off

figure(6)
contourf(xV,yV,levelsMan(end)*[1,0;0,-1],levelsMan)
colorbar
colormap jet
shading interp
hold on
contourf(incMesh,eccMesh,errTenF(:,:,6)-errTenB(:,:,6),levelsMan,'LineColor','none')
xlim([incRange(1),incRange(end)])
ylim([eccRange(1),eccRange(end)])
xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
hold off

%% Plot Absolute Errors - Brouwer
% Creat Mesh
[incMesh, eccMesh] = meshgrid(incRange,eccRange);
% Level Vectors for colorbar
levelsSma = linspace(0,max(errTenB(:,:,1),[],'all'),20);
levelsEcc = linspace(0,max(errTenB(:,:,2),[],'all'),20);
levelsInc = linspace(0,max(errTenB(:,:,3),[],'all'),20);
levelsRan = linspace(0,max(errTenB(:,:,4),[],'all'),20);
levelsAop = linspace(0,max(errTenB(:,:,5),[],'all'),20);
levelsMan = linspace(0,max(errTenB(:,:,6),[],'all'),20);
% Virtual colorfix mesh
[xV,yV] = meshgrid([1,2],[1,2]);


% Plot
figure(11)
contourf(xV,yV,levelsSma(end)*[1,0;0,1],levelsSma)
colorbar
colormap jet
shading interp
hold on
contourf(incMesh,eccMesh,errTenB(:,:,1),levelsSma,'LineColor','none')
xlim([incRange(1),incRange(end)])
ylim([eccRange(1),eccRange(end)])
xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
hold off

figure(12)
contourf(xV,yV,levelsEcc(end)*[1,0;0,1],levelsEcc)
colorbar
colormap jet
shading interp
hold on
contourf(incMesh,eccMesh,errTenB(:,:,2),levelsEcc,'LineColor','none')
xlim([incRange(1),incRange(end)])
ylim([eccRange(1),eccRange(end)])
xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
hold off

figure(13)
contourf(xV,yV,levelsInc(end)*[1,0;0,1],levelsInc)
colorbar
colormap jet
shading interp
hold on
contourf(incMesh,eccMesh,errTenB(:,:,3),levelsInc,'LineColor','none')
xlim([incRange(1),incRange(end)])
ylim([eccRange(1),eccRange(end)])
xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
hold off

figure(14)
contourf(xV,yV,levelsRan(end)*[1,0;0,1],levelsRan)
colorbar
colormap jet
shading interp
hold on
contourf(incMesh,eccMesh,errTenB(:,:,4),levelsRan,'LineColor','none')
xlim([incRange(1),incRange(end)])
ylim([eccRange(1),eccRange(end)])
xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
hold off

figure(15)
contourf(xV,yV,levelsAop(end)*[1,0;0,1],levelsAop)
colorbar
colormap jet
shading interp
hold on
contourf(incMesh,eccMesh,errTenB(:,:,5),levelsAop,'LineColor','none')
xlim([incRange(1),incRange(end)])
ylim([eccRange(1),eccRange(end)])
xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
hold off

figure(16)
contourf(xV,yV,levelsMan(end)*[1,0;0,1],levelsMan)
colorbar
colormap jet
shading interp
hold on
contourf(incMesh,eccMesh,errTenB(:,:,6),levelsMan,'LineColor','none')
xlim([incRange(1),incRange(end)])
ylim([eccRange(1),eccRange(end)])
xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
hold off

%% Plot Absolute Errors - Fourier
% Creat Mesh
[incMesh, eccMesh] = meshgrid(incRange,eccRange);
% Level Vectors for colorbar
levelsSma = linspace(0,max(errTenF(:,:,1),[],'all'),20);
levelsEcc = linspace(0,max(errTenF(:,:,2),[],'all'),20);
levelsInc = linspace(0,max(errTenF(:,:,3),[],'all'),20);
levelsRan = linspace(0,max(errTenF(:,:,4),[],'all'),20);
levelsAop = linspace(0,max(errTenF(:,:,5),[],'all'),20);
levelsMan = linspace(0,max(errTenF(:,:,6),[],'all'),20);
% Virtual colorfix mesh
[xV,yV] = meshgrid([1,2],[1,2]);


% Plot
figure(21)
contourf(xV,yV,levelsSma(end)*[1,0;0,1],levelsSma)
colorbar
colormap jet
shading interp
hold on
contourf(incMesh,eccMesh,errTenF(:,:,1),levelsSma,'LineColor','none')
xlim([incRange(1),incRange(end)])
ylim([eccRange(1),eccRange(end)])
xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
hold off

figure(22)
contourf(xV,yV,levelsEcc(end)*[1,0;0,1],levelsEcc)
colorbar
colormap jet
shading interp
hold on
contourf(incMesh,eccMesh,errTenF(:,:,2),levelsEcc,'LineColor','none')
xlim([incRange(1),incRange(end)])
ylim([eccRange(1),eccRange(end)])
xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
hold off

figure(23)
contourf(xV,yV,levelsInc(end)*[1,0;0,1],levelsInc)
colorbar
colormap jet
shading interp
hold on
contourf(incMesh,eccMesh,errTenF(:,:,3),levelsInc,'LineColor','none')
xlim([incRange(1),incRange(end)])
ylim([eccRange(1),eccRange(end)])
xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
hold off

figure(24)
contourf(xV,yV,levelsRan(end)*[1,0;0,1],levelsRan)
colorbar
colormap jet
shading interp
hold on
contourf(incMesh,eccMesh,errTenF(:,:,4),levelsRan,'LineColor','none')
xlim([incRange(1),incRange(end)])
ylim([eccRange(1),eccRange(end)])
xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
hold off

figure(25)
contourf(xV,yV,levelsAop(end)*[1,0;0,1],levelsAop)
colorbar
colormap jet
shading interp
hold on
contourf(incMesh,eccMesh,errTenF(:,:,5),levelsAop,'LineColor','none')
xlim([incRange(1),incRange(end)])
ylim([eccRange(1),eccRange(end)])
xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
hold off

figure(26)
contourf(xV,yV,levelsMan(end)*[1,0;0,1],levelsMan)
colorbar
colormap jet
shading interp
hold on
contourf(incMesh,eccMesh,errTenF(:,:,6),levelsMan,'LineColor','none')
xlim([incRange(1),incRange(end)])
ylim([eccRange(1),eccRange(end)])
xlabel('$\rm{Inclination} \left[deg\right]$','interpreter','latex','fontsize',12)
ylabel('$\rm{Eccentricity}$','interpreter','latex','fontsize',12)
hold off