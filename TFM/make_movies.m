function movies = make_movies(path)
load('RedBlueColormaps.mat')
load('BlackColormaps.mat')
load(fullfile(path, 'Allresults2.mat'));

%% define video parameters
length_video = size(Dx,3);
length_x = size(Dx,1);
length_y = size(Dx,2);
[PosX,PosY]=meshgrid(1:length_y,1:length_x);

h1=figure('units','Normalized','position',[0.45 0.05 0.4 0.4],'Name','Stress magnitude');
h2=figure('units','Normalized','position',[0.02 0.52 0.4 0.4],'Name','Displacements');

%% Displacement video
max_disp=max(Dx,[],'all');
for i=1:length_video
displacement = sqrt(Dx(:,:,i).^2+Dy(:,:,i).^2);
figure(h2) % Displacements
imagesc(displacement),colormap(blackcmap),colorbar, caxis([0 max_disp]), axis off, daspect([1 1 1]), %.*BWgrid
hold on;

if i == 1
    quiver(PosX(1:4:end,1:4:end),PosY(1:4:end,1:4:end),20e6*Dx(1:4:end,1:4:end,i),20e6*Dy(1:4:end,1:4:end,i),'r','AutoScale','off');
else
    quiver(PosX(1:4:end,1:4:end),PosY(1:4:end,1:4:end),20e6*Dx(1:4:end,1:4:end,i),20e6*Dy(1:4:end,1:4:end,i),'r','AutoScale','off');
end

if ~exist(fullfile(path,'movies'),'dir')
    mkdir(fullfile(path,'movies'))
end

figurepath=cat(2,path,'\movies\displacement',num2str(i),'.tif');
print(h2,fullfile(figurepath),'-dtiff','-r100');
A=imread(figurepath);
imwrite(A,cat(2,path,'\movies\displacement.tif'),'WriteMode','append');
delete (figurepath)
end

%% Traction force video
max_traction=max(Tx,[],'all');

for i=1:length_video
traction = sqrt(Tx(:,:,i).^2+Ty(:,:,i).^2);
figure(h1);
imagesc(traction),colormap(jet),colorbar, caxis([0 max_traction]), axis off, daspect([1 1 1]), %.*BWgrid
hold on;

if i == 1
    quiver(PosX(1:4:end,1:4:end),PosY(1:4:end,1:4:end),2e-2*Tx(1:4:end,1:4:end,i),2e-2*Ty(1:4:end,1:4:end,i),'r','AutoScale','off');
else
    quiver(PosX(1:4:end,1:4:end),PosY(1:4:end,1:4:end),2e-2*Tx(1:4:end,1:4:end,i),2e-2*Ty(1:4:end,1:4:end,i),'r','AutoScale','off');
end
if ~exist(fullfile(path,'movies'),'dir')
    mkdir(fullfile(path,'movies'))
end

figurepath=cat(2,path,'\movies\traction',num2str(i),'.tif');
print(h1,fullfile(figurepath),'-dtiff','-r100');
A=imread(figurepath);
imwrite(A,cat(2,path,'\movies\traction.tif'),'WriteMode','append');
delete (figurepath)
end

end