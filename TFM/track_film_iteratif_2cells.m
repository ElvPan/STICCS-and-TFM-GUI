%Changes on 11/04/2016: - introduce overlap between tracked windows
%change on 26/10/2016: iterative PIV before tracking
%change on 12/03/2018: adapted for 2 cells
% change on 16/04/2018: embedded in automated analysis, changed to be a function

function track = track_film_iteratif_2cells(beads_file,beads_path,initial_file,initial_path, brightfield_file, brightfield_path, mask_left_path, mask_left_file, mask_right_path, mask_right_file)
%% Load films

load('BlackColormaps.mat')


%[tempfilename,temppath]=uigetfile(fullfile(path,'*.*'),'Sequence d''images billes a traiter');
if beads_file
    path = beads_path;
    filename=beads_file;
    info = imfinfo(fullfile(path,filename),'tif');
    Nb_frames=numel(info);
    I=imread(fullfile(path,filename),'tif');
    
    cl=class(I);
    stressed=zeros([size(I) Nb_frames], cl);
    stressed(:,:,1)=I;
    for frame=2:Nb_frames
        stressed(:,:,frame)=imread(fullfile(path,filename),frame);
    end
end

%[tempfilename,temppath]=uigetfile(fullfile(path,'*.*'),'Image bille nulle');
if initial_file
    path=initial_path;
    nonstressed=imread(fullfile(initial_path,initial_file));
end

%[tempfilename,temppath]=uigetfile(fullfile(path,'*.*'),'Sequence d''images de cellules (facult)');
if brightfield_file
    I=imread(fullfile(brightfield_path,brightfield_file),'tif');
    cellule=zeros([size(I) Nb_frames], class(I));
    cellule(:,:,1)=contrast(I,0.1,0.1);
    for frame=2:Nb_frames
        cellule(:,:,frame)=contrast(imread(fullfile(brightfield_path,brightfield_file),frame),0.1,0.1);
    end
else
    cellule=stressed;
end

%[tempfilename,temppath]=uigetfile(fullfile(path,'*.*'),'Sequence de masques cellule 1');
if mask_left_file
    I=imread(fullfile(mask_left_path,mask_left_file),'tif');
    mask1=zeros([size(I) Nb_frames], class(I));
    mask1(:,:,1)=I;
    no_images=numel(imfinfo(cat(2,mask_left_path,mask_left_file)));
    for frame=2:Nb_frames
        if no_images==1
            mask1(:,:,frame)=imread(fullfile(mask_left_path,mask_left_file));
        else
            mask1(:,:,frame)=imread(fullfile(mask_left_path,mask_left_file),frame);
        end
    end
end

% [tempfilename,temppath]=uigetfile(fullfile(path,'*.*'),'Sequence de masques cellule 2');
if mask_right_file
    I=imread(fullfile(mask_right_path,mask_right_file),'tif');
    mask2=zeros([size(I) Nb_frames], class(I));
    mask2(:,:,1)=I;
    no_images=numel(imfinfo(cat(2,mask_right_path,mask_right_file)));
    for frame=2:Nb_frames
        if no_images==1
            mask2(:,:,frame)=imread(fullfile(mask_right_path,mask_right_file));
        else
            mask2(:,:,frame)=imread(fullfile(mask_right_path,mask_right_file),frame);
        end
    end
end
%% Begin analysis
%%%%%%%% Parameters to play with %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pix = 4.26e-8;      % pixel size in meter
E=1670;             % Young modulus in Pa
nu=0.5;             % Poisson ratio
alphadef=1e-19;     % coef of regularization should vary with the young modulus

%initial window size
window=128 ;        %

%number of iterations
iter=0;             % after each interation the window size is divided by 2
overlapPIV=64;      % in pixels (applied to first window)               %32  %64
overlapTrack=32;    % in pixels (applied to last window befor tracking)  %16  %32
interval=8;
pas=interval*pix;

%Tracking parameters 
featsize=6;         % minimal particle size in pixel
barrg=15;           % maximal particle size in pixel ?
barcc=0.3;          % eccentricity (0 = round 1 = not round at all)
IdivRg=0;         % minimal ratio intensity/radius
masscut=350000;     % minimal intensity: the higher, the fewer particles are detected
maxd=5;             % maximal displacement of one particle between the 2 images %5
r_neighbor = 0;     % minimum distance to the nearest neighbor in pixel
corr_th = 0.5;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h1=figure('units','Normalized','position',[0.02 0.05 0.9 0.83],'Name','Tracking');
h2=figure('units','Normalized','position',[0.02 0.05 0.4 0.4],'Name','Stress vectors');
h3=figure('units','Normalized','position',[0.45 0.05 0.4 0.4],'Name','Stress magnitude');
h4=figure('units','Normalized','position',[0.02 0.52 0.4 0.4],'Name','Displacements');

[n,s,p_corr,regx,regy]=registration_ssinterp(nonstressed,stressed(:,:,1),cellule(:,:,1));
[dimy,dimx]=size(n);
period=window-overlapPIV;
nbx=floor((dimx-window)/period)+1;
nby=floor((dimy-window)/period)+1;
[posx,posy]=meshgrid((0:nbx-1)*period+floor((dimx-(nbx-1)*period)/2)+1,(0:nby-1)*period+floor((dimy-(nby-1)*period)/2)+1);
nodeX=floor(period*nbx/interval);
nodeY=floor(period*nby/interval);
[gridX,gridY]=meshgrid((0:nodeX-1)*interval+floor((dimx-nbx*period)/2)+1+interval/2,(0:nodeY-1)*interval+floor((dimy-nby*period)/2)+1+interval/2);

Mask1=false([size(n) Nb_frames]);
Mask2=false([size(n) Nb_frames]);
Dx=zeros([size(gridX) Nb_frames]);
Dy=zeros([size(gridX) Nb_frames]);
Tx=zeros([size(gridX) Nb_frames]);
Ty=zeros([size(gridX) Nb_frames]);
Brightfield=uint8(zeros([size(n) Nb_frames]));

%Calculated global parameters
% Pour le doublet de cellules
MOY2cells=zeros(1,Nb_frames);% contrainte moyenne
Total2cells=zeros(1,Nb_frames);% Somme en norme des tractions cellules-substrat
Vect2cells=zeros(1,Nb_frames);% somme vectorielle des tractions cellules-subtrat: a l'�quilibre, devrait etre egale a zero (systeme isole)
Ec2cells=zeros(1,Nb_frames);%?nergie contractile
Fcellcell=zeros(1,Nb_frames);% force cellule-cellule
DFcellcell=zeros(1,Nb_frames);% estimation de l'incertitude sur force cellule-cellule
Moment=zeros(1,Nb_frames);%1st order moment of the traction distribution
Angle=zeros(1,Nb_frames);% angle d'orientation principale par rapport � la verticale
Polar=zeros(1,Nb_frames);%polarization degree, calculated from the eigenvalues of traction 1st order moment (l1-l2)/(l1+l2)

%Pour chaque cellule du doublet
MOYcell1=zeros(1,Nb_frames);% contrainte moyenne
MOYcell2=zeros(1,Nb_frames);
Tmaxcell1=zeros(1,Nb_frames);%maximum stress
Tmaxcell2=zeros(1,Nb_frames);%maximum stress
Totalcell1=zeros(1,Nb_frames);% Somme en norme des tractions cellule1-substrat
Totalcell2=zeros(1,Nb_frames);% Somme en norme des tractions cellule2-substrat
Cell_ECMratio1=zeros(1,Nb_frames);%rapport Fcell-cell/Fcell-ECM pour cellule 1 =Totalcell1/Fcellcell
Cell_ECMratio2=zeros(1,Nb_frames);%rapport Fcell-cell/Fcell-ECM pour cellule 2
Momentcell1=zeros(1,Nb_frames);% moment contractile net cell1
Momentcell2=zeros(1,Nb_frames);% moment contractile net cell2
%Parametres rajoutes (version Artur)

Forcex_cell1=zeros(1,Nb_frames);%sum of x component of forces (in magnitude) for cell 1
Forcey_cell1=zeros(1,Nb_frames);%sum of y component of forces (in magnitude) for cell 1
Forcex_cell2=zeros(1,Nb_frames);%sum of x component of forces (in magnitude) for cell 2
Forcey_cell2=zeros(1,Nb_frames);%sum of y component of forces (in magnitude) for cell 2
Deplx_cell1=zeros(1,Nb_frames);%sum of x component of displacements (in magnitude) for cell 1
Deply_cell1=zeros(1,Nb_frames);%sum of x component of displacements (in magnitude) for cell 1
Deplx_cell2=zeros(1,Nb_frames);%sum of x component of displacements (in magnitude) for cell 1
Deply_cell2=zeros(1,Nb_frames);%sum of x component of displacements (in magnitude) for cell 1



for ff=1:Nb_frames 
         
    disp(['Analyse: image ',num2str(ff)]);
    if stressed(:,:,ff)==0
        disp('Null image!!!')
       return
    end
    
    [n,s,p_corr,m1_corr,regx,regy]=registration_ssinterp_pattern(nonstressed,stressed(:,:,ff),cellule(:,:,ff),mask1(:,:,ff));
    [m2_corr]=corr_mask(regx,regy,mask2(:,:,ff));
    restex=regx-round(regx);
    restey=regy-round(regy);
    if size(n,1)<dimy
        temp=zeros(dimy-size(n,1),size(n,2));
         n=[n;temp];
         s=[s;temp];
         mask1_corr=logical([m1_corr;temp]);
         mask2_corr=logical([m2_corr;temp]);
         bf=[p_corr;temp];
    elseif size(n,1)>dimy
        mask1_corr=logical(m1_corr(1:dimy,:));
        mask2_corr=logical(m2_corr(1:dimy,:));
        bf=p_corr(1:dimy,:);
    else 
        mask1_corr=logical(m1_corr);
        mask2_corr=logical(m2_corr);
        bf=p_corr;
    end
    if size(n,2)<dimx
        temp=zeros(size(n,1),dimx-size(n,2));
         n=[n temp];
         s=[s temp];
         tempmask=false(size(mask1_corr,1),dimx-size(mask1_corr,2));
         mask1_corr=[mask1_corr tempmask];
         mask2_corr=[mask2_corr tempmask];
         bf=[bf uint16(tempmask)];
    elseif size(n,2)>dimx
        mask1_corr=mask1_corr(:,1:dimx);
        mask2_corr=mask2_corr(:,1:dimx);
        bf=bf(:,1:dimx);
    end
mask_corr=mask1_corr|mask2_corr;
deplaceX=zeros(10000,1);
deplaceY=zeros(10000,1);
positionX=zeros(10000,1);
positionY=zeros(10000,1);
compteur=0;
mauvaises_fenetres=0;
[M,N]=size(s);

disp(['Nombre de fenetres: ',num2str(nby*nbx)]);
for i=1:nbx  
for j=1:nby
         mauvais_track=false;
       
       %nombre de fen�tres
       Nwin=2.^(2*(0:iter));  
       %tailles des fen�tres
       mi=window./(2.^(0:iter));
       %tableau de positions (centre de fen�tres): � chaque �tape le nb de positions est
       %multipli� par 4
       posix=cell(1,iter+1);
       posiy=cell(1,iter+1);
       posix{1}=posx(1,i);
       posiy{1}=posy(j,1);
       % tableau d'offsets
       offsetx=cell(1,iter+1);
       offsety=cell(1,iter+1);
       offsetx{1}=0;
       offsety{1}=0;
     
       
       for k=0:iter
           %d�finit les offsets et les positions pour chaque it�ration
           offsetx{k+1}=zeros(2^k);
           offsety{k+1}=zeros(2^k);
           %recopie les valeurs d'avant
           for l=1:2^(k-1)
               for m=1:2^(k-1)
               offsetx{k+1}((l-1)*2+1:l*2,(m-1)*2+1:m*2)=offsetx{k}(l,m)*ones(2);
               offsety{k+1}((l-1)*2+1:l*2,(m-1)*2+1:m*2)=offsety{k}(l,m)*ones(2);

               end
           end
           [px,py]=meshgrid((-mi(k+1)/2*(2^k-1):mi(k+1):mi(k+1)/2*(2^k-1)),(-mi(k+1)/2*(2^k-1):mi(k+1):mi(k+1)/2*(2^k-1)));
           posix{k+1}=posx(1,i)+px;
           posiy{k+1}=posy(j,1)+py;

           for nw=1:Nwin(k+1) 
              
               
               [sub1,sub2]=ind2sub([2^k 2^k],nw);
               %1. interpolation de n
               if k>0
            xi=(1:mi(k+1))-offsetx{k+1}(sub1,sub2)+posix{k+1}(sub1,sub2)-mi(k+1)/2-1;
            yi=(1:mi(k+1))-offsety{k+1}(sub1,sub2)+posiy{k+1}(sub1,sub2)-mi(k+1)/2-1;
            n_win = interp2(double(n),xi',yi);
               else
            n_win=n(posiy{k+1}(sub1,sub2)-mi(k+1)/2:posiy{k+1}(sub1,sub2)+mi(k+1)/2-1,posix{k+1}(sub1,sub2)-mi(k+1)/2:posix{k+1}(sub1,sub2)+mi(k+1)/2-1);
               end
           %2. mesure du d�calage
          
         s_win=s(posiy{k+1}(sub1,sub2)-mi(k+1)/2:posiy{k+1}(sub1,sub2)+mi(k+1)/2-1,posix{k+1}(sub1,sub2)-mi(k+1)/2:posix{k+1}(sub1,sub2)+mi(k+1)/2-1); 
%            n_win=n(posy(j,1)-window/2:posy(j,1)+window/2-1,posx(1,i)-window/2:posx(1,i)+window/2-1);
%            s_win=s(posy(j,1)-window/2:posy(j,1)+window/2-1,posx(1,i)-window/2:posx(1,i)+window/2-1); 
        [offx,offy,Cmax]=decalage(s_win,n_win);
         deltaI=min((max(n_win(:))-min(n_win(:)))/mean(n_win(:)),(max(s_win(:))-min(s_win(:)))/mean(s_win(:)));
         if deltaI<5||Cmax<0.5%||sqrt(offx^2+offy^2)>mi(k+1)/3
             offx=0;
             offy=0;
         end

         %remplissage du tableau de valeurs
         if k==0
         offsetx{1}=offx;
         offsety{1}=offy; 
         else
            offsetx{k+1}(sub1,sub2)=offx+ offsetx{k}(ceil(sub1/2),ceil(sub2/2));
            offsety{k+1}(sub1,sub2)=offy+ offsety{k}(ceil(sub1/2),ceil(sub2/2));
%          [sub1,sub2]=ind2sub([2^k 2^k],nw);
%          offsetx(1+(sub1-1)*2^(iter-k):(sub1)*2^(iter-k),1+(sub2-1)*2^(iter-k):(sub2)*2^(iter-k))=offsetx(1+(sub1-1)*2^(iter-k):(sub1)*2^(iter-k),1+(sub2-1)*2^(iter-k):(sub2)*2^(iter-k))+offx;
%          offsety(1+(sub1-1)*2^(iter-k):(sub1)*2^(iter-k),1+(sub2-1)*2^(iter-k):(sub2)*2^(iter-k))=offsety(1+(sub1-1)*2^(iter-k):(sub1)*2^(iter-k),1+(sub2-1)*2^(iter-k):(sub2)*2^(iter-k))+offy;
         end
           
         %apr�s la derni�re it�ration, faire le tracking pour chaque
         %fenetre �largie
         if k==iter
             
            
         % interpolation de n
         win=mi(k+1)+overlapTrack;
         ystart=posiy{k+1}(sub1,sub2)-win/2;
         xstart=posix{k+1}(sub1,sub2)-win/2;
         yend=posiy{k+1}(sub1,sub2)+win/2-1;
         xend=posix{k+1}(sub1,sub2)+win/2-1;
         if ystart<1, ystart=1;end
         if xstart<1, xstart=1;end
         if yend>M, yend=M;end
         if xend>N, xend=N;end
            xi=(1:win)-offsetx{k+1}(sub1,sub2)+xstart-1;
            yi=(1:win)-offsety{k+1}(sub1,sub2)+ystart-1;
            n_interp = interp2(double(n),xi',yi);
           s_win=s(ystart:yend,xstart:xend); 
            
         %tracking
              mauvais_track=0; 
       if strcmp(cl,'uint16')
       %n_win=uint16(n_interp);
       n_win= uint16(65535/(max(n_interp(:))-min(n_interp(:)))*(n_interp-min(n_interp(:))));
       s_win= uint16(65535/double(max(s_win(:))-min(s_win(:)))*double(s_win-min(s_win(:))));
       elseif strcmp(cl,'uint8')
       n_win= uint8(255/(max(n_interp(:))-min(n_interp(:)))*(n_interp-min(n_interp(:))));
       s_win= uint8(255/double(max(s_win(:))-min(s_win(:)))*double(s_win-min(s_win(:))));%n_win=uint8(n_interp);
       end
       %Tracking sur les petites fen?tres
        nM = feature2D(n_win,1,featsize,masscut);
 %figure(5),imshow(n_win,[]), hold on, plot(nM(:,1),nM(:,2),'go'), hold off           
        if numel(nM)==1||isempty(nM)
            figure(16)
            imshow(n_win,[])
            mauvais_track=true;
        else
        %Rejection process
        X= nM(:,5)>barcc;
        nM(X,1:5)=0;
        X= nM(:,4)>barrg;
        nM(X,1:5)=0;
        X= nM(:,3)./nM(:,4)<IdivRg;
        nM(X,1:5)=0;
        nM=nM(nM(:,1)~=0,:);
        if numel(nM)==1||isempty(nM)
            figure(16)
            imshow(n_win,[])
            mauvais_track=true;
        end
        end
        sM = feature2D(s_win,1,featsize,masscut);     
        if numel(sM)==1||isempty(sM)
            figure(16)
            imshow(s_win,[])
            mauvais_track=true;
        else
        %Rejection process
        X= sM(:,5)>barcc;
        sM(X,1:5)=0;
        X= sM(:,4)>barrg;
        sM(X,1:5)=0;
        X= sM(:,3)./sM(:,4)<IdivRg;
        sM(X,1:5)=0;
        sM=sM(sM(:,1)~=0,:);
        if numel(sM)==1||isempty(sM)
            figure(16)
            imshow(s_win,[])
            mauvais_track=true;
        end
        end
        if not(mauvais_track)
        Mtrack=[nM(:,1:2),ones(size(nM,1),1),ones(size(nM,1),1); sM(:,1:2),2*ones(size(sM,1),1),2*ones(size(sM,1),1)];
        [lub] = trackmem(Mtrack,maxd,2,2,0);        
        if lub==-1
            figure(15)
            ax2(1,2,1), imshow(s_win,[]), hold on, plot(sM(:,1),sM(:,2),'go'), hold off
            ax2(1,2,2), imshow(n_win,[]), hold on, plot(nM(:,1),nM(:,2),'go'), hold off
            mauvais_track=true;
        end
        end
        if mauvais_track
        mauvaises_fenetres=mauvaises_fenetres+1;
        else
        nb_part=max(lub(:,5));
      
        deplaceX(compteur+1:compteur+nb_part)=lub((lub(:,4)==2),1)-lub((lub(:,4)==1),1)+offsetx{k+1}(sub1,sub2)+restex;
         deplaceY(compteur+1:compteur+nb_part)=lub((lub(:,4)==2),2)-lub((lub(:,4)==1),2)+offsety{k+1}(sub1,sub2)+restey;
         positionX(compteur+1:compteur+nb_part)=lub((lub(:,4)==2),1)+(posix{k+1}(sub1,sub2)-win/2)-1;
         positionY(compteur+1:compteur+nb_part)=lub((lub(:,4)==2),2)+(posiy{k+1}(sub1,sub2)-win/2)-1;
        compteur=compteur+nb_part;
        end
         end
           end
       end
end   
end
disp([num2str(compteur-mauvaises_fenetres) ' features tracked']);
disp([num2str(mauvaises_fenetres) ' mauvaises fenetres']);
deplaceX(compteur+1:10000)=[];
deplaceY(compteur+1:10000)=[];
positionX(compteur+1:10000)=[];
positionY(compteur+1:10000)=[];
%remove duplicates
Mpos=[positionX positionY];
[A,ind]=sortrows(Mpos);
deplaceX=deplaceX(ind);
deplaceY=deplaceY(ind);
tolerance=2;
diff=A(1:end-1,:)-A(2:end,:);
duplicate=find((diff(:,1).^2+diff(:,2).^2)<tolerance);
A(duplicate,:)=[];
positionX=A(:,1);
positionY=A(:,2);
deplaceX(duplicate)=[];
deplaceY(duplicate)=[];
disp(['Particules trackees (non duplicates): ', num2str(length(deplaceX))])

%figure;
nb_beads = size(positionX);
deplaceXh = deplaceX;
deplaceYh = deplaceY;
for i = 1:nb_beads        
    index_neighbors = find(sqrt((positionX-positionX(i)).^2+(positionY-positionY(i)).^2)<r_neighbor);
    allbutone = find(index_neighbors~=i);
    index_neighbors = index_neighbors(allbutone);
    r_corr = mean((deplaceXh(i)*deplaceXh(index_neighbors)+deplaceYh(i)*deplaceYh(index_neighbors))./(sqrt(deplaceXh(i)^2+deplaceYh(i)^2)*sqrt(deplaceXh(index_neighbors).^2+deplaceYh(index_neighbors).^2)));
    if r_corr < corr_th
        deplaceX(i) = mean(deplaceX(index_neighbors));
        deplaceY(i) = mean(deplaceY(index_neighbors));
       % quiver(positionX,positionY,3*deplaceX,3*deplaceY,'r','AutoScale','off');hold on;
       % plot(positionX(i),positionY(i),'*');
    end
end

figure(h1) % Tracking
Idisp=zeros([size(s) 3]);
sn=double(s)/double(max(s(:)));
nn=double(n)/double(max(n(:)));
Idisp(:,:,1)=nn;
Idisp(:,:,2)=sn;
Idisp(:,:,3)=nn;
imshow(Idisp,[],'InitialMagnification','fit'),colormap(gray)
title(['Image  ',num2str(ff),' sur ', num2str(Nb_frames)]);
hold on
quiver(positionX,positionY,deplaceX,deplaceY,'r');
%hold on
% for i=6
%     for j=2
%         hold on
% plot([posx(1,i)-window/2 posx(1,i)-window/2 posx(1,i)+window/2 posx(1,i)+window/2 posx(1,i)-window/2],[posy(j,1)-window/2 posy(j,1)+window/2 posy(j,1)+window/2 posy(j,1)-window/2 posy(j,1)-window/2],'g')
%     end
% end
hold off
drawnow
%print(h1,fullfile(path,['cell',num2str(ff),'_tracking.tif']),'-dtiff','-r600'); 
%print(h1,fullfile(path,'beads',['image',num2str(ff),'.tif']),'-dtiff','-r300');

if ~exist(fullfile(path,'figureTFM'),'dir')
    mkdir(fullfile(path,'figureTFM'))
end

figurepath=cat(2,path,'\figureTFM\Tracking',num2str(ff),'.tif');
print(h1,fullfile(figurepath),'-dtiff','-r100');
A=imread(figurepath);
imwrite(A,cat(2,path,'\figureTFM\Tracking.tif'),'WriteMode','append');
delete (figurepath)


%% Displacements

% nodeX=floor(window*nbx/interval);
% nodeY=floor(window*nby/interval);
% [gridX,gridY]=meshgrid((0:nodeX-1)*interval+floor((dimx-nbx*window)/2)+1+interval/2,(0:nodeY-1)*interval+floor((dimy-nby*window)/2)+1+interval/2);
deplacementXa=pix*griddata(positionX,positionY,deplaceX,gridX,gridY);
deplacementYa=pix*griddata(positionX,positionY,deplaceY,gridX,gridY);
deplacementXa(isnan(deplacementXa))=0;
deplacementYa(isnan(deplacementYa))=0;


%smoothing
wsmoothing=2;
[x,y]=meshgrid((-4*wsmoothing:4*wsmoothing),(-4*wsmoothing:4*wsmoothing));
Gauss=exp(-2*(x.^2+y.^2)/wsmoothing^2);
h=Gauss/sum(Gauss(:));
deplacementX=imfilter(deplacementXa,h);
deplacementY=imfilter(deplacementYa,h);

deplacement = sqrt(deplacementX.^2+deplacementY.^2);

% end

%% Traction 
[Tractionx,Tractiony,mu,theta]=regularized(deplacementX,deplacementY,E,nu,pas,alphadef);
Tmagn=sqrt(Tractionx.^2+Tractiony.^2);


%% Calculs sous le doublet de cellules
disp('---- Calculs sur le doublet de cellules -----')

SE = strel('disk',1);
BWpattern=imdilate(mask_corr,SE);
%BWpattern=mask_corr;
BWgrid=logical(BWpattern(gridY(:,1),gridX(1,:)));

Txcrop=Tractionx.*BWgrid;
Tycrop=Tractiony.*BWgrid;

%calcul de l'?nergie contractile
u=Txcrop.*deplacementX+Tycrop.*deplacementY;  %produit scalaire traction par deplacement
U=pas^2/2*sum(u(:));
disp(['contractile energy: ',num2str(U), ' J']);
Ec2cells(ff)=U;


% image de traction map avec decoupe du pattern
Tmagn2=sqrt(Txcrop.^2+Tycrop.^2);
% figure('Name','Contour des traction 2')
% contourf(Tmagn2(end:-1:1,:)),colorbar

% contrainte moyenne sous le doublet
TmagnCROPMOY=sum(Tmagn2(:))/bwarea(BWgrid);
disp(['mean stress: ',  num2str(TmagnCROPMOY),' Pa']);
MOY2cells(ff)=TmagnCROPMOY;

%Traction totale sous les 2 cellules
Sgrid = pas^2;%surface d'une maille en metre au carr?
TmagnCROPTOT=sum(Tmagn2(:))*Sgrid;% somme des forces de traction en norme
disp(['Total des tractions en norme: ',  num2str(TmagnCROPTOT),' N']);
Total2cells(ff)=TmagnCROPTOT;

%Somme vectorielle (syst?me isol? ? l'?quilibre: devrait =0)
Vect2cells(ff)=sqrt((sum(Txcrop(:)))^2+(sum(Tycrop(:)))^2)*Sgrid;
unbalance=Vect2cells(ff)/TmagnCROPTOT;
disp(['Pourcentage de force non equilibree: ',  num2str(unbalance*100),' %']);

%Main contraction direction for the cell doublet (it is not valid to
%calculat it on individual cells since they are not at equilibrium)
%1.clean the mask

s=regionprops(mask_corr,'Area','Centroid');
[~,ind]=max([s.Area]);
 x0=s(ind).Centroid(:,1);
 y0=s(ind).Centroid(:,2);
Mxx=sum(sum((gridX-x0).*Txcrop))*pix*Sgrid^2;
Mxy=sum(sum((gridX-x0).*Tycrop))*pix*Sgrid^2;
Myx=sum(sum((gridY-y0).*Txcrop))*pix*Sgrid^2;
Myy=sum(sum((gridY-y0).*Tycrop))*pix*Sgrid^2;
M=[Mxx -Mxy;-Myx Myy];
Moment(ff)=trace(M);
[V,D]=eig(M);
if isreal(D)
    DD=max(-D);
%on identifie la direction principale (indice ind)
[Dmax, ind]=max(DD);
ang=180/pi*atan(real(V(2,ind))/real(V(1,ind)));
%angle par rapport la verticale
if ang>=0
Angle(ff)=90-ang;
else
Angle(ff)= -90-ang;
end
L1=D(ind,ind);
%autre direction
autre=mod(ind,2)+1;
L2=D(autre, autre);
else% non diagonalizable
    Mxyn=(Mxy+Myx)/2;
    Myxn=Mxyn;
    M=[Mxx -Mxyn;-Myxn Myy];
    [V,D]=eig(M);
DD=max(-D);
%on identifie la direction principale (indice ind)
[Dmax, ind]=max(DD);
ang=180/pi*atan(real(V(2,ind))/real(V(1,ind)));
%angle par rapport la verticale
if ang>=0
Angle(ff)=90-ang;
else
Angle(ff)= -90-ang;
end
L1=DD(ind);
%autre direction
autre=mod(ind,2)+1;
L2=DD(autre);
end
Polar(ff)=(L1-L2)/(L1+L2);
disp(['Angle of main contraction compared to the vertical: ',  num2str(Angle(ff))]);

disp(['Polarisation degree of the cell doublet (0=isotropic,1=uniaxial): ',  num2str(Polar(ff))]);

%% Calculs pour chaque cellule du doublet et force cell-cell
disp('---- Calculs sur les cellules individuelles -----')
BW=mask1_corr;
BWc1=imdilate(BW,SE);
BWgrid1=logical(BWc1(gridY(:,1),gridX(1,:)));
BW=mask2_corr;
BWc2=imdilate(BW,SE);
BWgrid2=logical(BWc2(gridY(:,1),gridX(1,:)));
Txcrop1=Tractionx.*BWgrid1;
Tycrop1=Tractiony.*BWgrid1;
Txcrop2=Tractionx.*BWgrid2;
Tycrop2=Tractiony.*BWgrid2;


Tnorm=sqrt(Txcrop1.^2+Tycrop1.^2);
MOYcell1(ff)=sum(Tnorm(:))/bwarea(BWgrid1);% contrainte moyenne
Totalcell1(ff)=sum(Tnorm(:))*Sgrid;% somme des forces de traction en norme
%Maximum stress
ntot=bwarea(BWgrid1);
thres=max(Tnorm(:));
while numel(find(Tnorm>=thres))<0.05*ntot
    thres=thres-20;
end
Tmaxcell1(ff)=mean(Tnorm(find(Tnorm>=thres)));

Tnorm=sqrt(Txcrop2.^2+Tycrop2.^2);
MOYcell2(ff)=sum(Tnorm(:))/bwarea(BWgrid2);
Totalcell2(ff)=sum(Tnorm(:))*Sgrid;% somme des forces de traction en norme
%Maximum stress
ntot=bwarea(BWgrid2);
thres=max(Tnorm(:));
while numel(find(Tnorm>=thres))<0.05*ntot
    thres=thres-20;
end
Tmaxcell2(ff)=mean(Tnorm(find(Tnorm>=thres)));

disp(['Contrainte moyenne sous cellule 1: ',  num2str(MOYcell1(ff)),' Pa']);
disp(['Contrainte moyenne sous cellule 2: ',  num2str(MOYcell2(ff)),' Pa']);

disp(['Contrainte max sous cellule 1: ',  num2str(Tmaxcell1(ff)),' Pa']);
disp(['Contrainte max sous cellule 2: ',  num2str(Tmaxcell2(ff)),' Pa']);

disp(['Total des tractions en norme sous cellule 1: ',  num2str(Totalcell1(ff)),' N']);
disp(['Total des tractions en norme sous cellule 2: ',  num2str(Totalcell2(ff)),' N']);

disp('---- Force cellule-cellule -----')
Fcellcell(ff)=sqrt((sum(Txcrop1(:))-sum(Txcrop2(:)))^2+(sum(Tycrop1(:))-sum(Tycrop2(:)))^2)*Sgrid/2;
disp(['Force cellule-cellule : ',  num2str(Fcellcell(ff)),' N']);
DFcellcell(ff)=sqrt((sum(Txcrop1(:))+sum(Txcrop2(:)))^2+(sum(Tycrop1(:))+sum(Tycrop2(:)))^2)*Sgrid/2;
disp(['Incertitude : ',  num2str(DFcellcell(ff)),' N']);
Cell_ECMratio1(ff)=sqrt((sum(Txcrop1(:)))^2+(sum(Tycrop1(:)))^2)*Sgrid/Totalcell1(ff);
Cell_ECMratio2(ff)=sqrt((sum(Txcrop2(:)))^2+(sum(Tycrop2(:)))^2)*Sgrid/Totalcell2(ff);
disp(['Rapport F cell-cell/F cell-ECM pour cellule 1: ',  num2str(Cell_ECMratio1(ff)*100),' %']);
disp(['Rapport F cell-cell/F cell-ECM pour cellule 2: ',  num2str(Cell_ECMratio2(ff)*100),' %']);
%Calcul des moments
s=regionprops(mask1_corr,'Area','Centroid');
[~,ind]=max([s.Area]);
 x1=s(ind).Centroid(:,1);
 y1=s(ind).Centroid(:,2);
s=regionprops(mask2_corr,'Area','Centroid');
[~,ind]=max([s.Area]);
 x2=s(ind).Centroid(:,1);
 y2=s(ind).Centroid(:,2);
 Mxx1=sum(sum((gridX-x1).*Txcrop1))*pix*Sgrid;
Mxx2=sum(sum((gridX-x2).*Txcrop2))*pix*Sgrid;
Myy1=sum(sum((gridY-y1).*Tycrop1))*pix*Sgrid;
Myy2=sum(sum((gridY-y2).*Tycrop2))*pix*Sgrid;
Mxy1=sum(sum((gridX-x1).*Tycrop1))*pix*Sgrid;
Mxy2=sum(sum((gridX-x1).*Tycrop2))*pix*Sgrid;
Myx1=sum(sum((gridY-y1).*Txcrop1))*pix*Sgrid;
Myx2=sum(sum((gridY-y2).*Txcrop2))*pix*Sgrid;
M1=[Mxx1 -Mxy1;-Myx1 Myy1];
M2=[Mxx2 -Mxy2;-Myx2 Myy2];

%TBD
Momentcell1(ff)=trace(M1);
Momentcell2(ff)=trace(M2);
disp(['contractile moment cell 1: ',  num2str(Momentcell1(ff)),' N.m']);
disp(['contractile moment cell 2: ',  num2str(Momentcell2(ff)),' N.m']);

Forcex_cell1(ff)=sum(abs(Txcrop1(:)))*Sgrid;%sum of x component of forces (in magnitude) for cell 1
Forcey_cell1(ff)=sum(abs(Tycrop1(:)))*Sgrid;%sum of y component of forces (in magnitude) for cell 1
Forcex_cell2(ff)=sum(abs(Txcrop2(:)))*Sgrid;%sum of x component of forces (in magnitude) for cell 2
Forcey_cell2(ff)=sum(abs(Tycrop2(:)))*Sgrid;%sum of y component of forces (in magnitude) for cell 2
Deplcropx1=deplacementX.*BWgrid1;
Deplcropy1=deplacementY.*BWgrid1;
Deplcropx2=deplacementX.*BWgrid2;
Deplcropy2=deplacementY.*BWgrid2;
Deplx_cell1(ff)=sum(abs(Deplcropx1(:)));%sum of x component of displacements (in magnitude) for cell 1
Deply_cell1(ff)=sum(abs(Deplcropy1(:)));%sum of x component of displacements (in magnitude) for cell 1
Deplx_cell2(ff)=sum(abs(Deplcropx2(:)));%sum of x component of displacements (in magnitude) for cell 1
Deply_cell2(ff)=sum(abs(Deplcropy2(:)));%sum of x component of displacements (in magnitude) for cell 1

    
%% Display

figure(h4) % Displacements
imagesc(deplacement),colormap(blackcmap),colorbar,axis off, daspect([1 1 1]),
title(['Displacement - Frame ',num2str(ff)]);%.*BWgrid, caxis([0 max_disp]), 
%save(fullfile(path,['cell',num2str(ff),'_Displacement.dat']),'deplacement','-ascii','-tabs')

figure(h3)
contourf(Tmagn(end:-1:1,:)),colormap(jet),colorbar%,colormap(jet)
axis equal
set(gca,'Visible','off');
hold on
[c,h]=contour(BWgrid(end:-1:1,:),[0.5 0.5]);
set(h,'LineWidth',1,'LineColor','w')
hold on
title(['Traction - Frame ',num2str(ff)]);
hold off

%
figure(h2)% Forces on BF image
imshow(bf,[],'InitialMagnification','fit')
title(['Image ', num2str(ff)])
hold on
quiver(gridX(1:2:end,1:2:end),gridY(1:2:end,1:2:end),8*Tractionx(1:2:end,1:2:end),8*Tractiony(1:2:end,1:2:end),2,'r','AutoScale','off')

figurepath=cat(2,path,'\figureTFM\BF_polar_degree',num2str(ff),'.tif');
print(h2,fullfile(figurepath),'-dtiff','-r100');
A=imread(figurepath);
imwrite(A,cat(2,path,'\figureTFM\BF_polar_degree.tif'),'WriteMode','append');
delete (figurepath)
%%
%save (fullfile(path,[filename(1:end-4),'_results', num2str(ff),'.mat']),'p_corr','BW','BWgrid','deplacement*','Traction*','eigenvalue*');
%cellule(1:size(p_corr,1),1:size(p_corr,2),ff)=p_corr;
Brightfield(:,:,ff)=uint8(round(double(bf-min(bf(:)))/double(max(bf(:))-min(bf(:)))*255));
Mask1(:,:,ff)=mask1_corr;
Mask2(:,:,ff)=mask2_corr;
Dx(:,:,ff)=deplacementX;
Dy(:,:,ff)=deplacementY;
Tx(:,:,ff)=Tractionx;
Ty(:,:,ff)=Tractiony;

%caculate contractile Energy in fJ
EcfJ= Ec2cells.*10^15;

%calculate mean values
meanEc = mean(Ec2cells, 'all'); 
meanEcfJ = mean(EcfJ, 'all');
meanPolar = mean(Polar, 'all');
meanMOY2cells = mean(MOY2cells, 'all');
meanMoment = mean(Moment, 'all');
end
%%
save (fullfile(path,'Allresults2.mat'),'Mask*','Dx','Dy','Tx','Ty','grid*','Brightfield');

%Save calculated parameters
fid=fopen(fullfile(path,'Resultats2cells.txt'),'w+');
fprintf(fid,['Ec \t Pmoy \t Ttot \t  T_imbalance \t Fcell-cell \t DFcell-cell \t Moment \t Angle \t Polar degree \t Pmoy-cell1 \t Pmax_cell1 \t' ...
'Pmoy-cell2 \t Pmax_cell2 \t T-cell1 \t Ttot-cell2 \t Rcell1-ECM \t Rcell2-ECM \t Mmt-cell1 \t Mmt-cell2 \t Force x cell1 \t Force y cell1 \t' ...
'Force x cell2 \t Force y cell2 \t Displ x cell1 \t Displ y cell1 \t Displ x cell2 \t Displ y cell2 \n']);
fclose(fid);
dlmwrite(fullfile(path,'Resultats2cells.txt'),[Ec2cells',MOY2cells',Total2cells',Vect2cells',Fcellcell',DFcellcell',Moment',Angle',Polar', ...
    MOYcell1',Tmaxcell1',MOYcell2',Tmaxcell2',Totalcell1',Totalcell2',Cell_ECMratio1',Cell_ECMratio2',Momentcell1',Momentcell2',Forcex_cell1',...
    Forcey_cell1',Forcex_cell2',Forcey_cell2',Deplx_cell1',Deply_cell1',Deplx_cell2',Deply_cell2'],'-append','delimiter','\t')
disp('Results saved in: ')
disp(fullfile(path,'Resultats2cells.txt'))


% Export to Excel
T = table(Ec2cells', EcfJ', MOY2cells', Moment', Polar', (Vect2cells./Total2cells*100)', 'VariableNames', {'contractile energy (J)','contractile energy (fJ)','average stress (Pa)','Moment','Polar','Out of Equilibrium'});
X = table(meanEc, meanEcfJ, meanMOY2cells ,meanMoment, meanPolar );
writetable(T,fullfile(path,'results.xlsx'),'Sheet',1,'Range','A1','WriteVariableNames', true);
writetable(X,fullfile(path,'results.xlsx'),'Sheet',2,'Range','A1');

%graphs
pos=(1:Nb_frames);
graphs=figure('units','Normalized','Position', [0.01,0.03,0.99,0.9],'Name','Doublet de cellules');
ax2(2,3,1),plot(pos,MOY2cells,'o-'),title('average stress (Pa)');
ax2(2,3,2), plot(pos,Momentcell1,'-og'), hold on, plot(pos,Momentcell2,'-om'),hold off,title('contractile moment(Nm)')
ax2(2,3,3),plot(pos,Vect2cells./Total2cells*100,'d-g'),title('Out of equilibrium');
ax2(2,3,4), plot(pos,Ec2cells, '--o'),title('contractile energy  (J)');
ax2(2,3,5), plot(pos,Polar,'-sc'),  title('Degre of polarisation')

print(graphs,cat(2,path,'graphs1'),'-dpng');

save (fullfile(path,[filename(1:end-4),'_param','.mat']),'grid*','interval','pas');
fid=fopen(fullfile(path,[filename,'_Param.txt']),'w');
fprintf(fid,'Registration (x direction): %g \n',regx);
fprintf(fid,'Registration (y direction): %g \n',regy);
fprintf(fid,'PIV window size: %g \n \n',window);
fprintf(fid,'Number of PIV iterations (added to the 1st): %g \n \n',iter);
fprintf(fid,'Overlap PIV: %g \n \n',overlapPIV);
fprintf(fid,'Overlap tracking: %g \n \n',overlapTrack);

fprintf(fid,'Tracking parameters: \n');
fprintf(fid,'Feature size: %g \n',featsize);
fprintf(fid,'Min intensity(masscut): %g\n',masscut);
fprintf(fid,'Max radius: %g\n',barrg);
fprintf(fid,'Max eccentricity: %g\n',barcc);
fprintf(fid,'Min ratio intensity/radius: %4.1f\n',IdivRg);
fprintf(fid,'Max displacement in pixels: %g\n\n',maxd);
fprintf(fid,'Pixel size (in meter): %g\n',pix);
fprintf(fid,'Data interval (for displacement and force) in pixels: %g\n',interval);
fprintf(fid,'Young modulus (in Pa): %g\n',E);
fprintf(fid,'Poisson ratio: %g\n',nu);
fprintf(fid,'Regularization parameter: %g\n',alphadef);
fclose(fid);

end