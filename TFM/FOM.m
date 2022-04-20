function ForceOrientationMap = FOM(path)
    load(fullfile(path, 'Allresults2.mat'))
    size_tension_tensor = size(Ty(1,1,:));
    nb_frames = size_tension_tensor(3);
    %nb_frames = 1;
    for i=1:nb_frames    
        close all;
        Angle=atan(abs(Ty(:,:,i)./Tx(:,:,i)));
        TMagn=sqrt(Tx(:,:,i).^2+Ty(:,:,i).^2);


        AmpT=TMagn/(max(TMagn(:))-50); % this adjusts the intensity. There is a saturation with a lot of values at 1 on purpose, since it is nicer to look at.
        ind=find(AmpT>1);
        AmpT(ind)=1;

        angle_255=floor(Angle/max(Angle(:))*255)+1;         % Norm the angle to be between 1 and 256

        DIMS=size(angle_255);
        cm=jet(256);

        R=reshape(cm(angle_255,1),DIMS);                    % Convert to RGB colors
        G=reshape(cm(angle_255,2),DIMS);
        B=reshape(cm(angle_255,3),DIMS);

        RGB=zeros(DIMS(1),DIMS(2),3);
        RGB(:,:,1)=R.*AmpT;
        RGB(:,:,2)=G.*AmpT;
        RGB(:,:,3)=B.*AmpT;                                         
           
        h=figure;
        imshow(RGB,[],'InitialMagnification','fit')
        colormap(jet)
        colorbar;  caxis([0 90])% Display FOM
        c = colorbar;
        c.Label.String = 'Angle of traction force in degree';
        axis equal
        set(gca,'Visible','off');
        title(['frame: ',num2str(i)]);
        
        figurepath=cat(2,path,'figureTFM\FOM',num2str(i),'.tif');   % Save Film
        print(h,fullfile(figurepath),'-dtiff','-r100');
        A=imread(figurepath);
        imwrite(A,cat(2,path,'figureTFM\FOM.tif'),'WriteMode','append');
        delete (figurepath)
    end
end