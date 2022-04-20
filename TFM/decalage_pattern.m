function [deltax,deltay,Cmax]=decalage_pattern(im1,im2)
%Attention!! la taille des 2 images doivent être identiques
Nx=size(im1,2);
Ny=size(im1,1);
if rem(Nx,2)
    Nx=Nx-1;
    im1=im1(1:Ny,1:Nx);
    im2=im2(1:Ny,1:Nx);   
end
if rem(Ny,2)
     Ny=Ny-1;
     im1=im1(1:Ny,1:Nx);
    im2=im2(1:Ny,1:Nx);
end

A1=real(ifft2(fft2(im1) .* fft2(rot90(im1,2))));
A2=real(ifft2(fft2(im2) .* fft2(rot90(im2,2))));
M1=max(A1(:));
M2=max(A2(:));
C = real(ifft2(fft2(im1) .* fft2(rot90(im2,2))));
C2=fftshift(C);
%imagesc(C2),colorbar
[Vmax,imax]=max(C2(:));
Cmax=Vmax/sqrt(M1*M2);
[I,J]=ind2sub(size(C),imax);
while (abs(I-Ny/2)>Ny/4)||(abs(J-Nx/2)>Nx/4)
    left=min(I-1,30);
    right=min(Ny-I,30);
    left2=min(J-1,30);
    right2=min(Nx-J,30);
    C2(I-left:I+right,J-left2:J+right2)=0;
    [Vmax,imax]=max(C2(:));
    Cmax=Vmax/sqrt(M1*M2);
    [I,J]=ind2sub(size(C),imax);
disp('Attention! il existe un max de correlation très décalé');
end
%figure,imagesc(C2),colorbar
if (I==size(C,1)) || (I==1)
    Isub=I;
else
%mesure décalage précision subpixel
    Isub=I+(log(C2(I-1,J))-log(C2(I+1,J)))/(2*(log(C2(I-1,J))+log(C2(I+1,J))-2*log(C2(I,J))));
end
if (J==size(C,2)) || (J==1)
    Jsub=J;
else
Jsub=J+(log(C2(I,J-1))-log(C2(I,J+1)))/(2*(log(C2(I,J-1))+log(C2(I,J+1))-2*log(C2(I,J))));
end
deltax=Jsub-Nx/2; 
deltay=Isub-Ny/2;
% si deltay>0 im2 décalée vers le haut par rapport à im1
% si deltax>0 im2 décalée vers la gauche par rapport à im1