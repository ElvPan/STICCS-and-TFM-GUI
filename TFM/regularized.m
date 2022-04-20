function [Tractionx,Tractiony,mu,theta]=regularized(deplacementx,deplacementy,E,nu,ech,alpha)

%ech=interval*pix;
[M,N]=size(deplacementx);
%padding with zeros
m=2;
while (2^m < M)||(2^m < N )
    m=m+1;
end
M2=2^m;
N2=M2;
% n=2;
% while (2^n < N)
%     n=n+1;

% end
% N2=2^n;

uxt=fft2(deplacementx,M2,N2);
uyt=fft2(deplacementy,M2,N2);
%suppression de la translation
uxt(1,1)=0;
uyt(1,1)=0;

Kx=[(2*pi/ech)/N2*[0:N2/2],-(2*pi/ech)/N2*(N2-[N2/2+1:N2-1])];
Ky=[(2*pi/ech)/M2*[0:M2/2],-(2*pi/ech)/M2*(M2-[M2/2+1:M2-1])];
[kx, ky]=meshgrid(Kx,Ky);
k=sqrt(kx.^2+ky.^2);
Txt=zeros(M2,N2);
Tyt=zeros(M2,N2);

for i=1:M2
    for j=1:N2
        if (i==M2/2+1)||(j==N2/2+1)  %Fréquences de Nyquist
            Gt=2*(1+nu)/(E*k(i,j)^3)*[(1-nu)*k(i,j)^2+nu*ky(i,j)^2  0; 0 (1-nu)*k(i,j)^2+nu*kx(i,j)^2] ;
            Tt=(Gt'*Gt+alpha*eye(2))^-1*Gt'*[uxt(i,j); uyt(i,j)];
            Txt(i,j)=Tt(1);
            Tyt(i,j)=Tt(2);
        elseif ~((i==1) && (j==1))
        Gt=2*(1+nu)/(E*k(i,j)^3)*[(1-nu)*k(i,j)^2+nu*ky(i,j)^2  -nu*kx(i,j)*ky(i,j); -nu*kx(i,j)*ky(i,j) (1-nu)*k(i,j)^2+nu*kx(i,j)^2] ;
        Tt=(Gt'*Gt+alpha*eye(2))^-1*Gt'*[uxt(i,j) ; uyt(i,j)]; 
        Txt(i,j)=Tt(1);
        Tyt(i,j)=Tt(2);
       
        end
    end
end

        
Tx=ifft2(Txt);
Ty=ifft2(Tyt);
Tractionx=real(Tx([1:M],[1:N]));
Tractiony=real(Ty([1:M],[1:N]));

%calcul des moments contractile(selon Butler)
Mxx=-i*(Txt(1,2)+Txt(1,N2))/(2*kx(1,2));
Myy=-i*(Tyt(2,1)+Tyt(M2,1))/(2*ky(2,1));
Mxy=-i/2*((Tyt(1,2)+Tyt(1,N))/(2*kx(1,2))+(Txt(2,1)+Txt(M,1))/(2*ky(2,1)));
%net contractile moment
mu=real(Mxx+Myy);
%angle de l'axe principal en deg
theta=real(180/pi*atan(2*Mxy/(Mxx-Myy))/2);
