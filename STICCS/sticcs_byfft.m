open sticsfunction [timecorrnorm] = sticcs_byfft(imgser1,imgser2,aveim1,aveim2,tauLimit);

%July 2015
% Elvis Pandzic
% Calculates the full time correlation function given 3D array of image series
%tic
sizet=size(imgser1,3);
sizey=size(imgser1,2);
sizex=size(imgser1,1);

fftimser1=fftn(imgser1,[size(imgser1,1) size(imgser1,2) size(imgser1,3)]);clear imgser1
fftimser2=fftn(imgser2,[size(imgser2,1) size(imgser2,2) size(imgser2,3)]);clear imgser2

timecorrnorm=double(zeros(sizex,sizey,tauLimit));
corrfn1=conj(fftimser1).*fftimser1;
corrfn2=conj(fftimser2).*fftimser2;
corrfn3=conj(fftimser1).*fftimser2;
corrfn4=conj(fftimser2).*fftimser1;
clear fftimser1 fftimser2;
timecorr3=ifft(corrfn3,[],3);clear corrfn3
timecorr4=ifft(corrfn4,[],3);clear corrfn4
timecorr1=ifft(corrfn1,[],3);clear corrfn1
timecorr2=ifft(corrfn2,[],3);clear corrfn2


for i=1:tauLimit
% 
  timecorrnorm(:,:,i,1)=abs(fftshift(ifft2(timecorr1(:,:,i))))/(sizex*sizey*sizet*(aveim1^2))-1;
  timecorrnorm(:,:,i,2)=abs(fftshift(ifft2(timecorr2(:,:,i))))/(sizex*sizey*sizet*(aveim2^2))-1;
  timecorrnorm(:,:,i,3)=abs(fftshift(ifft2(timecorr3(:,:,i))))/(sizex*sizey*sizet*aveim1*aveim2)-1;
  timecorrnorm(:,:,i,4)=abs(fftshift(ifft2(timecorr4(:,:,i))))/(sizex*sizey*sizet*aveim1*aveim2)-1;
% 
%   timecorrnorm(:,:,i,1)=timecorrnorm(:,:,i,1)-min(min(timecorrnorm(:,:,i,1)))+10^-10;
%   timecorrnorm(:,:,i,2)=timecorrnorm(:,:,i,2)-min(min(timecorrnorm(:,:,i,2)))+10^-10;
%   timecorrnorm(:,:,i,3)=timecorrnorm(:,:,i,3)-min(min(timecorrnorm(:,:,i,3)))+10^-10;
%   timecorrnorm(:,:,i,4)=timecorrnorm(:,:,i,4)-min(min(timecorrnorm(:,:,i,4)))+10^-10;
%   
% timecorrnorm(:,:,i,1)=abs(fftshift( ifft2(timecorr1(:,:,i)./timecorr1(:,:,1) )));
%   timecorrnorm(:,:,i,2)=abs(fftshift( ifft2(timecorr2(:,:,i)./timecorr2(:,:,1) )));
%   timecorrnorm(:,:,i,3)=abs(fftshift( ifft2(timecorr3(:,:,i)./(timecorr1(:,:,1).*timecorr2(:,:,1)).^0.5 )));;
%   timecorrnorm(:,:,i,4)=abs(fftshift( ifft2(timecorr4(:,:,i)./(timecorr1(:,:,1).*timecorr2(:,:,1)).^0.5 )));;

end




%toc