function [timecorrnorm] = stics_byfft(imgser,aveim,tauLimit);
%tic 
% August 23 2010
% Elvis Pandzic
% Calculates the full time correlation function given 3D array of image series
sizet=size(imgser,3);
sizey=size(imgser,2);
sizex=size(imgser,1);
if nargin==2
fractau=0.1;
tauLimit=round(fractau*sizet);
end
%tic

fftimser=fftn(imgser,[size(imgser,1) size(imgser,2) size(imgser,3)]);
clear imgser
timecorrnorm=double(zeros(sizex,sizey,tauLimit));
corrfn=conj(fftimser).*fftimser;clear fftimser;
timecorr=ifft(corrfn,[],3);clear corrfn


for i=1:tauLimit
    timecorrnorm(:,:,i)=abs(fftshift(ifft2(timecorr(:,:,i))))/(sizex*sizey*sizet*(aveim^2))-1;
end

%toc