function [invconv,h2]=smoothSeries(data,sigmafilt,type,posneg);

%define filter
if strcmp(type,'log')
h2 = fspecial('log', [size(data,1) size(data,2)], sigmafilt);
elseif strcmp(type,'g')
h2 = fspecial('gaussian', [size(data,1) size(data,2)], sigmafilt);
elseif strcmp(type,'lap')
h2 = fspecial('laplacian', sigmafilt);
elseif strcmp(type,'disk')
h2 = fspecial('disk', sigmafilt);
elseif strcmp(type,'average')
 h2 = fspecial('average', sigmafilt);
end

if strcmp(posneg,'n')
    h2=-h2;
end

if strcmp(type,'log') || strcmp(type,'g')
 ffth2=fft2(h2);
for i=1:size(data,3)
fftim=fft2(data(:,:,i));
fftprod=fftim.*(ffth2);
invconv(:,:,i)=fftshift(ifft2(fftprod));
% if strcmp(type,'log')
% invconv(:,:,i)=invconv(:,:,i).*(invconv(:,:,i)>0)+10^-12;
% end
end
else
  invconv=imfilter(data,h2);  
end
