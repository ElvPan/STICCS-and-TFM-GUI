function [BW3]=binarize_pattern_isodata(t_corr)
pattern=contrast(t_corr,0.5,0.5);
cl=class(t_corr);
%pattern=pattern=double(intmax(cl))/double(max(t_corr(:))-min(t_corr(:)))*double(t_corr-min(t_corr(:)));;
%imshow(pattern,[]);
th=intmax(cl)/2;
pix_high=pattern(pattern>th);
pix_low=pattern((pattern<th)|(pattern==th));
mean_high=mean(pix_high);
mean_low=mean(pix_low);
av=(mean_high+mean_low)/2;
while abs(th-av)>200;
    th=av;
    pix_high=pattern(pattern>th);
    pix_low=pattern((pattern<th)|(pattern==th));
    mean_high=mean(pix_high);
    mean_low=mean(pix_low);
    av=(mean_high+mean_low)/2;

end
patternBIN=pattern>th;%+6000;
%figure(1),imshow(patternBIN,[],'InitialMagnification','fit')


%% Amélioration de l'image binaire
%élimine les zones blanches sur le bord
BW1 = imclearborder(patternBIN);
%figure, imshow(BINnobord), title('cleared border image');
%figure, imshow(BW1);
SE = strel('disk',8);
BW=imdilate(BW1,SE);
%figure, imshow(BW);
SE = strel('disk',5);
BW2= imerode(BW,SE);
%figure, imshow(BW2);

%élimine les zones noires au milieu du pattern
BW3 = imfill(BW2, 'holes');
CC=bwconncomp(BW3);
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);
for i=1:CC.NumObjects
    if i~=idx
BW3(CC.PixelIdxList{i}) = 0;
    end
end
%figure, imshow(BW3);

% figure('Name','Image BINARISE')
% imagesc(BW3*50+pattern),colormap(gray)



