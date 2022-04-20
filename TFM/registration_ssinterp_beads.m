function [n_corr,s_corr,regx,regy]=registration_ssinterp_beads(nonstressed,stressed)
%corrige le décalge au pixel 2images: billes stressées, non
%stressées, 
%soustraction de background
background = imopen(stressed,strel('disk',12));
stressed2=imsubtract(stressed,background);
background = imopen(nonstressed,strel('disk',12));
nonstressed2=imsubtract(nonstressed,background);
%augmentation du contraste
n=contrast(nonstressed2,0.1,0.1);
s=contrast(stressed2,0.1,0.1);
[regx,regy]=decalage(n, s);

deltax = round(regx);
deltay = round(regy);

%correction du décalage
if (deltax>0) && (deltay>0)
        n_corr=n(deltay+1:size(nonstressed,1),deltax+1:size(nonstressed,2));
        s_corr=s(1:size(stressed,1)-deltay,1:size(stressed,2)-deltax);
    elseif (deltax>0) && (deltay<=0)
        n_corr=n(1:size(nonstressed,1)-abs(deltay),deltax+1:size(nonstressed,2));
        s_corr=s(abs(deltay)+1:size(stressed,1),1:size(stressed,2)-deltax);
    elseif (deltax<=0) && (deltay>0)
        n_corr=n(deltay+1:size(nonstressed,1),1:size(nonstressed,2)-abs(deltax));
        s_corr=s(1:size(stressed,1)-deltay,abs(deltax)+1:size(stressed,2));
    else
        s_corr=s(abs(deltay)+1:size(stressed,1),abs(deltax)+1:size(stressed,2));
        n_corr=n(1:size(nonstressed,1)-abs(deltay),1:size(nonstressed,2)-abs(deltax));
end