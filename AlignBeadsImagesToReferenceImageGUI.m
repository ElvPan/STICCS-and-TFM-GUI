function [alignedSeries,displaceX,displaceY]=AlignBeadsImagesToReferenceImageGUI(ax,ref,series,mask,CtoCdist,MinInt,SzFeat);
%% find boundary of mask and define which points aro outside (good to use)...
[B,L] = bwboundaries(mask,'noholes');
%% finding beads in each frame and in ref image
clear xcm2 ycm2 outs2
% ref channel
out=pkfnd(ref,MinInt,SzFeat);
in = inpolygon(out(:,1),out(:,2),B{1,1}(:,1),B{1,1}(:,2));
outInd=find(~in);
outB=out(outInd,:);
% now series
%%
cla(ax)
ylim(ax,[0,1])
xlim(ax,[0,1])
ph = patch(ax,[0 0 0 0],[0 0 1 1],[0.67578 1 0.18359]); %greenyellow
th = text(ax,1,1,'Aligning Images To A Reference Image...0%','VerticalAlignment','bottom','HorizontalAlignment','right');

for i=1:size(series,3)
out2=pkfnd(squeeze(series(:,:,i)),MinInt,SzFeat);
in = inpolygon(out2(:,1),out2(:,2),B{1,1}(:,1),B{1,1}(:,2));
out2Ind=find(~in);
out2B=out2(out2Ind,:);
D12=pdist2(outB,out2B);
pairs=0;
clear id12
for j=1:size(D12,1)
    firstInd=find(D12(j,:)==min(D12(j,:)));
    id12(j)=firstInd(1);
    if D12(j,id12(j))>CtoCdist
        id12(j)=NaN;
    else
        pairs=pairs+1;
    end
end

indgood=find(~isnan(id12));
ind1=id12(indgood);
clear dist dispx dispy
for j=1:length(indgood)
    dist(j)=D12(indgood(j),ind1(j));
    dispx(j)=out2B(ind1(j),1)-outB(indgood(j),1);
    dispy(j)=out2B(ind1(j),2)-outB(indgood(j),2);
end
displaceX(i)=median(dispx);
displaceY(i)=median(dispy);
% now define rigid transform and translate data 
trans=[-displaceX(i) -displaceY(i)];
tforms = rigid2d([1 0 ; 0 1],trans);
Rmoving1 = imref2d(size(squeeze(series(:,:,i))),1,1);
RegisteredVolume = imwarp(squeeze(series(:,:,i)),Rmoving1,tforms,'bicubic','OutputView',Rmoving1);
alignedSeries(:,:,i)=RegisteredVolume;
%i
ph.XData = [0 i/size(series,3)  i/size(series,3) 0];
    th.String = sprintf('Aligning Images To A Reference Image...%.0f%%',round(i/size(series,3)*100));
    drawnow %update graphics
end





