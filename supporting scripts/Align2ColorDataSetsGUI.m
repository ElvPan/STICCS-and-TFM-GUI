function [imch2]=Align2ColorDataSetsGUI(ax,data1,data2)

% aligning images from channel 2 w.r.t. images in channel 1

[optimizer, metric] = imregconfig('multimodal');
optimizer.InitialRadius = 0.009;
optimizer.Epsilon = 1.5e-4;
optimizer.GrowthFactor = 1.01;
optimizer.MaximumIterations = 300;
% first find the most common (average) transformation between two channels

ave1=mean(data1,3);
ave2=mean(data2,3);
tform = imregtform(ave2,ave1, 'rigid', optimizer, metric);

cla(ax)
ylim(ax,[0,1])
xlim(ax,[0,1])
ph = patch(ax,[0 0 0 0],[0 0 1 1],[0.67578 1 0.18359]); %greenyellow
th = text(ax,1,1,'Aligning Channels...0%','VerticalAlignment','bottom','HorizontalAlignment','right');
for i=1:size(data2,3)
    fixed=double(data1(:,:,i));
   moving=double(data2(:,:,i));
    %tform = imregtform(moving, fixed, 'rigid', optimizer, metric);
   movingRegistered = imwarp(moving,tform,'OutputView',imref2d(size(fixed)));
   imch2(:,:,i)=movingRegistered;
   ph.XData = [0 i/size(data2,3)  i/size(data2,3) 0];
   th.String = sprintf('Aligning Channels...%.0f%%',round(i/size(data2,3)*100));
   drawnow %update graphics
end