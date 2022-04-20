function [data,meta]=readOMEAlternativeGUI(ax2,fullname);

meta = GetOMEData(fullname);
reader = bfGetReader(fullname);

iPlane = reader.getIndex(0,0,0) + 1;
im=single(bfGetPlane(reader, iPlane));
data=zeros(size(im,1),size(im,2),meta.SizeZ,meta.SizeT,meta.SizeC);

cla(ax2)
ylim(ax2,[0,1])
xlim(ax2,[0,1])
ph = patch(ax2,[0 0 0 0],[0 0 1 1],[0.67578 1 0.18359]); %greenyellow
th = text(ax2,1,1,'Loading Data...0%','VerticalAlignment','bottom','HorizontalAlignment','right');
totframes=meta.SizeZ*meta.SizeT*meta.SizeC;
j=1;
for zload=1: meta.SizeZ
    for tload=1: meta.SizeT
        for cload=1: meta.SizeC
            iPlane = reader.getIndex(zload - 1, cload -1,tload - 1) + 1;
            im=single(bfGetPlane(reader, iPlane));
            data(:,:,zload,tload,cload) = im;
            ph.XData = [0 j/totframes  j/totframes 0];
            th.String = sprintf('Loading Data...%.0f%%',round(j/totframes*100));
            drawnow %update graphics
            j=j+1;
        end
    end
end

data=squeeze(data);
