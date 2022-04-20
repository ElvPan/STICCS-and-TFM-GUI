function [data,meta]=readOMEAlternative(fullname);

meta = GetOMEData(fullname);
reader = bfGetReader(fullname);

iPlane = reader.getIndex(0,0,0) + 1;
im=single(bfGetPlane(reader, iPlane));
data=zeros(size(im,1),size(im,2),meta.SizeZ,meta.SizeT,meta.SizeC);



for zload=1: meta.SizeZ 
    for tload=1: meta.SizeT
        for cload=1: meta.SizeC 
      iPlane = reader.getIndex(zload - 1, cload -1,tload - 1) + 1;
      im=single(bfGetPlane(reader, iPlane));
      data(:,:,zload,tload,cload) = im;
        end
    end
end

data=squeeze(data);
