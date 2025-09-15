function M = PrepVideo(data,pTile,type)

    if (nargin < 3)
        type = "MPEG-4"
    end
    scale   =   prctile(abs(data(:)), pTile);
    data    =   double(abs(data)/scale);
    data(abs(data)>1)   =   1;
    
    for t = 1:size(data, 3)
        if (type ~= "MPEG-4")
            M(t)=im2frame(repmat(data(:,:,t),1,1,3));
            M(t).cdata = rgb2gray(M(t).cdata);
        else
             M(t)=im2frame(repmat(data(:,:,t),1,1,3));
        end
    
    end
end