function res = ApplyPCACompression(kdata,datasize,NCha,compressionMatrix)
    if (ndims(kdata) ~=4)
    kdata = reshape(kdata,[prod(datasize), NCha]);
    res = kdata * compressionMatrix;
    res = reshape(res,[datasize,size(res,2)]);
    else
        NRep = size(kdata,4);
        res = zeros([datasize size(compressionMatrix,2) NRep]);
        for r = 1:NRep
            tmp = reshape(kdata(:,:,:,r),[prod(datasize), NCha]);
            tmp = tmp * compressionMatrix;
            res(:,:,:,r) = reshape(tmp,[datasize,size(tmp,2)]);
        end
    end
end