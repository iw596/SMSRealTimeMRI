
% Function which takes a note that the size of the dimension of the data must be divisible
% by the binsize
function [binnedData] = binData(data,binSize,dim)
% First check if data (at dim) is divisible by binsize
if (mod(size(data,dim),binSize) ~=0 )
    error("size of data is not divisible by binSize")
end
NBins = size(data,dim)./binSize;
NDims = ndims(data);
NewNDims = NDims + 1;

dimSize = size(data);


binnedData = reshape()





end