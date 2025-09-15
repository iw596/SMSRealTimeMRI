
% Takes matrix data [x,y,..,coils] and applies
% covaraince matrix to it
function res = ApplyCovarianceMatrix(data,covMat)
 orig_size = size(data);
    ncoils = size(covMat,1);
    nelements = prod(orig_size)/ncoils;
    
    if ~(ncoils == orig_size(end))
        error('Number of coils in decorrelation matrix does not match the data');
    end
    
    data = permute(reshape(data,nelements,ncoils),[2 1]);
    res = reshape(permute(covMat*data,[2 1]),orig_size);
end