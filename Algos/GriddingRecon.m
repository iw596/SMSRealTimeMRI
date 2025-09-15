%% Function to perform coil by coil gridding of data using adjoint FT (passed in as an anonymous function)

function res = GriddingRecon(data,FTH,imsize)
    
    if (ndims(data) == 4)
        % If there are multiple time frames use a for loop to go through them
        Nt = size(data,4);
        res = zeros([imsize size(data,3) Nt]);
        for t = 1:Nt
            res(:,:,:,t) = FTH(data(:,:,:,t));
        end
    else
        res = FTH(data(:,:,:));
    end
end