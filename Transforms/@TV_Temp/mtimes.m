
% This function calculates the temporal finite difference of matrix b using
% operator a. b has the dimensions x-y-t or x-y-z-t
function res = mtimes(a,b)

    if a.adjoint
        % Adjoint takes us from TV domain to image domain
        res = adjDt(b);
    else
        % This takes us from the image domain to the temporal TV domain
         if (ndims(b) == 4)
            res = b(:,:,:,[2:end,end]) - b;
         else
             res = b(:,:,[2:end,end]) - b;
         end
    end

end



function y = adjDt(x)
    if (ndims(x) == 4)
        y =  x(:,:,:,[1,1:end-1]) -x;
        y(:,:,:,1) = -x(:,:,:,1);
        y(:,:,:,end) = x(:,:,:,end-1);
    else
        y =  x(:,:,[1,1:end-1]) -x;
        y(:,:,1) = -x(:,:,1);
        y(:,:,end) = x(:,:,end-1);
    end

end


% Original adjDz function
% function y = adjDz(x)
% y= x(:,:,[1,1:end-1]) - x;
% y(:,:,1) = -x(:,:,1);
% y(:,:,end) = x(:,:,end-1);