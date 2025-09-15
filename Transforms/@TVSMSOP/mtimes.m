% mtime applys either the adjoint or forward TV SMS transform to the data
% a is the TV object and b is the data in the form x-y-sli-t

function res = mtimes(a,b)
if a.adjoint
    Nt = size(b,5); % Number of time points
    NSli = size(b,4);
    nx = size(b,1);
    ny = size(b,2);
    res = zeros(nx,ny,NSli,Nt);
    % Maps transform domain to image space
    for t = 1:Nt
        for sli = 1:NSli
            res(:,:,sli,t) = adjD(squeeze(b(:,:,:,sli,t)));
        end
    end
else
    % Maps image space to transform domain
    imgDims = zeros([size(b)]);
    res = zeros(size(imgDims,1),size(imgDims,2),2,size(imgDims,3),size(imgDims,4));
    for t = 1:size(imgDims,4)
        for i = 1:size(imgDims,3)
            res(:,:,:,i,t) = D(squeeze(b(:,:,i,t)));
        end
    end
end

end




    
