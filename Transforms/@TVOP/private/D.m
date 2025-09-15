function res = D(image)

%
% res = D(image)
%
% image = a 2D image
%
% This function computes the finite difference transform of the image
%
% Related functions:
%       adjD , invD 
%
%
% (c) Michael Lustig 2005

[sx,sy,nt] = size(image);
res = zeros(sx,sy,2,nt);
for t = 1:nt
    Dx = image([2:end,end],:,t) - image(:,:,t);
    Dy = image(:,[2:end,end],t) - image(:,:,t);
    res(:,:,1,t) = Dx;
    res(:,:,2,t) = Dy;
end

%res = [sum(image(:))/sqrt(sx*sy); Dx(:);  Dy(:)]; 
% if (nt == 1)
% res = cat(3,Dx,Dy);
% else
%     res = cat(4,Dx,Dy);
% end

