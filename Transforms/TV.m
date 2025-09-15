classdef TV
    % T This class implements a spatial finite difference transfom in the x
    % and y directions. The proximal operator of this transform is the soft-thresholding function 
    
    properties
        adjoint  = 0;
    end
    
    methods
        function obj = TV()
        end
        
        % Implements the forward and adjoint transformas
        function res  = mtimes(obj,b)
            if (obj.adjoint ==0)
                % Maps image b to total variation domain
                if (ismatrix(b))
                    res = obj.D(b);

                elseif (ndims(b) == 3)
                    res = zeros([size(b,[1 2]),2 size(b,3)]);
                    % Apply TV to each frame individually
                    for t = 1:size(b,3)
                        res(:,:,1,t,:) = obj.D(b(:,:,t));
                    end
                else
                     res = zeros([size(b),2]);
                     % Apply TV to each frame and slice individually
                     for t = 1:size(b,4)
                         for sli = 1:size(b,3)
                             res(:,:,sli,t,:) = obj.D(b(:,:,sli,t));
                         end
                     end
                
                end

            else
                % Maps total variation domain b to image domain
                if (ndims(b) == 3)
                    res = obj.adjD(b);
                    return;

                elseif (ndims(b) == 4)
                    res = zeros([size(b,[1 2 4])]);
                    % Apply Adjoint TV to each frame individually
                    for t = 1:size(res,3)
                        res(:,:,t) = obj.adjD(b(:,:,:,t));
                    end
                    return;

                else
                     %res = zeros([size(b,[1 2 4 5])]);
                     res = obj.adjD(b);
                     % Apply Adjoint TV to each frame and slice individually
                     % for t = 1:size(res,4)
                     %     for sli = 1:size(res,3)
                     %        res(:,:,sli,t) = obj.adjD(b(:,:,:,sli,t));
                     %     end
                     % end
                     return;
                end

            end
        end
        
        function res = ctranspose(obj)
            obj.adjoint = xor(obj.adjoint,1);
            res = obj;
        end

        function res = prox(obj,b,lambda)
            % Applies proximal operator to data b, proximal operator for TV
            % transform is soft thresholding. The thresholding level is
            % controlled by lambda
            s   =   sqrt(sum(abs(b).^2,5));
            res   =   (max(abs(s)-lambda,0)./(s+eps)).*b;
        end

    end
    
    methods(Hidden)
        % Applies TV transform to 2D image
        function res = D(obj,image)
            [sx,sy,~] = size(image);
            res = zeros(sx,sy,2);
            Dx = image([2:end,end],:) - image;
            Dy = image(:,[2:end,end]) - image;
            res(:,:,1) = Dx;
            res(:,:,2) = Dy;
        end

        function res = adjD(obj,y)
            res = obj.adjDx(squeeze(y(:,:,:,:,1))) + obj.adjDy(squeeze(y(:,:,:,:,2)));
            return;
        end
    end
    methods(Static)
        function res = adjDy(x)
            res = x(:,[1,1:end-1],:,:) - x;
            res(:,1,:,:) = -x(:,1,:,:);
            res(:,end,:,:) = x(:,end-1,:,:);
        end

        function res = adjDx(x)
            res = x([1,1:end-1],:,:,:) - x;
            res(1,:,:,:) = -x(1,:,:,:);
            res(end,:,:,:) = x(end-1,:,:,:);
        end
    end
        
end




