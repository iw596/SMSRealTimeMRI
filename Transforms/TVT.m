classdef TVT
    % T This class implements a spatial finite difference transfom in the temporal direction
    % The proximal operator of this transform is the soft-thresholding function 
    
    
    properties
        adjoint  = 0;
    end
    
    methods
        function obj = TVT()

        end
        
        function res = mtimes(obj,b)
            if obj.adjoint
                % Adjoint takes us from the temporal TV domain to image domain
                res = obj.adjDt(b);
            else
                % This takes us from the image domain to the temporal TV domain
                 if (ndims(b) == 4)
                    res = b(:,:,:,[2:end,end]) - b;
                 else
                     res = b(:,:,[2:end,end]) - b;
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
            res = (abs(b)-lambda).*b./abs(b+eps).*(abs(b)>lambda);
        end
    end

    methods(Static)
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
    end
end

