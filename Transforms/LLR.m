% Class to implement Locally-low rank sparsifying transform, a lot of this
% code was taken from T2 shuffling recon github page: https://github.com/utcsilab/t2shuffling-support/blob/master/src/utils/llr_thresh.m
classdef LLR
    % This forms a set of blocks for an image of dimensions nx X ny X nt
    %winSize: block size for LLR (two integers). 
    properties
        winsize = [];
        adjoint = 0;
        imgSize = [];
        zpad = []; % image size required afer zeropadding
        bx;
        by;
        zpadx;
        zpady;
    end
    
    methods
        % Contructor method 
        function obj = LLR(winsize)
            %LLR Construct an instance of this class
            %   Detailed explanation goes here
            obj.winsize = winsize;
        end
        
        function res = mtimes(obj,b)     
            % MTimes for LLR is simply multiplication with Identity matrix,
            % so just pass the data through
            if obj.adjoint == 0
               % Maps b in image domain to image patches
               res = b;
               return;
            else
                res = b;
                return;
            end
        end

        % Proximal operator for LLR breaks images into patches and uses SVD
        % to decompose each patch , b is the input data of dimensions
        % nx x ny x t or nx x ny x sli x t
        function res = prox(obj,b,lambda)
            res = zeros(size(b));
            if (ndims(b) == 4)
                % Apply prox to each slice independently
                for sli =1:size(b,3)
                    res(:,:,sli,:) = obj.patches(squeeze(b(:,:,sli,:)),lambda);
                end
            else
                res = obj.patches(b,lambda);
            end
        end
       
        function res = ctranspose(obj)
            obj.adjoint = xor(obj.adjoint,1);
            res = obj;
        end
    end
    methods(Hidden)
        function res = patches(obj,img,lambda)
            [nx, ny, nt] = size(img); % In original code nt was K
            %spatial dimensions of image block
            Wx = obj.winsize(1);
            Wy = obj.winsize(2);
            L = nx * ny / Wx / Wy; % Number of blocks
            % Random shifting to reduce block artefacts
            shift_idx = [randi(Wx), randi(Wy) 0];
            img = circshift(img, shift_idx);

            img_LLR = zeros(Wx*Wy, L, nt);
            for ii=1:nt
                img_LLR(:,:,ii) = im2col(img(:,:,ii), [Wx,Wy], 'distinct');
            end
            img_LLR = permute(img_LLR, [1 3 2]);
            s_LLR = zeros(nt, L);
            % threshold singular values
            for ii=1:L
                [U, S, V] = svd(img_LLR(:,:,ii), 'econ');
                sval = diag(S);
                s2 = obj.SoftThresh(sval, lambda);
                img_LLR(:,:,ii) = U*diag(s2)*V';
            end
            % Reshape thresholded patches back into image
            img_thresh = zeros(size(img));
            for ii=1:nt
                img_thresh(:,:,ii) = col2im(img_LLR(:,ii,:), [Wx, Wy], [nx, ny], 'distinct');
            end
            img_thresh = circshift(img_thresh, -1 * shift_idx);
            res = img_thresh;

        end
    end

    methods(Static)
        function res = SoftThresh(x, lambda)
            % Soft thresholding of x 
            res = x.*(abs(x) - lambda)./abs(x);
            res(abs(x) <= lambda) = 0;
        end

    end
end

