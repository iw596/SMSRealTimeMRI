classdef GASENSEOperator
    %UNTITLED14 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        NUFFT   = [];% NUFFT is a vector of NUFFT operators for each frame
        dcf = []; % Density compensation function
        maps = [];
        imgSize = [];
        dataSize = [];
        adjoint = 0;
        Nt = 1; % Number of time points
    end

    methods
        function obj = GASENSEOperator(NUFFT,dcf,maps,imgSize,dataSize,Nt)
            obj.NUFFT = NUFFT;
            obj.dcf = dcf;
            if isempty(maps)
                maps   =   ones([imgSize 1]);
            end

            if isempty(dcf)
                obj.dcf   =   ones([dataSize size(maps,3) Nt]);
            else
                obj.dcf = dcf;
            end
            obj.maps = maps;
            obj.imgSize = imgSize;
            obj.dataSize = dataSize;
            obj.Nt = Nt;
        end

        function res = mtimes(obj,b)     
            if obj.adjoint == 0
                % Maps b in the image domain to the kspace domain
                % Just to be safe we reshape our data
                b = reshape(b,[obj.imgSize obj.Nt]);
                res = zeros([obj.dataSize size(obj.maps,3) obj.Nt]);
                for t = 1:obj.Nt
                    res(:,:,:,t) = obj.NUFFT(t)* (obj.maps.*b(:,:,t));
                    %res(:,:,:,t) = res(:,:,:,t).*obj.dcf(:,:,:,t);
                end
            else
                % Maps b in the kspace domain to the image domain
                b = reshape(b,[obj.dataSize size(obj.maps,3) obj.Nt]);
                res = zeros([obj.imgSize obj.Nt]);
                for t = 1:obj.Nt
                    res(:,:,t) = sum(conj(obj.maps).*(obj.NUFFT(t)'*(b(:,:,:,t))),3);
                end
            end
        end

        function res = ctranspose(obj)
            obj.adjoint = xor(obj.adjoint,1);
            res = obj;
        end
        
        % mtimes 2 implements AHA*b
        function res = mtimes2(obj,b)
            res = obj'*(obj*b);
        end
        function est = iter(obj, d, optfn, tol, iters, L)

            %   Performs symmetric iterative recon using built-ins
            %   Input d should be shaped like the output of mtimes
            %   Solves normal equation

            if nargin < 3
                optfn   =   @minres;
            end
            if nargin < 4
                tol     =   [];
            end
            if nargin < 5
                iters   =   100;
            end


            if isequal(size(d), [obj.dataSize]) || (numel(d) == prod(obj.dataSize))
                d = obj'*d;
            end

            [est, ~, relres, iter] =   optfn(@(x,mode) reshape(mtimes2(obj, reshape(x,[obj.imgSize obj.Nt])),[],1), reshape(d,[],1), tol, iters, [], []);

           % est =   reshape(est, obj.msize);
            %est =   reshape(est, obj.imgSize);
            fprintf(1, 'Exit after %i iterations, residual: %G\n', iter, relres);
        end
        end

        
        methods (Static)
            function b = fftfn(b,dims)
                for i = dims
                    b   =   fftshift(fft(ifftshift(b, i), [], i), i)/sqrt(size(b,i));
                end
            end
        end
end