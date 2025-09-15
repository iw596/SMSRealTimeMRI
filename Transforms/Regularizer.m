
% This class provides a generic template for a regularizer and allows
% arbitary regularizers to be store together in a matrix
classdef Regularizer
    %UNTITLED12 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        T = []; % Regularizer that has been selected
        adjoint = 0;
        name = [];
        imsize = [];
        imSize_dyd = [];
    end

    methods
        function obj = Regularizer(name,imsize,winsize)
            obj.name = name;
            if (streq(name,"TV") == 1)
                obj.T = TV();
            elseif (streq(name,"TVT") == 1)
                obj.T = TVT();
            elseif (streq(name,"FD") == 1)
                obj.T = FD([1:2]); % Finite difference transform across 1st and 2nd dimensions
            elseif (streq(name,"Wavelet") == 1)
                if (nargin < 2)
                    error("No image size given")
                end
                obj.T = Wavelet(imsize,3,'db2');
            elseif (streq(name,"LLR") == 1)
                if (nargin < 3)
                    winsize = [8,8];
                    disp("No block size information passed in, assuming size of [8,8]")
                end
                obj.T = LLR(winsize);
            else
                error("No valid regularizer chosen...")
            end
            
        end

        % Applies transform or adjoint transform to data input b
        % we assume that the image b will always have the shape [nx ny nsli nt]
        function res = mtimes(obj,b) 
            
            if obj.adjoint == 0
                res = obj.T * b;
                return;
            else
                % Maps b in the transform domain to the image domain
                res = obj.T' *b;
                return;
            end


        end

        function res = ctranspose(obj)
            obj.adjoint = xor(obj.adjoint,1);
            res = obj;
        end
        
        % This method calculates the proximal operator of the regularizer
        function res = prox(obj,a,lambda)
            res = obj.T.prox(a,lambda);
        end

    end
end