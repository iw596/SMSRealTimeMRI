classdef FD
    %FD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dims = [];
        adjoint  = 0;

    end
    
    methods
        function obj = FD(dims)
            %FD Construct an instance of this class
            %   Detailed explanation goes here
            obj.dims = dims;
        end
        
        function Q = mtimes(A,B)
            % Performs the actual calculation of the finite difference.
            %   Called when using A*x or A'*x
            %   Input: A is usually the finite_difference_operator
            %          B is the image if A is not adjoint.
            %          B is the finite diff. of an image if A is adjoint
            %   Output: Finite difference or image itself, if A is not or
            %           is adjoint, respectively.

            if isa(A, 'FD')
                if A.adjoint==1
                    if isvector(B)
                        Q = B - circshift(B,-1);
                    else
                        s = size(B);
                        Q = zeros(s(2:end));
                        for id = 1:length(A.dims)
                            Q = Q + squeeze(B(id,:,:,:,:)) - circshift(squeeze(B(id,:,:,:,:)),-1,A.dims(id));
                        end
                    end
                else
                    if isvector(B)
                        Q = B - circshift(B,1);
                    else
                        Q = [];
                        for id = 1:length(A.dims)
                            Q = cat(1, Q, reshape(B - circshift(B,1,A.dims(id)), [1 size(B)]));
                        end
                    end
                end
            else   % now B is the operator and A is the vector
                Q = mtimes(B',A')';
            end
        end
        function res = ctranspose(obj)
            obj.adjoint = xor(obj.adjoint,1);
            res = obj;
        end

        function res = prox(obj,input,lambda)
            % Applies proximal operator to data b, proximal operator for TV
            % transform is soft thresholding. The thresholding level is
            % controlled by lambda
            abs_input = abs(input);

            sign = input ./ abs_input;
            sign(isnan(sign))=0;
            mag = abs_input - lambda;
            mag = (abs(mag) + mag) / 2;
            res = sign.*mag;

        end
    end
end

