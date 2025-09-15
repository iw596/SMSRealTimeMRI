classdef GPUNUFFT
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        k = []; % k-space trajectory
        imsize = [];
        datasize = [];
        NCha = 1;
        params = [];
        adjoint = 0;
    end

    methods
        function obj = GPUNUFFT(k,imsize,datasize,NCha,osr,kernelwidth)
            % Constructor method takes trajectory, imagesize, oversampling
            % ratio and kernel width
            traj(1,:) = real(k(:));
            traj(2,:) = imag(k(:));
            obj.k = single(traj)';
            obj.imsize = uint32(imsize);
            obj.datasize = datasize;
            obj.NCha = NCha;
            obj.params.kernelWidth = (kernelwidth);
            obj.params.osr = (osr);
        end
        % Data is of the form NCol x NLin x NCha or 
        % img width x im height x NCha
        function res = mtimes(obj,data)  
             if obj.adjoint == 1
                 % Adjoint takes k-space to image domain
                 % First vectorise input to get in form of NPoints x NCha
                 data = reshape(data,[prod(obj.datasize) obj.NCha]);
                 tmp = zeros([2,prod(obj.datasize),obj.NCha]);
                 tmp(1,:,:) = real(data);
                 tmp(2,:,:) = imag(data);
                 outdata = MexAdjoint2D(single(tmp),obj.k,obj.imsize,obj.params);
                 outdata = complex(outdata(1,:),outdata(2,:));
                 res = reshape(outdata,[obj.imsize,obj.NCha]);
             else
                 % Image domain to k-space
                 data = reshape(data,[prod(obj.imsize) obj.NCha]);
                 tmp = zeros([2,prod(obj.imsize),obj.NCha]);
                 tmp(1,:,:) = real(data);
                 tmp(2,:,:) = imag(data);
                 outdata = MexForward2D(single(tmp),obj.k,obj.imsize,obj.params);
                 outdata = complex(outdata(1,:),outdata(2,:));
                 res = reshape(outdata,[obj.datasize,obj.NCha]);

             end
        end

        function res = ctranspose(obj)
            obj.adjoint = xor(obj.adjoint,1);
            res = obj;
        end
    end
end