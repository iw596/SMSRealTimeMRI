% This class implements basic SENSE style reconstruction of 2D data using
% coil sensitvity maps, through using anonymous functions to represent the
% forward and back FT's both Cartesian and Non-Cartesian trajectories can
% be used

classdef SENSEOperator


    properties
        FT   = [];% FT is an anonymous function which maps image space to kspace
        FTH = []; % FTH is anonymous function which maps kspace to the image domain;
        maps = [];
        imgSize = [];
        dataSize = [];
        mask = [];
        adjoint = 0;
        Nt = 1; % Number of time points
    end

    methods
        function obj = SENSEOperator(FT,FTH,maps,mask,imgSize,dataSize,Nt)
            obj.FT = FT;
            obj.FTH = FTH;
            obj.maps = maps;
            obj.imgSize = imgSize;
            obj.dataSize = dataSize;
            obj.Nt = Nt;
            
            if (isempty(mask))
                obj.mask = ones(size(dataSize));
            else
                obj.mask = mask;
            end
        end

        function res = mtimes(obj,b)     
            if obj.adjoint == 0
                % Maps b in the image domain to the kspace domain
                % Just to be safe we reshape our data
                b = reshape(b,[obj.imgSize obj.Nt]);
                res = zeros([obj.dataSize size(obj.maps,3) obj.Nt]);
                for t = 1:obj.Nt
                    res(:,:,:,t) = obj.mask.*obj.FT(obj.maps.*b(:,:,t));
                end
                %for ch = 1:size(obj.maps,3)
                   %res(:,:,ch) = obj.FT(obj.maps(:,:,ch).*b);
                %end
                return;
            else
                % Maps b in the kspace domain to the image domain
                b = reshape(b,[obj.dataSize size(obj.maps,3) obj.Nt]);
                res = zeros([obj.imgSize obj.Nt]);
                for t = 1:obj.Nt
                    res(:,:,t) = sum(conj(obj.maps).*obj.FTH(obj.mask.*b(:,:,:,t)),3);
                end
                % for ch = 1:size(obj.maps,3)
                %     res =  res + conj(obj.maps(:,:,ch)).*obj.FTH(squeeze(b(:,:,ch)));
                % end
                return;
            end
        end

        function res = ctranspose(obj)
            obj.adjoint = xor(obj.adjoint,1);
            res = obj;
        end
    
    end
end