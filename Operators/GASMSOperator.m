% This class implements basic SENSE style reconstruction of 2D Simultaneous multislice data using
% coil sensitvity maps, through using anonymous functions to represent the
% forward and back FT's both Cartesian and Non-Cartesian trajectories can
% be used.

classdef GASMSOperator 
    properties
        FT = []; % FT is a vector of NUFFT operators
        dcf = [];
        maps = []; % Maps have the dimension x-y-slc-cha
        phaseCycling = []; % CAIPIRINHA phase cycling
        imgSize = []; % Imgs have the dimensions x-y-slc-t
        dataSize = []; % This has the dimensions kx-ky-t
        adjoint = 0;
        Nt = 1; % Number of time points
    end
    methods
        function obj = GASMSOperator(FT,dcf,maps,phaseCycling,imgSize,dataSize,Nt)
            obj.FT = FT;
            if  isempty(dcf)
                obj.dcf = 1;
            else
                obj.dcf = dcf;
            end
            obj.maps = maps;
            obj.phaseCycling = phaseCycling;
            obj.imgSize = imgSize;
            obj.dataSize = dataSize;
            obj.Nt = Nt;
        end

        function res = mtimes(obj,b)     
            if obj.adjoint == 0
                NSli = obj.imgSize(3);
                NCha = size(obj.maps,4);
                Nt = obj.Nt;
                % Maps b in the image domain to the kspace domain
                % Just to be safe we reshape our data it should look like
                % x-y-sli-t
                b = reshape(b,[obj.imgSize Nt]);
                res = zeros([obj.dataSize,size(obj.maps,4),Nt]);
                % Recontruct each frame indepdently (can this be a
                % parfor??)
                for t = 1:Nt
                    dataFrame = squeeze(b(:,:,:,t));
                    tmp = zeros([obj.dataSize(1:2) size(obj.maps,4)]);
                    for sli = 1:NSli
                        tmp = tmp + obj.phaseCycling(:,:,sli,t).*(obj.FT(t)*(squeeze(obj.maps(:,:,sli,:)).*dataFrame(:,:,sli)));
                    end
                    res(:,:,:,t) = tmp;
                end
                return;
            else
                % Maps b in the kspace domain to the image domain
                b = reshape(b,[obj.dataSize size(obj.maps,4) obj.Nt]);
                res = zeros([obj.imgSize obj.Nt]);
                NSli = obj.imgSize(3);
                for t = 1:obj.Nt
                    for slc = 1:NSli
                         tmp = sum(conj(squeeze(obj.maps(:,:,slc,:))).*(obj.FT(t)'*(conj(obj.phaseCycling(:,:,slc,t)) .* (b(:,:,:,t).* obj.dcf))),3);
                         res(:,:,slc,t) = tmp;
                    end
                end
                return;
            end
        end

        function res = ctranspose(obj)
            obj.adjoint = xor(obj.adjoint,1);
            res = obj;
        end
    
    end
end