% This class implements basic SENSE style reconstruction of 2D Simultaneous multislice data using
% coil sensitvity maps, through using anonymous functions to represent the
% forward and back FT's both Cartesian and Non-Cartesian trajectories can
% be used.

classdef SMSOperator 
    properties
        FT   = [];% FT is an anonymous function which maps image space to kspace
        FTH = []; % FTH is anonymous function which maps kspace to the image domain;
        maps = []; % Maps have the dimension x-y-slc-cha
        phaseCycling = []; % CAIPIRINHA phase cycling
        imgSize = []; % Imgs have the dimensions x-y-slc-t
        dataSize = []; % This has the dimensions kx-ky-t
        adjoint = 0;
        Nt = 1; % Number of time points
    end
    methods
        function obj = SMSOperator(FT,FTH,maps,phaseCycling,imgSize,dataSize,Nt)
            obj.FT = FT;
            obj.FTH = FTH;
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
                        tmp = tmp + obj.phaseCycling(:,:,sli).*obj.FT(squeeze(obj.maps(:,:,sli,:)).*dataFrame(:,:,sli));
                    end
                    res(:,:,:,t) = tmp;
                    % for ch = 1:NCha
                    %     tmp = zeros(obj.dataSize(1:2));
                    %     for sli = 1:NSli
                    %         tmp = tmp + obj.phaseCycling(:,:,sli).*obj.FT(obj.maps(:,:,sli,ch).*dataFrame(:,:,sli));
                    %     end
                    %     res(:,:,ch,t) = tmp;
                    % end
                end
                return;
            else
                % Maps b in the kspace domain to the image domain
                Nt = obj.Nt;
                b = reshape(b,[obj.dataSize size(obj.maps,4) Nt]);
                res = zeros([obj.imgSize Nt]);
                NSli = obj.imgSize(3);
                for t = 1:Nt
                    for slc = 1:NSli
                         %tmp = zeros([obj.imgSize(1:2)]
                         tmp = sum(conj(squeeze(obj.maps(:,:,slc,:))).*obj.FTH(conj(obj.phaseCycling(:,:,slc)) .* b(:,:,:,t)),3);
                         res(:,:,slc,t) = tmp;
                    end

                    % for slc = 1:NSli
                    %     tmp = zeros([obj.imgSize(1:2)]);
                    %     for ch = 1:size(obj.maps,4)
                    %         tmp =  tmp + conj(obj.maps(:,:,slc,ch)).*obj.FTH(conj(obj.phaseCycling(:,:,slc)) .* b(:,:,ch,t));
                    %     end
                    %     res(:,:,slc,t) = tmp;
                    % end
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