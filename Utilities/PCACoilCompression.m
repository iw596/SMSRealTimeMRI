% Perform SVD-coil compression for given kdata
% number of virtual coils is given by nvc
% Huang F, Vijayakumar S, Li Y, Hertel S, Duensing GR. A software channel compression technique
function compressionMatrix = PCACoilCompression(kdata,nvc,plotEnergy)
    if (nargin < 3)
     plotEnergy = false;
    end
    [U,S,V] = svd(kdata,"econ");
    compressionMatrix = V(:,1:nvc);

    if (plotEnergy == true)
        figure,plot((diag(S)).*100,'LineWidth',2);
        h = gca; % Get axis to modify.
        h.FontSize = 18; 
        xlabel("Number of Channels",'fontweight','bold','fontsize',25); 
        ylabel("Singular Values (%)",'fontweight','bold','fontsize',25);
    end
end


