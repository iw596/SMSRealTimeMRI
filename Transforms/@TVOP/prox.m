%% Implements the L1 proximal operator of finite difference stransform
%% this is soft-thresholding
function res = prox(obj,b,lambda)
    res= (abs(b)-lambda).*b./abs(b+eps).*(abs(b)>lambda); 
end

