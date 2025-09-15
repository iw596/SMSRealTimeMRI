%%% Function to perform Alternating Direction Method of Multipliers Algorithm
%%% Require A: Forward Operator, y: K-space data T: Regularizer objects, imsize: image dimensions [x,y,z,Nframes]
function res = ADMM(A,y,T,lambda,NIters,NInnerIters,imsize,plotIterations,iterationPath)
     
    if nargin < 9
        saveIterations = false
    else
        saveIterations = true;
    end
    
    if nargin< 8
        plotIterations = true;
        saveIterations = false;
    end
    b = reshape(A'*(y),imsize);
    xOld = 0.*b;
    x = xOld;
    z = {};
    u = {};

    zOld = {};
    uOld = {};
    for r = 1:length(T)
        zOld{r} = T(r) * xOld;
        uOld{r} = zOld{r}
    end
  
    rho = ones(length(T),1);

    mu = 1.2; % ADMM residual rescaling tolerance
    tauMax = 10.0;
    AHA = @(x)(A'*((A*x)));
    for i = 1:NIters
        %% X-update
        lhsCalc = @(x)lhs(x,AHA,T,rho,imsize);
        rhs = zeros(size(vec(b),1),1);
        for r = 1:length(T)
            rhs = rhs + vec(b + rho(r) * (T(r)'*(zOld{r} - uOld{r})));
        end
        x = symmlq(lhsCalc,rhs, 1e-4, NInnerIters, [], [], xOld(:));
        x = reshape(x,imsize);
        for r = 1:length(T)
            normFx = 0;
            normz = 0;
            normu = 0;
            pRes = 0;
            dRes = 0;
            Fx = T(r) * x;
            Fxu = Fx + uOld{r};

            %% z-update
            z{r} = T(r).prox(Fxu,lambda(r)./rho(r));

            %% U-update
            u{r} = Fxu - z{r};

            %% Calculate some parameters for adaptive rho
            normFx = normFx + squaredNorm(Fx);
            normz = normz +  squaredNorm(z{r});
            normu = normu + squaredNorm(T(r)'*u{r});
            pRes = pRes + squaredNorm((Fx - z{r}));
            dRes = dRes + squaredNorm(T(r)'*(z{r} - zOld{r}));
            normFx = sqrt(normFx);
            normz = sqrt(normz);
            normu = sqrt(normu);
            pRes = sqrt(pRes)/max(normFx,normz);
            dRes = sqrt(dRes)./normu;
            ratio = sqrt(pRes./dRes);

            %% Adaptive rho
            if (ratio < 1)
                tau = max(1/tauMax,1/ratio);
            else
                tau = min(tauMax,ratio);
            end

            if (pRes > mu.*dRes)
                rho(r) = rho(r) * tau;
                u{r} =  u{r} ./ tau;

            elseif (dRes > mu.*pRes)
                rho(r) = rho(r) ./ tau;
                u{r} = u{r} .* tau;
            end
        end

        %% Update old u and old z
        for r = 1:length(T)
            uOld{r} = u{r};
            zOld{r}= z{r};
        end
        xOld = x;
        rho
        if (plotIterations == true)
            figure(44);
            if (ndims(x) ==3)
                imagesc(abs(x(:,:,10))); colormap(gray);axis off; axis image; brighten(0.3); drawnow;
            else
                imagesc(abs(x(:,:,3,10))); colormap(gray);axis off; axis image; brighten(0.3); drawnow;
            end
        end
        
        if (saveIterations == true)
            flName = sprintf("ADMM_Iteration_%d",i);
            flName = iterationPath + "\" + flName;
            save(flName, "xOld");
        end


    end
    res = x;
end


function res = lhs(x,AHA,T,rho,imgdims)
    x = reshape(x,imgdims);
    res = zeros(size(x));
    for r = 1:length(T)
        res = res + reshape(AHA(x),imgdims) + rho(r).*(T(r)'*(T(r)*x));
    end
    res = vec(res);

end