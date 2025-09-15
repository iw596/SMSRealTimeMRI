function [traj,angles] = GenerateRadialTrajectory(params)
NPoints = params.NCol;
NSpokes = params.NLin;
k = zeros(NPoints,NSpokes);
angles = zeros(NSpokes,1);
%rho = linspace(-0.5,0.5 - 1/NPoints,NPoints);
rho = linspace(-0.5,0.5,NPoints);
%rho = linspace(-NPoints/2,NPoints/2-1,NPoints);
%if (~isempty(params.radialMode))
 %   radialMode = params.radialMode;
%end
radialMode = cell2mat(params.WIPParameters(7));
%mode = params.mode;
if (isempty(radialMode))
    if (mod(NSpokes,2) == 0)
       turn = pi/NSpokes;
       for i = 0:NSpokes-1
            angles(i+1) = turn * i;
       end
    else
       turn = 2*pi/NSpokes;
       for i = 0:NSpokes-1
            angles(i+1) = turn * i;
       end
    end
elseif(radialMode == 1)
    disp("Generating Segmented trajectory")
    NSegments = (params.WIPParameters{10});
    NShots = params.NLin/NSegments;
    segmentOffset = 2*pi/NSpokes;
    shotAngle = 2*pi/NShots;
    count = 1;
    for i = 0:NSegments-1
        for j = 0:NShots-1
            angles(count) = j*shotAngle + i*segmentOffset; 
            count = count + 1;
        end
    end
elseif (radialMode == 2)
    % GR trajectory
        disp("Generating standard non continous GA trajectory")

    GR = (1.0 + sqrt(5.0)) / 2.0;
    count = 1;
    for (i = 0:NSpokes-1)
        angle = i*180.0 / GR;
        angle = mod(angle, 360.0);
        %Convert degrees to radians
        angle = angle * (pi / 180.0);
        angles(i+1) = angle;
    end
elseif (radialMode == 4)
    % Continous GA trajectory, this is more involved to calculate...
        disp("Generating continous GA trajectory")

    GR = (1.0 + sqrt(5.0)) / 2.0;
    count = 1;
    k = repmat(k,params.NRep);
    NSpokes = NSpokes * params.NRep;
    for (i = 0:NSpokes-1)
        angle = i * 180.0/GR;
        angle = mod(angle,360.0);
        %Convert degrees to radians
        angle = angle * (pi / 180.0);
        angles(i+1) = angle;

    end

elseif (radialMode == 5)
    disp("Generating adapted GA trajectory for SMS acquisition")
    NSli = params.RefNSli;
    GR = (1.0 + sqrt(5.0))./2.0;
    count = 1;
    % Recalculate NSpokes and size of angles array
    NSpokes = NSpokes*params.NRep;
    angles = zeros(NSpokes,1);
    for i = 0:NSpokes-1
        angle = (i*(180.0./GR))./NSli;
        angle = mod(angle,360.0);
        %Convert degrees to radians
        angle = angle * (pi / 180.0);
        angles(count) = angle;
        count = count + 1;
    end
end
traj = rho'*exp(1j*angles)';
    
end