function  res = TV_Temp()

%res = TV_Temp()
%
% Implements a difference operator along the time dimensin for dynamic MRI
%
% Original by Ricardo Otazo 2008 Ed
% Edited to work across multiple slices by IW 2023


res.adjoint = 0;
res = class(res,'TV_Temp');

