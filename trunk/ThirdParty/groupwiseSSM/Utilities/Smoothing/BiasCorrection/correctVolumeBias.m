function B = correctVolumeBias(V, val)

D = single(squeeze(V));

if(nargin <= 2)
    [B,U]=BCFCM3D(D,[-300 70 250]',struct('maxit',5,'epsilon',1e-5,'sigma',1));
else
    [B,U]=BCFCM3D(D,val,struct('maxit',5,'epsilon',1e-5,'sigma',1));
end


