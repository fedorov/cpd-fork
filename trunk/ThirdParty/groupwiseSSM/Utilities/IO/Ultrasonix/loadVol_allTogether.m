function [volume, datatype, numVols, fpV, h, w, ss, degPerFrm]  = loadVol_allTogether(datapath)
fid= fopen(datapath, 'r');
volinfo = fread(fid, 7, 'int');

datatype = volinfo(1);  
numVols = volinfo(2);
fpV = volinfo(3);
w = volinfo(4);
h = volinfo(5);
ss = volinfo(6);
degPerFrm = volinfo(7);
volData = fread(fid, inf, 'uchar=>uchar');
volume = reshape(volData, w, h, fpV*numVols);