add_bcpd_paths;

%% Read in the raw file
im_size = [352 373 644];
im_spacing = [0.168015 0.156006 0.156006];
fp=fopen('C:\data\CPD\UWO\1.raw','r');
I1=fread(fp,im_size(1)*im_size(2)*im_size(3),'uint8');
I=reshape(I1,im_size(1),im_size(2),im_size(3));

origin = [0 0 0];

%% Re-arrange the raw image such that axial slices are along the Z-axis
J = zeros(size(I,2),size(I,3),size(I,1));
for i=1:size(I,1)
    J(:,:,i) = squeeze(I(i,:,:));
end

%% Write out the stradwin file for segmentation
scale = 0.1; step=10;
write_stradwin(J(:,:,1:step:end) ,'C:\data\CPD\UWO\','1', im_spacing(2)*scale, im_spacing(3)*scale, im_spacing(1)*step*scale, origin);