% given the images and contours and desired filename (with directory), makes a stradwin compatible file

function write_stradwin(I ,dir, filename, x_spacing, y_spacing, z_spacing, origin)

%input contours: [3xm double] [3xn double]....(number of slices), m,n points in
%each slice . IN PIXELS
%x_scale,y_scale,spacing: in cm

% x_scale, y_scale: pixel size in cm
% spacing: distance between each slice in cm

sw_file=fopen(strcat(dir,filename,'.sw'),'wt');

fprintf(sw_file, 'RES_BUF_FRAMES %d\n', size(I,3));
fprintf(sw_file, 'RES_BUF_WIDTH %d\n', size(I,2));
fprintf(sw_file, 'RES_BUF_HEIGHT %d\n', size(I,1));
fprintf(sw_file, 'RES_POS_REC true\n');
fprintf(sw_file, 'RES_BUF_RF false\n');
fprintf(sw_file, 'RES_END_HEADER\n');
fprintf(sw_file, 'RES_XSCALE %.6f\n',x_spacing);
fprintf(sw_file, 'RES_YSCALE %.6f\n',y_spacing);
fprintf(sw_file, ['RES_BIN_IM_FILENAME ', filename,'.sxi\n']);

for i=1:size(I,3)
    fprintf(sw_file, ['IM 0 ' num2str(origin(1)) ' ' num2str(origin(2)) ' %6.5f 0 0 0 \n'], (i-1)*z_spacing + origin(3));
end
fclose(sw_file);

sxi_file = fopen(strcat(dir,filename,'.sxi'), 'w');

for i=1:size(I,3)
    for j=1:size(I,1)
        for k=1:size(I,2)
            fwrite(sxi_file, I(j,k,i));
        end
    end
end

fclose(sxi_file);