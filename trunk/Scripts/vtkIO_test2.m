add_bcpd_paths;

%% Folder, file and prefix containers
directory = '';
dPrefix = '';
dMRSuffix = '';
noPatients = 0;

%% Read the file list from the text file
flist = fopen('..\data\prostate\MR_TRUS.txt','r');
if (flist < 0)
    error(['Cannot open the file ''', flist,'''']);
end

%% Read the directory
if( strcmp( fgetl(flist), 'Directory') )
    directory = fgetl(flist);
    if (directory(end) ~='\')
        directory(end+1) ='\';
    end
else
    error(['Unable to read the directory where images are located''', flist,'''']);
end

fgetl(flist);

%% Read the prefix
if( strcmp( fgetl(flist), 'Patient directory prefix') )
    dPrefix = fgetl(flist);
else
    error(['Unable to read the directory prefix''', flist,'''']);
end

fgetl(flist);

%% Read the MR suffix
if( strcmp( fgetl(flist), 'Registered MR directory suffix ') )
    dMRSuffix = fgetl(flist);
else
    error(['Unable to read the suffix for registered MR images''', flist,'''']);
end

fgetl(flist);

%% Read patient folder numbers
if( strcmp( fgetl(flist), 'Folders') )
    noPatients = str2num(fgetl(flist));
else
    error(['Unable to read the number of patients''', flist,'''']);
end    

%% Read patient folders
dPatients = zeros(noPatients,3);

for i=1:noPatients
    dPatients(i,:) = fgetl(flist);
end

%% Create a VTK file for the MR and TRUS images in each directory
origin = [0 0 0]; % data lacks origin info
for i=1:noPatients
    
    % Read the TRUS volume info
    fUSinfo = dir([directory dPrefix dPatients(i,:) '\*3DQ.nfo']);
    fUSinfo = fUSinfo.name;
    fUSstr = fopen( [directory dPrefix dPatients(i,:) '\' fUSinfo], 'r');
    str = fgetl(fUSstr);
    wUS = sscanf(str, 'width: %u');
    str = fgetl(fUSstr);
    hUS = sscanf(str, 'height: %u');
    str = fgetl(fUSstr);
    dUS = sscanf(str, 'numframes: %u');
    str = fgetl(fUSstr);
    sxUS = sscanf(str, 'xvoxelsize: %f');
    str = fgetl(fUSstr);
    syUS = sscanf(str, 'yvoxelsize: %f');
    str = fgetl(fUSstr);
    szUS = sscanf(str, 'zvoxelsize: %f');
    
    % Read the TRUS volume intensities
    fUSraw = dir([directory dPrefix dPatients(i,:) '\*3DQ.raw']);
    fUSraw = fUSraw.name;
    fUSstr = fopen( [directory dPrefix dPatients(i,:) '\' fUSraw], 'r');
    
    US=fread(fUSstr,wUS*hUS*dUS,'uint8');
    US=reshape(US,wUS,hUS,dUS);
    
    % Set origin to be the center of the US intensity volume
    origin = -0.5*[sxUS*wUS syUS*hUS szUS*dUS];
    writeImageDataInVTK(US, origin, ['..\data\prostate\P' dPatients(i,:) '_US_Volume.vtk'], [sxUS syUS szUS] );
    
    % Read the MR volume info
    fMRinfo = dir([directory dPrefix dPatients(i,:) '\*8bits.nfo']);
    fMRinfo = fMRinfo.name;
    fMRstr = fopen( [directory dPrefix dPatients(i,:) '\' fMRinfo], 'r');
    str = fgetl(fMRstr);
    wMR = sscanf(str, 'width: %u');
    str = fgetl(fMRstr);
    hMR = sscanf(str, 'height: %u');
    str = fgetl(fMRstr);
    dMR = sscanf(str, 'numframes: %u');
    str = fgetl(fMRstr);
    sxMR = sscanf(str, 'xvoxelsize: %f');
    str = fgetl(fMRstr);
    syMR = sscanf(str, 'yvoxelsize: %f');
    str = fgetl(fMRstr);
    szMR = sscanf(str, 'zvoxelsize: %f');
    spMR = [sxMR syMR szMR];
    
    % Read the MR volume intensities
    fMRraw = dir([directory dPrefix dPatients(i,:) '\*8bits.raw']);
    fMRraw = fMRraw.name;
    fMRstr = fopen( [directory dPrefix dPatients(i,:) '\' fMRraw], 'r');
    
    MR=fread(fMRstr,wMR*hMR*dMR,'int8');
    MR=reshape(MR,wMR,hMR,dMR);
    
    % Set origin to be the center of the MR intensity volume
    origin = -0.5*[sxMR*wMR syMR*hMR szMR*dMR];
    
    % Resample the MR such that the x-axis is normal to the axial plane
    %reMR = zeros(size(MR,2), size(MR,1), size(MR,3)); %1 2 3; 3 2 1; 2 1 3
    order = [2 1 3];
    re_sp_MR = [spMR(order(1)) spMR(order(2)) spMR(order(3))];
    re_or_MR = [origin(order(1)) origin(order(2)) origin(order(3))];
    
    reMR = permute(MR, order);
    
    for j=1:size(reMR,3)
        reMR(:,:,j) = squeeze(flipdim(reMR(:,:,j),1));
    end
    
    writeImageDataInVTK(reMR, re_or_MR, ['..\data\prostate\P' dPatients(i,:) '_MR_Volume.vtk'], re_sp_MR );
    
    % Read the segmented MR and resample it in US coordinates
    R = rotation_matrix(0, 0, -pi/2, 'radians');
    fMRinfo = dir([directory dPrefix dPatients(i,:) '\*_MR_Prostate.vtk']);
    fMRinfo = fMRinfo.name;
    [v,f] = readPolyDataInVTK( [directory dPrefix dPatients(i,:) '\' fMRinfo] );
    v = v*R;
    writePolyDataInVTK(v,f, ['..\data\prostate\P' dPatients(i,:) '_MR_Seg.vtk']);
    
    % Copy the segmented TRUS and copy it to the data folder
    fUSinfo = dir([directory dPrefix dPatients(i,:) '\*_US_Prostate.vtk']);
    fUSinfo = fUSinfo.name;
    copyfile([directory dPrefix dPatients(i,:) '\' fUSinfo],['..\data\prostate\P' dPatients(i,:) '_US_Seg.vtk']);
    
end