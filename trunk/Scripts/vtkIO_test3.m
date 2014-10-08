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
    spMR = [sxUS syUS szUS];
    
    % Read the MR volume intensities
    fMRraw = dir([directory dPrefix dPatients(i,:) dMRSuffix '\*_DFD.raw']);
    fMRraw = fMRraw.name;
    fMRstr = fopen( [directory dPrefix dPatients(i,:) dMRSuffix '\' fMRraw], 'r');
    
    MR=fread(fMRstr,wUS*hUS*dUS,'uint8');
    MR=reshape(MR,wUS,hUS,dUS);
        
    % Set origin to be the center of the MR intensity volume
    origin = -0.5*[sxUS*wUS syUS*hUS szUS*dUS];
    
    order = [2 1 3];
    reMR = permute(MR, order);
    
    for j=1:size(reMR,3)
        reMR(:,:,j) = squeeze(flipdim(reMR(:,:,j),1));
    end
    
    re_sp = [spMR(order(1)) spMR(order(2)) spMR(order(3))];
    re_or = [origin(order(1)) origin(order(2)) origin(order(3))];
  
    writeImageDataInVTK(reMR, re_or, ['..\data\prostate\P' dPatients(i,:) '_MR_Reg_Volume.vtk'], re_sp );
    
    % Write out the US too
    fUSraw = dir([directory dPrefix dPatients(i,:) dMRSuffix '\*_US.raw']);
    fUSraw = fUSraw.name;
    fUSstr = fopen( [directory dPrefix dPatients(i,:) dMRSuffix '\' fUSraw], 'r');
    
    US=fread(fUSstr,wUS*hUS*dUS,'uint8');
    US=reshape(US,wUS,hUS,dUS);
    
    reUS = permute(US, order);
    
    for j=1:size(reUS,3)
        reUS(:,:,j) = squeeze(flipdim(reUS(:,:,j),1));
    end
  
    writeImageDataInVTK(reUS, re_or, ['..\data\prostate\P' dPatients(i,:) '_US_Reg_Volume.vtk'], re_sp );
end
