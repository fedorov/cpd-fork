function write_fcsv( Points, Fiducial_Prefix, FileName )

fid = fopen(FileName, 'w');
fprintf(fid, '# Markups fiducial file version = 4.3\n');
fprintf(fid, '# CoordinateSystem = 0\n');
fprintf(fid, '# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n');

for i=1:size(Points,1)
    fprintf(fid, ['vtkMRMLMarkupsFiducialNode_' num2str(i) ',' num2str(Points(i,1)) ',' num2str(Points(i,2)) ',' num2str(Points(i,3)) ',0,0,0,1,1,1,0,' Fiducial_Prefix '-' num2str(i) ',,vtkMRMLScalarVolumeNode2\n' ]);
end

fclose(fid);

end

