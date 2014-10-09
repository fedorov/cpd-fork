function writeFEMvtk(fem,u,fileName,apply)
  fid = fopen(fileName,'w');
  % header
  fprintf(fid,'# vtk DataFile Version 2.0\nFEM result\nASCII\nDATASET UNSTRUCTURED_GRID\n');
  
  % points
  numPoints = fem.getNumNodes;
  fprintf(fid,'POINTS %i float\n',numPoints);
  for i=1:numPoints
      point = fem.getNode(i);
      if apply == false
        fprintf(fid,'%f %f %f\n',point(1),point(2),point(3));      
      else
        fprintf(fid,'%f %f %f\n',point(1)+u(i,1),point(2)+u(i,2),point(3)+u(i,3));
      end
  end
  
  fprintf(fid,'\n');
  
  % cells
  numCells = fem.getNumElements;
  fprintf(fid,'CELLS %i %i\n', numCells, numCells*5);
  for i=1:numCells
      cell = fem.getElement(i);
      indices = cell.getNodeIdxs;
      fprintf(fid,'4 %i %i %i %i\n',indices(1)-1, indices(2)-1, indices(3)-1, indices(4)-1);
  end
  
  fprintf(fid,'\n');
  
  % cell types
  fprintf(fid,'CELL_TYPES %i\n',numCells);
  for i=1:numCells
      fprintf(fid,'10\n');
  end
  
  fprintf(fid,'\n');

  % point data
  fprintf(fid,'POINT_DATA %i\nVECTORS displacements float\n',numPoints);
  for i=1:numPoints
      if apply == false
        fprintf(fid,'%f %f %f\n',u(i,1),u(i,2),u(i,3));
      else
        fprintf(fid,'%f %f %f\n',-u(i,1),-u(i,2),-u(i,3));
      end
  end

  fprintf(fid,'\n');
    
  % cell data
  %fprintf(fid,'CELL_DATA %i\nFIELD FieldData 1\nLabel 1 %i int\n',numCells,numCells);
  %for i=1:numCells
  %    fprintf(fid,'0\n');
  %end
  
  fclose(fid);
end