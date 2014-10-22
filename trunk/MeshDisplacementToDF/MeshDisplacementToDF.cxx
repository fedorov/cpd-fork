#include <itkImage.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include "itkImageToVTKImageFilter.h"
#include "itkVTKImageToImageFilter.h"
#include "itkImageDuplicator.h"

#include <iostream>

#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkDataSet.h"
#include "vtkFieldData.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkTetra.h"
#include "vtkProbeFilter.h"
#include "vtkCellArray.h"

#include <MeshDisplacementToDFCLP.h>

typedef int InputPixelType;
typedef itk::Vector<float,3> OutputPixelType;
typedef itk::Image<InputPixelType,3> InputImageType;
typedef itk::Image<OutputPixelType,3>  OutputImageType;
typedef itk::ImageFileReader<InputImageType> ReaderType;
typedef itk::ImageFileWriter<OutputImageType> WriterType;
typedef itk::ImageFileWriter<InputImageType> MaskWriterType;

typedef itk::ImageToVTKImageFilter<InputImageType> ToVTKConvertType;
typedef itk::VTKImageToImageFilter<InputImageType> FromVTKConvertType;

typedef itk::ImageDuplicator<InputImageType> DuplicateType;

int main(int argc, char* argv[]){
  PARSE_ARGS;

  // read input image
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(refImageName.c_str());
  reader->Update();
  InputImageType::Pointer ref = reader->GetOutput();
  std::cout << ref << std::endl;

  // read the input mesh
  vtkSmartPointer<vtkUnstructuredGridReader> meshReader = 
    vtkSmartPointer<vtkUnstructuredGridReader>::New();
  meshReader->SetFileName(meshName.c_str());
  meshReader->Update();
  vtkUnstructuredGrid *mesh = meshReader->GetOutput();
  vtkDataArray *meshDeformation = mesh->GetPointData()->GetArray("displacements");
  if(meshDeformation == NULL){
    std::cerr << "Failed to load deformation" << std::endl;
    return -1;
  }

  // allocate output image
  OutputImageType::Pointer df = OutputImageType::New();
  df->SetOrigin(ref->GetOrigin());
  df->SetSpacing(ref->GetSpacing());
  df->SetRegions(ref->GetBufferedRegion());
  df->SetDirection(ref->GetDirection());
  df->Allocate();

  // allocate output inverse image
  OutputImageType::Pointer dfInverse = OutputImageType::New();
  dfInverse->SetOrigin(ref->GetOrigin());
  dfInverse->SetSpacing(ref->GetSpacing());
  dfInverse->SetRegions(ref->GetBufferedRegion());
  dfInverse->SetDirection(ref->GetDirection());
  dfInverse->Allocate();

  DuplicateType::Pointer dup = DuplicateType::New();
  dup->SetInputImage(ref);
  dup->Update();
  InputImageType::Pointer mask = dup->GetOutput();
  mask->FillBuffer(0);

#if USE_vtkProbeFilter
  // prepare the point set
  ToVTKConvertType *toVtk = ToVTKConvertType::New();
  toVtk->SetInput(reader->GetOutput());  
  toVtk->Update();

  vtkSmartPointer<vtkProbeFilter> pf = vtkSmartPointer<vtkProbeFilter>::New();
  pf->SetSourceConnection(meshReader->GetOutputPort());
  pf->SetInputData(toVtk->GetOutput());
  pf->Update();
  vtkImageData *probed = pf->GetImageDataOutput();

  // iterate over the probed volume and transfer the displacements to ITK image
  vtkDataArray* imageDF = probed->GetPointData()->GetArray("Deformation");
  for(int i=0;i<probed->GetNumberOfPoints();i++){
    double* pt = probed->GetPoint(i);
    double* df = imageDF->GetTuple3(i);

    std::cout << pt[0] << "," << pt[1] << "," << pt[2] << std::endl;
    std::cout << df[0] << "," << df[1] << "," << df[2] << std::endl;

  }
#else
  vtkPoints *meshPoints = mesh->GetPoints();
  vtkCellArray *cells = mesh->GetCells();
  vtkIdType npts, *pts;
  int i, cellId=0;
  
  for(cells->InitTraversal();cells->GetNextCell(npts,pts);cellId++){
    assert(npts==4);
    vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();
    vtkCell* thisCell = mesh->GetCell(cellId);
    tetra->Initialize(npts, pts, meshPoints);
    vtkPoints *tetraPoints = tetra->GetPoints();

    InputImageType::IndexType bbMin, bbMax;
    bbMin[0] = ref->GetBufferedRegion().GetSize()[0];
    bbMin[1] = ref->GetBufferedRegion().GetSize()[1];
    bbMin[2] = ref->GetBufferedRegion().GetSize()[2];
    bbMax = ref->GetBufferedRegion().GetIndex();

    for(i=0;i<4;i++){
      InputImageType::PointType point;
      InputImageType::IndexType index;
      point[0] = -1.*tetraPoints->GetPoint(i)[0];
      point[1] = -1.*tetraPoints->GetPoint(i)[1];
      point[2] = tetraPoints->GetPoint(i)[2];
      bool isInside = ref->TransformPhysicalPointToIndex(point, index);
      if(!isInside){
        std::cerr << "Tetra point maps outside the reference image!" << std::endl;
        abort();
      }

      if(index[0]<bbMin[0])
        bbMin[0] = index[0];
      if(index[1]<bbMin[1])
        bbMin[1] = index[1];
      if(index[2]<bbMin[2])
        bbMin[2] = index[2];
      
      if(index[0]>bbMax[0])
        bbMax[0] = index[0];
      if(index[1]>bbMax[1])
        bbMax[1] = index[1];
      if(index[2]>bbMax[2])
        bbMax[2] = index[2];
    }

    bbMax[0]++;
    bbMax[1]++;
    bbMax[2]++;
    bbMin[0]--;
    bbMin[1]--;
    bbMin[2]--;

    InputImageType::RegionType tetraRegion;
    InputImageType::SizeType tetraRegionSize;
    tetraRegion.SetIndex(bbMin);
    tetraRegionSize[0] = bbMax[0]-bbMin[0];
    tetraRegionSize[1] = bbMax[1]-bbMin[1];
    tetraRegionSize[2] = bbMax[2]-bbMin[2];
    tetraRegion.SetSize(tetraRegionSize);
    tetraRegion.SetSize(tetraRegionSize);

    double dp0[3], dp1[3], dp2[3], dp3[3];
      
    memcpy(&dp0[0], meshDeformation->GetTuple3(thisCell->GetPointId(0)), 
        sizeof(double)*3);
    memcpy(&dp1[0], meshDeformation->GetTuple3(thisCell->GetPointId(1)), 
        sizeof(double)*3);
    memcpy(&dp2[0], meshDeformation->GetTuple3(thisCell->GetPointId(2)), 
        sizeof(double)*3);
    memcpy(&dp3[0], meshDeformation->GetTuple3(thisCell->GetPointId(3)), 
        sizeof(double)*3);

    itk::ImageRegionIteratorWithIndex<InputImageType> imageI(ref, tetraRegion);
    itk::ImageRegionIteratorWithIndex<OutputImageType> dfImageI(df, tetraRegion);

    imageI.GoToBegin(); dfImageI.GoToBegin();
    for(;!imageI.IsAtEnd();++imageI,++dfImageI){
      InputImageType::PointType point;
      ref->TransformIndexToPhysicalPoint(imageI.GetIndex(), point);
      double coord[3], closestPoint[3], pcoords[3], bc[4], dist;
      OutputImageType::PixelType dfPixel;
      int subId, isInside;
      coord[0] = -1.*point[0];
      coord[1] = -1.*point[1];
      coord[2] = point[2];

      isInside = tetra->EvaluatePosition(&coord[0], &closestPoint[0],
        subId, pcoords, dist, &bc[0]);

      if(!isInside)
        continue;

      dfPixel[0] = bc[0]*dp0[0]+bc[1]*dp1[0]+bc[2]*dp2[0]+bc[3]*dp3[0];
      dfPixel[1] = bc[0]*dp0[1]+bc[1]*dp1[1]+bc[2]*dp2[1]+bc[3]*dp3[1];
      dfPixel[2] = bc[0]*dp0[2]+bc[1]*dp1[2]+bc[2]*dp2[2]+bc[3]*dp3[2];

      OutputImageType::PixelType dfInverseVector;
      OutputImageType::PointType dfInversePoint;
      OutputImageType::IndexType dfInverseIndex;

      df->TransformIndexToPhysicalPoint(imageI.GetIndex(), dfInversePoint);
      dfInversePoint[0] += dfPixel[0];
      dfInversePoint[1] += dfPixel[1];
      dfInversePoint[2] += dfPixel[2];

      dfInverseVector[0] = -1.*dfPixel[0];
      dfInverseVector[1] = -1.*dfPixel[1];
      dfInverseVector[2] = -1.*dfPixel[2];
      df->TransformPhysicalPointToIndex(dfInversePoint, dfInverseIndex);

      dfInverse->SetPixel(dfInverseIndex, dfInverseVector);
      
      dfImageI.Set(dfPixel);
      mask->SetPixel(imageI.GetIndex(), 1);
    }
  }

#endif

  {
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(df);
  writer->SetFileName(dfImageName.c_str());
  writer->SetUseCompression(1);
  writer->Update();
  }

  {
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(dfInverse);
  writer->SetFileName(dfInverseImageName.c_str());
  writer->SetUseCompression(1);
  writer->Update();
  }

  if(0){
    MaskWriterType::Pointer writer = MaskWriterType::New();
    writer->SetInput(mask);
    writer->SetFileName("/Users/fedorov/Temp/mask.nrrd");
    writer->SetUseCompression(1);
    writer->Update();
  }

  return 0;
}
