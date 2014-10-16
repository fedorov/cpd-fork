#include <iostream>

#include "vtkPLYReader.h"
#include "vtkPLYWriter.h"
#include "vtkPolyData.h"
#include "vtkDataSet.h"
#include "vtkSmartPointer.h"

#include "vtkClipPolyData.h"
#include "vtkPlane.h"

#include <ClipApexBaseCLP.h>

int main(int argc, char* argv[]){
  PARSE_ARGS;

  // read the input mesh
  vtkSmartPointer<vtkPLYReader> meshReader = 
    vtkSmartPointer<vtkPLYReader>::New();
  meshReader->SetFileName(meshName.c_str());
  meshReader->Update();
  vtkSmartPointer<vtkPolyData> mesh = meshReader->GetOutput();

  double bounds[6];
  mesh->GetBounds(bounds);

  vtkSmartPointer<vtkClipPolyData> clipper1 = 
    vtkSmartPointer<vtkClipPolyData>::New();

  vtkSmartPointer<vtkPlane> plane1 =
    vtkSmartPointer<vtkPlane>::New();
  double pn1[3], po1[3];
  pn1[0] = 0;
  pn1[1] = 0;
  pn1[2] = 1;
  po1[0] = 0;
  po1[1] = 0;
  po1[2] = bounds[4]+offset1;
  plane1->SetNormal(pn1);
  plane1->SetOrigin(po1);

  vtkSmartPointer<vtkClipPolyData> clipper2 = 
    vtkSmartPointer<vtkClipPolyData>::New();

  vtkSmartPointer<vtkPlane> plane2 =
    vtkSmartPointer<vtkPlane>::New();
  double pn2[3], po2[3];
  pn2[0] = 0;
  pn2[1] = 0;
  pn2[2] = -1;
  po2[0] = 0;
  po2[1] = 0;
  po2[2] = bounds[5]-offset2;
  plane2->SetNormal(pn2);
  plane2->SetOrigin(po2);

  clipper1->SetClipFunction(plane1);
  clipper1->SetInputConnection(meshReader->GetOutputPort());

  clipper2->SetClipFunction(plane2);
  clipper2->SetInputConnection(clipper1->GetOutputPort());
  clipper2->Update();

  vtkSmartPointer<vtkPLYWriter> meshWriter = 
    vtkSmartPointer<vtkPLYWriter>::New();
  meshWriter->SetFileName(outputMeshName.c_str());
  meshWriter->SetInputConnection(clipper2->GetOutputPort());
  meshWriter->Update();
  
  return 0;
}
