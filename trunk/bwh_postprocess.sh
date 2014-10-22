mesher=/Users/fedorov/github/cpd/trunk/SurfaceMesher-build/lib/Slicer-4.3/cli-modules/QuadEdgeSurfaceMesher
clipper=/Users/fedorov/github/cpd/trunk/ClipApexBase-build/lib/Slicer-4.3/cli-modules/ClipApexBase
dfconverter=/Users/fedorov/github/cpd/trunk/MeshDisplacementToDF-build/lib/Slicer-4.3/cli-modules/MeshDisplacementToDF
data=/Users/fedorov/Documents/Projects/BRP/MR-US-registration
decimation=0.01

#for c in 9 10 12 14 16 17 18 19 20 21 22 23

femPrefix='rigid'

for c in 10 
do
  $dfconverter ${data}/Case${c}/SmoothReg/case${c}-US-smooth.nrrd /Users/fedorov/Documents/Projects/BRP/MR-US-registration/Case${c}/CPD_registration/case${c}_${femPrefix}_fem_result_mesh_Partial.vtk /Users/fedorov/Documents/Projects/BRP/MR-US-registration/Case10/CPD_registration/case${c}_${femPrefix}_fem_fixedToMoving_DF.nrrd /Users/fedorov/Documents/Projects/BRP/MR-US-registration/Case10/CPD_registration/case${c}_${femPrefix}_fem_movingToFixed_DF.nrrd

done
