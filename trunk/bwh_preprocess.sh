mesher=/Users/fedorov/github/cpd/trunk/SurfaceMesher-build/lib/Slicer-4.3/cli-modules/QuadEdgeSurfaceMesher
clipper=/Users/fedorov/github/cpd/trunk/ClipApexBase-build/lib/Slicer-4.3/cli-modules/ClipApexBase
data=/Users/fedorov/Documents/Projects/BRP/MR-US-registration
decimation=0.01

#for c in 9 10 12 14 16 17 18 19 20 21 22 23
for c in 10 
do
  $mesher ${data}/Case${c}/SmoothReg/case${c}-MR-smooth.nrrd --decimation ${decimation} ${data}/Case${c}/SmoothReg/case${c}-MR-smooth.ply
  $mesher ${data}/Case${c}/SmoothReg/case${c}-US-smooth.nrrd --decimation ${decimation} ${data}/Case${c}/SmoothReg/case${c}-US-smooth.ply
  $clipper ${data}/Case${c}/SmoothReg/case${c}-US-smooth.ply 10 10 ${data}/Case${c}/SmoothReg/case${c}-US-smooth-cut10.ply
  mkdir ${data}/Case${c}/CPD_registration
done
