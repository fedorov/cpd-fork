for c in 9 10
do
  for partial in 1 0
  do
    for rigid in 1 0
    do
      bash Scripts/bwh_launcher.sh $c $partial $rigid
    done
  done
done
