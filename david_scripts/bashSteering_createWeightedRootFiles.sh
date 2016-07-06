rm -r weighted_root_files
mkdir weighted_root_files
# mkdir weighted_root_files/logs

for inputDirectory in output_pythiasamples/user.ddob*; do
  # echo 1st loop 
  # echo $inputDirectory 
  for file in $inputDirectory/*; do
    # echo 2nd loop 
    # echo $(basename $file) 
    python createWeightedRootFiles.py --inputDS ${inputDirectory} --inputFile $(basename $file) &  # > weighted_root_files/logs/${inputDirectory}.log & a
    # echo 
  done
done 
