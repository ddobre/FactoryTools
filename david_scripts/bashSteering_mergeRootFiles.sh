rm -r merged_weighted_root
mkdir merged_weighted_root
total_string=" "
echo finished initializing
for file in weighted_root_files/user.ddobre*; do
  total_string+=" weighted_root_files/$(basename $file)"
done 
rm merged_weighted_root/root_output_merged.root
hadd merged_weighted_root/root_output_merged.root $total_string
