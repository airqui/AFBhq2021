for pol in "eL_pR" "eR_pL"
do
    hadd -f -j 4 results/jettag_cuts6_2f_hadronic_sample_${pol}_250GeV.root output/jettag_cuts6_2f_hadronic_sample_${pol}_file_*_250GeV.root
    hadd -f -j 4 results/jettag_cuts6_2f_hadronic_${pol}_250GeV.root output/jettag_cuts6_2f_hadronic_${pol}_file_*_250GeV.root
done
