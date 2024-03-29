for pol in "eL_pR" "eR_pL"
do
    for cut in {0..6}
    do
	hadd -f -j 4 results/selection_cuts${cut}_4f_ZZ_hadronic_${pol}_250GeV.root output/selection_cuts${cut}_4f_ZZ_hadronic_${pol}_*250GeV.root &
	hadd -f -j 4 results/selection_cuts${cut}_2f_hadronic_sample_${pol}_250GeV.root output/selection_cuts${cut}_2f_hadronic_sample_${pol}_*250GeV.root 
	hadd -f -j 4 results/selection_cuts${cut}_qqH_${pol}_250GeV.root output/selection_cuts${cut}_qqH_${pol}_*250GeV.root &
	hadd -f -j 4 results/selection_cuts${cut}_4f_WW_hadronic_${pol}_250GeV.root output/selection_cuts${cut}_4f_WW_hadronic_${pol}_*250GeV.root 
	hadd -f -j 4 results/selection_cuts${cut}_2f_hadronic_${pol}_250GeV.root output/selection_cuts${cut}_2f_hadronic_${pol}_*250GeV.root
    done
done



