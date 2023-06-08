for pol in "eL_pR" "eR_pL"
do
    for sample in "2f_hadronic" "2f_hadronic_sample" "qqH" "4f_ZZ_hadronic" "4f_WW_hadronic" "6f_ttbar"
    do
		for quark in 4 5 
		do
			hadd -f -j 4 results_500GeV_dNdx/SigBkg_pdg${quark}_${sample}_${pol}.root output/SigBkg_quark${quark}_${sample}_${pol}*root
		done
    done
done
