for process in "ISR"
do 
#"ww_240" "ISR_240" "HV_240_mDv120_mqv60" "HV_240_mDv100_mqv50" 
    for cut in 0 1 2 3 4 5 6 20
    do
	hadd -f -j 4 -k results_FCC/selection_cuts${cut}_${process}.root output/selection_cuts${cut}_*${process}*.root
	#hadd -f -j 4 -k results_cut5_6_emax80/selection_cuts${cut}_eetoDvDvbar_m62.5.root output/selection_cuts${cut}_eetoDv*62*.root
	#hadd -f -j 4 -k results_cut5_6_emax80/selection_cuts${cut}_eetoZv_m62.5.root output/selection_cuts${cut}_eetoZv*62*.root
	#hadd -f -j 4 -k results_FCC/QCDcorrelations_cuts${cut}_${process}.root output/QCDcorrelations_cuts${cut}_*${process}*.root
    done
done


