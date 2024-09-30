for cut in 0 6 20 
do
    #hadd -f -j 8 -k results_pythia_herwig/QCDcorrelations_cuts${cut}_HV_qv100GeV.root output/QCDcorrelations_cuts${cut}_HV_qv100GeV_*.root
    hadd -f -j 4 -k results_202409/QCDcorrelations_cuts${cut}_eetoDvDvbar_m80.root output/QCDcorrelations_cuts${cut}_*80*.root
    #hadd -f -j 4 -k results_tests202408/QCDcorrelations_cuts${cut}_eetoDvDvbar_m62.5.root output/QCDcorrelations_cuts${cut}_*Dv*62*.root
    #hadd -f -j 4 -k results_tests202408/QCDcorrelations_cuts${cut}_eetoZv_m62.5.root output/QCDcorrelations_cuts${cut}_*Zv*62*.root
    #hadd -f -j 8 -k results_pythia_herwig/QCDcorrelations_cuts${cut}_HV_qv50GeV.root output/QCDcorrelations_cuts${cut}_HV_qv50GeV_*.root
    #hadd -f -j 8 -k results_pythia_herwig/QCDcorrelations_cuts${cut}_HV_qv10GeV.root output/QCDcorrelations_cuts${cut}_HV_qv10GeV_*.root
    #hadd -f -j 8 -k results_pythia_herwig/QCDcorrelations_cuts${cut}_qqSM.root output/QCDcorrelations_cuts${cut}_qqSM_*.root
    #hadd -f -j 8 -k results_pythia_herwig/QCDcorrelations_cuts${cut}_VBF.root output/QCDcorrelations_cuts${cut}_VBF_*.root

    #hadd -f -j 12 -k output/QCDcorrelations_cuts${cut}_ISR_merg_0.root output/QCDcorrelations_cuts${cut}_ISR_file_0*.root
    #hadd -f -j 12 -k output/QCDcorrelations_cuts${cut}_ISR_merg_1.root output/QCDcorrelations_cuts${cut}_ISR_file_1*.root
    #hadd -f -j 12 -k output/QCDcorrelations_cuts${cut}_ISR_merg_2.root output/QCDcorrelations_cuts${cut}_ISR_file_2*.root
    #hadd -f -j 12 -k output/QCDcorrelations_cuts${cut}_ISR_merg_3.root output/QCDcorrelations_cuts${cut}_ISR_file_3*.root
    #hadd -f -j 12 -k output/QCDcorrelations_cuts${cut}_ISR_merg_4.root output/QCDcorrelations_cuts${cut}_ISR_file_4*.root
    #hadd -f -j 12 -k output/QCDcorrelations_cuts${cut}_ISR_merg_5.root output/QCDcorrelations_cuts${cut}_ISR_file_5*.root
    #hadd -f -j 12 -k output/QCDcorrelations_cuts${cut}_ISR_merg_6.root output/QCDcorrelations_cuts${cut}_ISR_file_6*.root
    #hadd -f -j 12 -k output/QCDcorrelations_cuts${cut}_ISR_merg_7.root output/QCDcorrelations_cuts${cut}_ISR_file_7*.root
    #hadd -f -j 8 -k results_pythia_herwig/QCDcorrelations_cuts${cut}_ISR.root output/QCDcorrelations_cuts${cut}_ISR_merg*
    #rm output/QCDcorrelations_cuts${cut}_ISR_merg*

done



