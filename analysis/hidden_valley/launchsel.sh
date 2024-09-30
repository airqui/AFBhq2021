
#sleep 3h

source ../init_ilcsoft.sh
#source eliminar.sh

#cp PreSelectionHidden.h ../common/.

#source test_selection_HV.sh ISR
#source test_selection_HV.sh HV_qv100GeV
#source test_selection_HV.sh HV_qv50GeV 
#source test_selection_HV.sh HV_qv10GeV 
#source test_selection_HV.sh qqSM 
#source test_selection_HV.sh qqbar_Pythia
#source test_selection_HV.sh qqbar_Herwig

#sleep 5m
#mkdir results_cut6_npfo22_nchpfo15
#source haddselection_2.sh results_cut6_npfo22_nchpfo15

#source launchCorrel.sh


#source test_selection_HV.sh eetoDvDvbar_mDV100
#source test_selection_HV.sh eetoDvDvbar_mDV62.5

#source test_selection_HV.sh HV_240_mDv120_mqv100
#source test_selection_HV.sh HV_240_mDv120_mqv60
#source test_selection_HV.sh HV_240_mDv120_mqv10
#source test_selection_HV.sh HV_240_mDv100_mqv50
#source test_selection_HV.sh ISR_240
source test_selection_HV.sh ISR

#source test_QCDcorrelations_HV.sh HV_240_mDv120_mqv100
#source test_QCDcorrelations_HV.sh HV_240_mDv120_mqv10
#source test_QCDcorrelations_HV.sh ISR_240

#source launchCorrel.sh

sleep 1m

source haddselection.sh
