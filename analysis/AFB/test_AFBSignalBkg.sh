process=$1
pol=$2
method=$3
bkg=$4
folder="/lustre/ific.uv.es/prj/ific/flc/ntuples-2020/"${process}"_"${pol}"/"
local=$PWD


for quark in 5 4
do
    counter=0
    
    for file in ${folder}/*
    do
	echo $file $counter
	name="00"$counter
	if [ $counter -gt 9 ]; then
	    name="0"$counter
	fi
	
	if [ $counter -gt 99 ]; then
	    name=$counter
	fi
#int AFBAnalysis( TString file="", TString process="2f_hadronic", TString pol="eL", int file_n=0, float Kv=35, int quark=4, int method=0){
	
	cat > ${local}/steer/afbsignal_quark${quark}_${process}_${pol}_${name}_${method}.sh <<EOF
source ${local}/../init_ilcsoft.sh
root -l ${local}/AFBSignalBkg.cc\(\"${file}\",\"${process}\",\"${pol}\",${counter},35,${quark},${method},${bkg}\) > ${local}/output/log_AFBSignalBkg_quark${quark}_${process}_${pol}_${name}_${method}_${bkg}
mv SigBkg_quark${quark}_${process}_${pol}_file_${name}_250GeV_cheatmethod_${method}.root ${local}/output/.
EOF
	
	cat > ${local}/steer/afbsignal_quark${quark}_${process}_${pol}_${name}_${method}.sub <<EOF
# Unix submit description file
# kt_xNAMEfile.sub -- simple Marlin job
executable              = ${local}/steer/afbsignal_quark${quark}_${process}_${pol}_${name}_${method}.sh
log                     = ${local}/log/afbsignal_quark${quark}_${process}_${pol}_${name}_${method}.log
output                  = ${local}/log/outfile_afbsignal_quark${quark}_${process}_${pol}_${name}_${method}.txt
error                   = ${local}/log/errors_afbsignal_quark${quark}_${process}_${pol}_${name}_${method}.txt
should_transfer_files   = Yes
when_to_transfer_output = ON_EXIT
queue 1
EOF
      	
	if [ -f ${local}/output/SigBkg_quark${quark}_${process}_${pol}_file_${name}_250GeV_cheatmethod_${method}.root ];
	then
            echo "Skip quark${quark}_${process}_${pol}_${name}_${method}_${method}"
	else	
	    condor_submit ${local}/steer/afbsignal_quark${quark}_${process}_${pol}_${name}_${method}.sub
	fi
	
	counter=$((counter+1))
    done
done



