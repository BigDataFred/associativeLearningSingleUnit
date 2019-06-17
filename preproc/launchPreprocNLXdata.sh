#!/bin/bash/

pID=$1
ExpMode=$2
#basepath=$3
#nworkers=$3

echo "dataset id is: "${pID}
echo "exp-mode is: "${ExpMode}

basepath="/media/rouxf/rds-share/Archive/MICRO/"${pID}"/"${ExpMode}"/"
savepath="/media/rouxf/rds-share/iEEG_DATA/MICRO/"${pID}"/"${ExpMode}"/"

dum="/media/rouxf/rds-share/iEEG_DATA/MICRO/"${pID}
if [ ! -d "${dum}" ];then
    mkdir "${dum}"
fi

dum="/media/rouxf/rds-share/iEEG_DATA/MICRO/"${pID}"/"${ExpMode}
if [ ! -d "${dum}" ];then
    mkdir "${dum}"
fi

sesh=$( ls $basepath)
#sesh=$( echo ${sesh} | grep -oi "201[0-9]-[0-9][0-9]-[0-9][0-9]_[0-9][0-9]-[0-9][0-9]-[0-9][0-9]")

#sesh="2016-12-09_14-49-14"

Destination="/home/rouxf/in/"${ExpMode}"/"

if [ ! -d "${Destination}" ];then
    mkdir -v "${Destination}"
fi

Source2="/home/rouxf/out/"${ExpMode}"/"

if [ ! -d "${Source2}" ];then
    mkdir -v "${Source2}"
fi

c=0
for i in ${sesh}
do
    c=$((c+1))
    Source=${basepath}${i}
    echo "session folder "${c}":"${i}

    if [ ! -d "${Destination}${i}" ];then
        mkdir "${Destination}${i}"
    fi

    if [ ! -d "${Source2}${i}" ];then
        mkdir "${Source2}${i}"
    fi

    #cp -rv ${Source} ${Destination}${i}
    cp -rv "${Source}" "${Destination}"
    
    cd "${Destination}${i}"
    for ncsFiles in *.ncs
    do
        echo "${ncsFiles}"
	fS=$( head -n 15 ${ncsFiles} | grep "SamplingFrequency" | awk '{print $2}')
	echo "${fS}"
	fS=${fS//[$'\001'-$'\037']}
	echo "${fS}"
	if [ $((${fS}+0)) -ge 32000 ] 
	then 
		echo "keeping file" ${ncsFiles}
	else 
		echo "deleting file" ${ncsFiles}
		rm -v "${ncsFiles}"

	 fi
    done
		
    cd /home/rouxf/in/
	
    #command="try;analyseEMtest2('${pID}','${ExpMode}','${i}',${nworkers});catch;delete(gcp);exit;end;"
    #command="try;analyze_EM_data2('${pID}','${ExpMode}');catch;delete(gcp);exit;end;"
    cp -vf ~/prj/Bham/code/mcode/project_EM/analyze_EM_data2.m ~/in/
    command="analyze_EM_data2('${pID}','${ExpMode}');"
    echo "${command}"
    matlab -nosplash -nodesktop -r ${command} 

    rm -rv "${Destination}${i}"

    if [ ! -d "${savepath}${i}" ];then
        mkdir "${savepath}${i}"
    fi

    #cp -rfv ${Source2}${i} ${savepath}${i} 
    cp -rfv "${Source2}${i}" "${savepath}"
    rm -rv "${Source2}${i}"

done

rm -r "${Destination}"

