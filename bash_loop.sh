#!/bin/bash
nparams="c('ZTOFD','GKWAKE','GKDRAG','GFRCRIT')"
params_exclude='"NULL" "GKDRAG" "GFRCRIT"'
#params_exclude="NULL"
vars="ua"
levels="c(5000,25000,85000)"
sqorder='"linear" "mixed"'
#sqorder="quadratic"
seasons="c('DJFM')" # c('DJFM','JJAS')"
centering='"TRUE" "FALSE"'
#seasons="DJFM"

for pex in $params_exclude ; do
	for season in $seasons ; do
		for order in $sqorder ; do
			for center in $centering ; do
		echo $pex $season $order $center
		logfile=logfiles/sloop_${pex}_${season}_${order}_${center}.txt
		rm -f $logfile
		Rscript new_fingerprint.R $vars $nparams $pex $levels $order $season ${center} >> $logfile 2>&1 &
	done
done
done
done
wait
#Rscript new_fingerprint.R $nparams $pex $levels $sqorder $seasons
