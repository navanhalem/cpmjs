paramfile=makefiles/params/params_$2_$3.txt
nsim=$4
time=$1

# ---------------------------------------------------------------------
# CODE:
# ---------------------------------------------------------------------

np=$( cat $paramfile | wc -l )

echo ".DELETE_ON_ERROR :"
echo "all : "
# Loop over the simulations and parameter combis

for p in $(seq 1 $np) ; do

	LAMBDACHEM=$( cat $paramfile | awk -v line=$p 'NR==line{print $1}')
	KILLINGTIME=$( cat $paramfile | awk -v line=$p 'NR==line{print $2}')
	BIASEDENTRY=$( cat $paramfile | awk -v line=$p 'NR==line{print $3}')
	EXPNAME=$( cat $paramfile | awk -v line=$p 'NR==line{print $4}')
	TYPE=3	

	if [ "$EXPNAME" == "sensitive" ] 
		then
			TYPE=1
	fi
	if [ "$EXPNAME" == "arresting" ] 
		then
			TYPE=2
	fi

	NAME=$EXPNAME-lambdachem$LAMBDACHEM-killingtime$KILLINGTIME-biasedentry$BIASEDENTRY

	# Now the recipes for the individual simulation tracks
	for sim in $(seq 1 $nsim) ; do

		# Ensure the loop can be easily stopped
		trap "exit 1" SIGINT SIGTERM

		# trackfiles
		FILE=data/$EXPNAME/$KILLINGTIME/$NAME-sim$sim.txt
		echo "$FILE : skin/epidermis0.js"
		echo -e "\t@"node \$\< $time 0 600 $LAMBDACHEM 0.025 $BIASEDENTRY $KILLINGTIME $TYPE 0.00002 20 1000 2 "> \$@"
		echo "all : "$FILE
	done
done
