for d in $(find data/arresting/* -maxdepth 3 -type d)
do
	for f in $d/*.txt
	do 
		cat $f | tail -1 | awk '{print $0, "arresting"}' >> data/arresting/arresting_summary.txt
	done 
done
