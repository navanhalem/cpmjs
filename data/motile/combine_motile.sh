for d in $(find data/motile/* -maxdepth 3 -type d)
do
	for f in $d/*.txt
	do 
		cat $f | tail -1 | awk '{print $0, "motile"}' >> data/motile/motile_summary.txt
	done 
done
