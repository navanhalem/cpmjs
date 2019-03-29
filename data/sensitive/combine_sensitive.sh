for d in $(find data/sensitive/* -maxdepth 3 -type d)
do
	for f in $d/*.txt
	do 
		cat $f | tail -1 | awk '{print $0, "sensitive"}' >> data/sensitive/sensitive_summary.txt
	done 
done
