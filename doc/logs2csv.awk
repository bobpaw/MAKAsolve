#!/usr/bin/awk -f

BEGIN {
	print("Ranks,LocalNodes,TotalNodes,Setup,SolverConstructor,Integration," \
		"Assembly,Solve,SolverTotal,Total")
}

/^ranks num_proc_nodes.*/ {
	getline
	print($1 "," $2 "," $3 "," $4 "," $5 "," $6 "," $7 "," $8 "," \
		$5 + $6 + $7 + $8 "," $4 + $5 + $6 + $7 + $8)
}
