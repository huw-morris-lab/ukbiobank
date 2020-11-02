for n in {0..151}
do
	echo "qsub -pe make 2 -cwd GWASscript_$n.r"
done > master_GWASscript_for_qsub.sh