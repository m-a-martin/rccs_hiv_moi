# replace values with actual values
i=$1
k=${i%".tex"}_vals.tex
scp $i $k
while read p; do
	var=$(echo $p | cut -d';' -f1)
	val=$(echo $p | cut -d';' -f2- | sed 's/"//g'| sed 's/\%/\\%/g')
	sed "s/\\\var{$var}/$val/g" $k > tmp.tex
	mv tmp.tex $k
done < $2

# manually replace some values
sed 's/1.0\times/1.0\\times/g' $k > tmp.tex
mv tmp.tex $k
sed 's/\\text{obs}/obs/g' $k > tmp.tex
mv tmp.tex $k
#sed 's/\\MI/M^{obs}/g' $k > tmp.tex
#mv tmp.tex $k
# replace namerefs with values
# don't think pandoc supports namerefs
# get a list of all namerefs
cat <(echo "Materials and methods;Materials and methods")\
	<(paste -d";" \
		<(sed -n -e '/Supporting information/,$p' $i | \
			grep label | sed 's/\\label{//g' | \
				sed 's/}//g' | sed 's/{//g')  \
		<(sed -n -e '/Supporting information/,$p' $i | \
			grep paragraph | sed 's/\\paragraph\*{//g' | \
				sed 's/}//g' | sed 's/{//g' | sed 's/\.//g')) \
	> tmp_nameref.scsv
while read p; do
	var=$(echo $p | cut -d';' -f1)
	val=$(echo $p | cut -d';' -f2-)
	sed "s/\\\nameref{$var}/$val/g" $k > tmp.tex
	mv tmp.tex $k
done < tmp_nameref.scsv
rm -rf tmp_nameref.scsv
# remove equation references in parenthesis 
sed "s/ (Eq[^)]*)//g" $k > tmp.tex
mv tmp.tex $k
# finally add references header
sed 's/\\bibliography/\\section{References}\n\\bibliography/g' $k > tmp.tex
mv tmp.tex $k
# remove booleans so figures print
sed "s/\\\ifthenelse{\\\boolean{includefigs}}{//g" $k | \
	sed "s/}{}//g" > tmp.tex
mv tmp.tex ${k%".tex"}_figs.tex
for g in $k ${k%".tex"}_figs.tex; do
	echo $g
	o=${g%".tex"}.docx
	pandoc --biblatex --bibliography=./references.bib -s $g -o $o --citeproc
	mv $o ${o/_vals/}
done
# clean up
#for g in $k ${k%".tex"}_figs.tex; do
#	rm $g
#done

echo "done here"
bash ../scripts/compile_manuscript.sh $i

echo "done all"

#latexdiff previous_versions/2024-10-23/rccs_hiv_multiple_infections_vals.tex rccs_hiv_multiple_infections_vals.tex \
#	> compare_2024-10-23_2025-01-30.tex

#lualatex compare_2024-10-23_2025-01-30.tex
#bibtex compare_2024-10-23_2025-01-30
#lualatex compare_2024-10-23_2025-01-30.tex
#lualatex compare_2024-10-23_2025-01-30.tex
#mv compare_2024-10-23_2025-01-30* previous_versions/compare_2024-10-23_2025-01-30/.
