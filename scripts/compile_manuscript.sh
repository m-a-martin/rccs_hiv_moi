i=$1
lualatex $i
bibtex ${i%".tex"}
lualatex $i
lualatex $i