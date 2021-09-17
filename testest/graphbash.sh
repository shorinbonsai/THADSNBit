

awk '{printf "%s %s\n", NR-1,$0}' input.txt > new.txt