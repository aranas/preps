#!/bin/bash
echo Please enter path:
read pathname
FILES="$pathname"results*
###
#for each subject transpose response file and sort responses according to items.
###
for f in $FILES
do
filename="${f%.*}"
awk -F ',' '{for (i=1; i<=NF; i++) {a[NR,i] = $i}} NF>p { p = NF} END{
	for(j=1; j<=p; j++){
		str=a[1,j]" "a[2,j];
		print str
	}
}' $f > "$filename"_converted.txt
sort -k1,1 "$filename"_converted.txt -o  "$filename"_converted.txt
done
FILESconv="$pathname"*_converted.txt
###
#combine all subject's response files together
###
for f in $FILESconv
do
join -a 2 -e "-" empty.txt $f > "$pathname"responses_all.txt
mv "$pathname"responses_all.txt empty.txt
done
mv empty.txt "$pathname"responses_all.txt
echo "" > empty.txt
#awk '$0 !~ /#/{arr[$1]=arr[$1]" " $2}END{for(i in arr) print i, arr[i]}' $FILESconv  > "$pathname"responses_all.txt
#sort -k1 -n "$pathname"responses_all.txt -o "$pathname"responses_all.txt
header="surveynum "$FILESconv
echo $header > headerfile.txt
cat headerfile.txt "$pathname"responses_all.txt > tmp.txt
mv tmp.txt "$pathname"responses_all.txt
#awk '{ a[FNR]= (a[FNR] ? a[FNR] FS : "") $2 } END { for(i=1;i<=FNR;i++) print a[i] }' $FILESconv > "$pathname"responses_allsorted.txt
