#!/bin/bash
#this script takes all question files and appends them to a question group for limesurvey.
#call it with arguments: $1=surveyid, $2=groupid, $3=grouporder, $4=questionid
#
#
echo Please enter name of stimulus file
read filename
echo Please enter number for  survey id
read sid
echo Please enter number for group id
read gid
echo Please enter number for group order
read gord
echo Please enter number for question id
read qid
echo Thank you! The file will be created
#loop through question files
qord=0
#create backup of template_group as it will be overwritten
cp template_group.lsg template_group.lsg-copy
count=1
pausepoint=80
while read line
	do
	echo $count
		q=$(echo $line | awk '{print $1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9}') 
		code=$(echo $line | awk '{print $12$13}')
		if [ $count -eq $pausepoint ]
		then
			cat template_pause.txt | sed -e "s:enterquestionidhere:$qid:;s:enterquestionorderhere:$qord:" > tmppause.txt	
			awk '/enterquestioncodehere/{while(getline line<"tmppause.txt"){print line}} //' template_group.lsg > tmp
			mv tmp template_group.lsg
			qid=$(($qid+1))
			qord=$(($qord+1))
		echo $count
		fi 
		cat template_question.txt | sed -e "s:enterquestionidhere:$qid:;s:entercodehere:$code:;s:enterquestionhere:$q:;s:enterquestionorderhere:$qord:" > tmpquestion.txt
		cat template_answer.txt | sed -e "s:enterquestionidhere:$qid:" > tmpanswer.txt
		cat template_attribute.txt | sed -e "s:enterquestionidhere:$qid:" > tmpattribute.txt
		awk '/enterquestioncodehere/{while(getline line<"tmpquestion.txt"){print line}} //' template_group.lsg > tmp
		mv tmp template_group.lsg
			awk '/enteranswercodehere/{while(getline line<"tmpanswer.txt"){print line}} //' template_group.lsg > tmp
		mv tmp template_group.lsg
		awk '/enterattributecodehere/{while(getline line<"tmpattribute.txt"){print line}} //' template_group.lsg > tmp
		mv tmp template_group.lsg
		count=$(($count+1))
		qord=$(($qord+2))
		qid=$(($qid+3))
	done <  $filename
		cat template_group.lsg | sed -e "s:entersurveyidhere:$sid:;s:entergroupidhere:$gid:;s:enterorderhere:$gord:" > newGroup.lsg
	# replace the overwritten template by the original template again
	mv template_group.lsg-copy template_group.lsg

