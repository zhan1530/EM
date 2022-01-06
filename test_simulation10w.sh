#! /bin/sh

input_dir='test10'
ref='test10/EF1813.fasta'
#count=0
#########################
##modified from here
NN=100
f1='test10/EF18130.fasta'
n1=(32500 32343 32354 32566 32661 32663 32483 32541 32757 32711 32436 32714 32710 32634 32595 32604 32593 32479 32366 32636 32680 32569 32629 32538 32660 32574 32599 32369 32349 32449 32630 32455 32531 32266 32564 32630 32388 32515 32513 32576 32463 32566 32587 32749 32611 32654 32400 32294 32822 32509 32403 32802 32405 32360 32361 32555 32429 32579 32497 32666 32649 32524 32673 32420 32668 32562 32668 32650 32534 32449 32599 32793 32638 32496 32414 32378 32255 32522 32466 32842 32522 32699 32470 32676 32184 32544 32685 32483 32275 32590 32618 32453 32208 32433 32678 32495 32492 32589 32565 32507)
f2='test10/EF18131.fasta'
n2=(15140 14991 15193 15005 15063 15062 15170 14993 15063 14936 14973 15068 15065 14967 15126 15152 14971 15019 14989 14931 15092 15120 15114 15030 14934 15035 15176 15108 15343 14908 15081 15239 15028 15067 15162 15083 15167 15020 15035 14825 15079 15135 14922 14882 14844 15007 15036 15117 14791 15205 15160 15023 15106 15156 15209 15116 15168 15058 15115 14998 15038 15170 15096 15308 15101 15064 15044 15190 14940 14926 14899 15088 14827 15087 15190 14941 15257 15138 15026 15067 15087 15109 15175 15149 15246 15227 15124 15132 14883 15236 14991 15079 15039 15259 15006 15077 15062 15107 15101 15159)
f3='test10/EF18132.fasta'
n3=(16145 16124 16028 16227 15954 15925 16115 16166 16036 16084 16259 15950 15866 16181 16079 16122 16214 16145 16265 16139 15978 15954 16053 16238 16212 16441 15817 16117 16268 16235 16103 15908 16136 16290 16171 16109 16168 16245 16111 16222 16280 16004 16241 16051 16214 16118 16122 15982 16240 16155 16179 15903 16060 16086 16070 16072 15942 16101 15932 16142 16037 16001 16077 16084 15921 15885 16124 16139 16005 16223 16157 15751 16126 16042 16179 15986 16090 15953 16006 16105 16255 16124 16134 16138 16295 15926 16064 16097 16246 16207 16227 16229 16262 16045 16065 16223 16154 16096 15791 16130)
f4='test10/EF18133.fasta'
n4=(8209 8169 8173 8204 8235 8203 8208 8189 8150 8324 8238 8221 8273 8192 8204 8013 8240 8335 8340 8316 8234 8091 8231 8112 8109 8225 8124 8283 8105 8213 8188 8314 8225 8242 8231 8343 8248 8060 8129 8172 8109 8252 8157 8239 8246 8102 8159 8399 8252 8040 8469 8202 8275 8437 8175 8164 8210 8197 8186 8231 8265 8336 8051 8177 8238 8310 8140 8100 8371 8175 8179 8335 8383 8292 8216 8205 8128 8364 8364 8134 8252 8355 8163 8073 8236 8430 8206 8218 8417 8252 8145 8037 8186 8148 8192 8296 8143 8214 8192 8200)
f5='test10/EF18134.fasta'
n5=(28006 28373 28252 27998 28087 28147 28024 28111 27994 27945 28094 28047 28086 28026 27996 28109 27982 28022 28040 27978 28016 28266 27973 28082 28085 27725 28284 28123 27935 28195 27998 28084 28080 28135 27872 27835 28029 28160 28212 28205 28069 28043 28093 28079 28085 28119 28283 28208 27895 28091 27789 28070 28154 27961 28185 28093 28251 28065 28270 27963 28011 27969 28103 28011 28072 28179 28024 27921 28150 28227 28166 28033 28026 28083 28001 28490 28270 28023 28138 27852 27884 27713 28058 27964 28039 27873 27921 28070 28179 27715 28019 28202 28305 28115 28059 27909 28149 27994 28351 28004)
##########################
Rs='nodirectbt.R'
Rno='directbt.R'
cppp='EMdirre.cpp'
cpp='EMnodirect.cpp'
py='cleanliang.py'
for iter in $(seq 1 ${NN});do
mkdir ${input_dir}/${iter}
	#############################
python ${input_dir}/program_rev.py -i $f1 -n ${n1[$iter]}
python ${input_dir}/program_rev.py -i $f2 -n ${n2[$iter]}
python ${input_dir}/program_rev.py -i $f3 -n ${n3[$iter]}
python ${input_dir}/program_rev.py -i $f4 -n ${n4[$iter]}
python ${input_dir}/program_rev.py -i $f5 -n ${n5[$iter]}
	###############################
cat ${input_dir}/*_left.txt > ${input_dir}/${iter}/combined_left.txt
cat ${input_dir}/*_right.txt > ${input_dir}/${iter}/combined_right.txt
rm ${input_dir}/*_left.txt ${input_dir}/*_right.txt
rm ${input_dir}/*location.txt 
bwa mem -w 0 -t 4 ${ref} ${input_dir}/${iter}/combined_left.txt ${input_dir}/${iter}/combined_right.txt > ${input_dir}/${iter}/aln.sam
    cp ${input_dir}/${py} ${input_dir}/${iter}/${py}
	cp ${input_dir}/${Rs} ${input_dir}/${iter}/${Rs}
    cp ${input_dir}/${Rno} ${input_dir}/${iter}/${Rno}
    cp ${input_dir}/${cpp} ${input_dir}/${iter}/${cpp}
    cp ${input_dir}/${cppp} ${input_dir}/${iter}/${cppp}
	cd ${input_dir}/${iter}
    python ${py}
    Rscript ${Rs}
    Rscript ${Rno}
	cd ../..
   
	rm ${input_dir}/${iter}/${Rs}
    rm ${input_dir}/${iter}/${Rno}
    rm ${input_dir}/${iter}/${cpp}
    rm ${input_dir}/${iter}/${cppp}
    rm ${input_dir}/${iter}/${py}
#rm ${input_dir}/${iter}/sim_aln.sam
    rm ${input_dir}/${iter}/fill1all
    rm ${input_dir}/${iter}/fill2all
    rm ${input_dir}/${iter}/fill3all
    rm ${input_dir}/${iter}/fill5all
	####################
	#####progress######
	###################
	#if ${iter}*100/${NN} == count;then
		#echo "######################"
		echo "${iter} iters have been done"
		#echo "################"
	#	count=count+10
	#fi
done

#while IFs= read -r line;do
#	for f in ${input_dir}/*fasta;do
#		echo "f $f N $line\n"
#		python program_rev.py -i $f -n $line
#	cat ${input_dir}/*${line}_left.txt > ${f}_${line}_combined_left.txt
#	cat ${input_dir}/*${line}_right.txt > ${f}_${line}_combined_right.txt
#	done
#done < num_list

