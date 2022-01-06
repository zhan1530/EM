#! /bin/sh

input_dir='test3'
ref='test3/EF1813.fasta'
#count=0
#########################
##modified from here
NN=100
f1='test3/EF18130.fasta'
n1=(4567 4615 4566 4600 4586 4590 4494 4561 4490 4571 4434 4582 4584 4407 4658 4592 4572 4613 4554 4573 4632 4519 4642 4473 4543 4545 4637 4656 4592 4671 4648 4525 4614 4478 4663 4546 4618 4615 4545 4535 4583 4682 4546 4622 4601 4539 4560 4592 4575 4667 4521 4493 4627 4616 4492 4612 4565 4618 4583 4588 4622 4671 4666 4604 4604 4648 4579 4549 4699 4603 4590 4588 4565 4497 4533 4708 4621 4623 4567 4521 4605 4636 4560 4637 4652 4688 4714 4531 4535 4643 4847 4654 4524 4459 4431 4491 4582 4616 4518 4532)

f2='test3/EF18131.fasta'
n2=(9256 9229 9316 9171 9380 9499 9317 9318 9399 9428 9375 9318 9385 9452 9241 9315 9393 9437 9283 9223 9294 9395 9236 9336 9364 9262 9231 9272 9311 9331 9342 9245 9251 9402 9314 9335 9348 9250 9412 9395 9354 9284 9349 9304 9286 9185 9407 9248 9357 9249 9417 9359 9313 9372 9378 9458 9425 9314 9316 9336 9254 9306 9371 9231 9360 9250 9289 9293 9249 9352 9392 9328 9362 9402 9438 9308 9335 9389 9342 9351 9266 9147 9294 9347 9268 9283 9264 9343 9386 9391 9234 9250 9374 9476 9354 9432 9306 9334 9404 9283)
f3='test3/EF18132.fasta'
n3=(6177 6156 6118 6229 6034 5911 6189 6121 6111 6001 6191 6100 6031 6141 6101 6093 6035 5950 6163 6204 6074 6086 6122 6191 6093 6193 6132 6072 6097 5998 6010 6230 6135 6120 6023 6119 6034 6135 6043 6070 6063 6034 6105 6074 6113 6276 6033 6160 6068 6084 6062 6148 6060 6012 6130 5930 6010 6068 6101 6076 6124 6023 5963 6165 6036 6102 6132 6158 6052 6045 6018 6084 6073 6101 6029 5984 6044 5988 6091 6128 6129 6217 6146 6016 6080 6029 6022 6126 6079 5966 5919 6096 6102 6065 6215 6077 6112 6050 6078 6185)

##########################
Rs='nodirectbt.R'
Rno='directbt.R'
cpp='EMnodirect.cpp'
cppp='gene_with_direct.cpp'
#py='cleanliang.py'
for iter in $(seq 1 ${NN});do
mkdir ${input_dir}/${iter}
	#############################
python ${input_dir}/program_rev.py -i $f1 -n ${n1[$iter]}
python ${input_dir}/program_rev.py -i $f2 -n ${n2[$iter]}
python ${input_dir}/program_rev.py -i $f3 -n ${n3[$iter]}
#python ${input_dir}/program_rev.py -i $f4 -n $n4[$iter]
#python ${input_dir}/program_rev.py -i $f5 -n $n5[$iter]
	###############################
cat ${input_dir}/*_left.txt > ${input_dir}/${iter}/combined_left.txt
cat ${input_dir}/*_right.txt > ${input_dir}/${iter}/combined_right.txt
rm ${input_dir}/*_left.txt ${input_dir}/*_right.txt
rm ${input_dir}/*location.txt
bwa index ${ref}
bwa mem -w 0 -t 4 ${ref} ${input_dir}/${iter}/combined_left.txt ${input_dir}/${iter}/combined_right.txt > ${input_dir}/${iter}/aln.sam
    #cp ${input_dir}/${py} ${input_dir}/${iter}/${py}
	cp ${input_dir}/${Rs} ${input_dir}/${iter}/${Rs}
    cp ${input_dir}/${Rno} ${input_dir}/${iter}/${Rno}
    cp ${input_dir}/${cppp} ${input_dir}/${iter}/${cppp}
    cp ${input_dir}/${cpp} ${input_dir}/${iter}/${cpp}
	cd ${input_dir}/${iter}
    #python ${py}
    Rscript ${Rs}
    Rscript ${Rno}
	cd ../..
   
	rm ${input_dir}/${iter}/${Rs}
    rm ${input_dir}/${iter}/${Rno}
    rm ${input_dir}/${iter}/${cpp}
    rm ${input_dir}/${iter}/${cppp}
    #rm ${input_dir}/${iter}/${py}
#rm ${input_dir}/${iter}/sim_aln.sam
   #rm ${input_dir}/${iter}/fill1all
   #rm ${input_dir}/${iter}/fill2all
   #rm ${input_dir}/${iter}/fill3all
   #rm ${input_dir}/${iter}/fill5all
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

