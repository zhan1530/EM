#! /bin/sh

input_dir='test5'
ref='test5/EF1813.fasta'
#count=0
#########################
##modified from here
NN=1
f1='test5/EF18130.fasta'
n1=(317 334 330 344 337 329 308 318 316 334 309 332 343 335 333 349 302 327 328 315 330 347 320 343 335 338 330 331 329 317 326 345 303 337 324 310 345 311 324 340 311 323 329 327 316 327 314 317 332 328 331 333 313 319 317 328 330 335 338 324 310 326 303 306 327 357 340 351 323 334 351 324 308 287 321 316 297 323 329 318 328 312 303 303 312 343 310 313 330 331 298 356 308 311 302 311 316 303 337 296)

f2='test5/EF18131.fasta'
n2=(139 145 167 140 132 143 122 157 155 156 171 139 142 172 149 134 165 138 174 153 151 145 154 144 142 160 169 151 138 155 142 162 163 136 161 164 151 146 153 168 145 151 139 161 149 152 150 134 135 126 141 145 142 172 162 150 142 154 162 150 148 146 149 157 140 121 154 142 132 153 136 142 144 180 152 170 153 146 145 168 149 150 146 180 142 144 134 165 144 149 157 162 140 158 153 125 159 140 149 164)
f3='test5/EF18132.fasta'
n3=(170 162 148 149 183 168 184 151 138 165 172 172 169 139 145 152 150 172 154 155 173 145 159 153 161 151 167 147 166 161 161 132 147 167 163 137 171 148 154 177 155 173 162 155 177 151 167 172 137 164 160 156 154 167 154 173 156 162 129 191 175 172 179 170 167 186 177 154 168 157 171 158 160 159 174 152 162 156 166 160 152 148 165 154 186 159 173 157 149 153 191 151 173 173 166 164 148 175 163 149)

f4='test5/EF18133.fasta'
n4=(90  77  97  83  91  85  76  88  97  73  74  70  63  80  78  83  84  84  73  79  87  78  76  93  74  78  85 94  70  87  81  75  81  86  81  91  98  81  75  83  74  89  83  76  72  86  88  78  86 101  85  69 109  80 72  76  82  80  84  76  67  59  84  75  99  80  77  87  82  91  66  72  87  80  67  74  98  90  73  70  80 84  73  81  80  78  85  85 108  81  73  78  84  84  79  83  83  81  73 105)
f5='test5/EF18134.fasta'
n5=(284 282 258 284 257 275 310 286 294 272 274 287 283 274 295 282 299 279 271 298 259 285 291 267 288 273 249 277 297 280 290 286 306 274 271 298 235 314 294 232 315 264 287 281 286 284 281 299 310 281 283 297 282 262 295 273 290 269 287 259 300 297 285 292 267 256 252 266 295 265 276 304 301 294 286 288 290 285 287 284 291 306 313 282 280 276 298 280 269 286 281 253 295 274 300 317 294 301 278 286)
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
bwa index ${ref}
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
    #    count=count+10
    #fi
done

