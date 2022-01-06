#! /bin/sh

input_dir='test5_1w'
ref='test5_1w/EF1813.fasta'
#count=0
#########################
##modified from here
NN=50
f1='test5_1w/EF18130.fasta'
n1=(3231 3181 3289 3223 3221 3281 3230 3285 3237 3236 3248 3220 3187 3237 3245 3233 3302 3241 3193 3205 3247 3245 3224 3266 3137 3306 3262 3219 3226 3271 3151 3239 3221 3232 3254 3212 3219 3265 3225 3278 3199 3238 3228 3200 3270 3238 3211 3291 3217 3290 3194 3249 3265 3300 3234 3188 3237 3343 3185 3170 3133 3250 3264 3297 3249 3273 3212 3255 3289 3230 3274 3182 3277 3299 3264 3312 3223 3264 3271 3269 3255 3281 3375 3213 3237 3331 3255 3267 3248 3258 3142 3163 3254 3175 3224 3237 3279 3321 3329 3156)

f2='test5_1w/EF18131.fasta'
n2=(1536 1527 1508 1561 1515 1466 1458 1493 1492 1526 1486 1534 1526 1464 1563 1526 1470 1450 1566 1511 1408 1542 1481 1446 1565 1466 1530 1511 1497 1553 1503 1488 1547 1544 1522 1498 1500 1535 1459 1474 1493 1496 1469 1539 1528 1454 1477 1514 1471 1472 1451 1555 1539 1476 1554 1523 1444 1478 1553 1552 1510 1515 1443 1447 1504 1482 1539 1445 1485 1485 1496 1474 1443 1575 1470 1492 1505 1506 1445 1463 1508 1430 1470 1537 1540 1498 1546 1526 1494 1448 1547 1499 1473 1515 1469 1551 1550 1455 1493 1564)
f3='test5_1w/EF18132.fasta'
n3=(1620 1602 1598 1608 1645 1610 1610 1580 1687 1605 1582 1697 1579 1561 1578 1580 1575 1684 1639 1664 1618 1595 1613 1564 1619 1640 1660 1667 1592 1595 1684 1578 1604 1584 1562 1593 1526 1589 1607 1586 1660 1585 1618 1620 1611 1621 1615 1618 1570 1666 1653 1602 1582 1558 1576 1559 1633 1653 1633 1655 1616 1543 1610 1651 1630 1566 1591 1678 1600 1598 1640 1581 1594 1570 1616 1586 1625 1575 1602 1643 1611 1623 1658 1595 1577 1633 1611 1604 1637 1674 1640 1628 1603 1663 1593 1591 1595 1570 1640 1652)
f4='test5_1w/EF18133.fasta'
n4=(812 841 836 872 805 786 824 812 809 872 801 850 813 815 827 806 796 827 814 782 822 796 794 905 849 809 781 837 851 799 787 847 805 781 881 869 900 816 822 840 821 815 836 791 728 843 847 823 813 813 884 806 789 846 826 770 833 821 825 800 834 835 867 825 840 844 854 776 820 884 788 860 815 774 839 864 844 873 851 851 824 798 827 844 903 799 736 759 820 802 814 793 795 802 820 818 795 821 748 839)
f5='test5_1w/EF18134.fasta'
n5=(2801 2849 2769 2736 2814 2857 2878 2830 2775 2761 2883 2699 2895 2923 2787 2855 2857 2798 2788 2838 2905 2822 2888 2819 2830 2779 2767 2766 2834 2782 2875 2848 2823 2859 2781 2828 2855 2795 2887 2822 2827 2866 2849 2850 2863 2844 2850 2754 2929 2759 2818 2788 2825 2820 2810 2960 2853 2705 2804 2823 2907 2857 2816 2780 2777 2835 2804 2846 2806 2803 2802 2903 2871 2782 2811 2746 2803 2782 2831 2774 2802 2868 2670 2811 2743 2739 2852 2844 2801 2818 2857 2917 2875 2845 2894 2803 2781 2833 2790 2789)

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
    #    count=count+10
    #fi
done

