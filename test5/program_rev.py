#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 9 12:36:11 2018


"""

import sys
import numpy
import re
import math
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type = str, required = True,
    help = "input fasta file")
    parser.add_argument("-n", "--num", type = int, required = True,
    help = "number of reads")
    args = parser.parse_args()
    return vars(args)

def SwitchStr(c):
    group = {
        'A':'T',
        'T':'A',
        'C':'G',
        'G':'C'
    }
    return group.get(c,None)


def WriteinTxt(output, sample, rand, seqnum, x):
    output += "@M00302:21:000000000-ANHLP:1:1101:" + str(rand) +":" + str(seqnum) + " " + str(x) + ":N:0:1"+"\n"
    output += sample + "\n"
    output += "+\n"
    output += "#8ACCGGGGGGGGGGGG8FAFGGGGFGFGC#:CFGGGGGGGGFGGGGGGGGGGGFFFGGGGCFGGGGG8EFGDGG\n"
    return output


def GetBeginningSequence(output, subtar, snum, rand, seqnum):
    sample = subtar[0:snum]
    output = WriteinTxt(output,sample,rand,seqnum,1)
    return output

def GetND(x1,x2,e,d):
    result = -1.0;
    try:
        result = e + math.sqrt(d*d)*math.sqrt((-2)*(math.log(x1)/math.log(math.e)))*math.cos(2*math.pi*x2)

    except ValueError:
        result = -1
    finally:
        return result


def GetEndSequence(output, subtar, snum, rand, seqnum):
    temp = []
    for i in subtar:
        temp.append(i)
    temp.reverse()
    subtemp = temp[0:snum]

    for i in range(0, snum):
        subtemp[i]=SwitchStr(subtemp[i])

    sample = ''.join(subtemp)

    output=WriteinTxt(output,sample,rand,seqnum,2)
    return output


def GetSeq(output,rand,snum,leng, index):
    first = rand+1
    last = rand+leng-snum
    output += "group"+str(index)+"   left: "+str(first)+"             right: "+str(last) +"\n"
    return output



try:
    c_args = parse_args()
    #filename = input("Please enter the file name: ")
    filename = c_args["input"]
    file = open(filename)
    ls = []
    for line in file:
        if not line.startswith('>'):
            ls.append(line.replace('\n',''))
    #print(ls[0])
    file.close()
    tar = ls[0]
    strlen = len(tar)
    print("The String length is " + str(strlen))

    #read group
    #numstr = input("The number of reads do you want: ")
    numstr = c_args["num"]
        #while not numstr.isdigit():
    #numstr = input("Wrong input, please input again: ")
    gnumber=int(numstr)


    #efiled = input("Please input mean: ")
    #while not efiled.isdigit():
    #    efiled = input("Wrong input, please input again: ")
    e=float(500)

    #read standard deviation
    #dfiled = input("Please input standard deviation(not variance): ")
    #while (not dfiled.isdigit()) or (dfiled == "0"):
    #    dfiled = input("Wrong input, please input again: ")
    d = float(50)

    #print("Generating " + str(gnumber) + " groups of random number for Normal distribution")
    shortnum = 0
    lengroup = numpy.random.normal(e,d,gnumber)

    for i in range(0, gnumber):
        temp = int(round(lengroup[i]))
        lengroup[i]=temp

        if i==0:
            shortnum = temp
        if temp < shortnum:
            shortnum = temp


    #print(lengroup)
    print("Generated "+ str(gnumber) + " groups of random number!")




    #read the number of samples in every group
    #snumstr = input("The number of DNA be read each end" + str(shortnum)+":")
    #while (not snumstr.isdigit()) or (int(snumstr)*2>shortnum):
    #    snumstr = input("Wrong input, please input again: ")
    snumber = int(75)
   # for j in range(0,10):
        #加循环注意空格 python！！！！
        #f1 = open('left0'+str(j)+'.txt','w')
        #f2 = open('right0'+str(j)+'.txt','w')
       # f3 = open('location0'+str(j)+'.txt','w')
    f1 = open(filename+str(numstr)+'_left.txt','w')
    f2 = open(filename+str(numstr)+'_right.txt','w')
    f3 = open(filename+str(numstr)+'_location.txt','w')


    seqnum = 1000
    #example: if length = 1500, in every group we get 500 characters, the lastnum is 1001

    output1 = ''
    output2 = ''
    output3 = ''
    for i in range(0,gnumber):
            lastnum = strlen - int(lengroup[i])+1
            seqnum +=1
        #if length = 1500, index is from 0 to 1499, and the randomNum is from 0 to 1000(lastnum - 1)
            rand1 = numpy.random.randint(0, lastnum)
            rand2 = numpy.random.randint(10000, 100000)
            subtar = tar[rand1:rand1+int(lengroup[i])]


            output1=GetBeginningSequence(output1, subtar, snumber, rand2, seqnum)
            output2=GetEndSequence(output2, subtar, snumber, rand2, seqnum)
            output3=GetSeq(output3,rand1,snumber,int(lengroup[i]),i+1)

    f1.write(output1)
    f2.write(output2)
    f3.write(output3)
    f1.close()
    f2.close()
    f3.close()
    print("Job is done!")


except OSError as err:
    print("Error:{0}", format(err))

except ValueError:
    print("Error Input")
