#!/usr/bin/python
from cyclic_utilities import *
from sys import argv
from os import system, popen
import string

if len(argv) <=1:
    print '\n'
    print 'USAGE: %s <silent file> -L <substring length> -output_seqbin -symmetrize \n'%argv[0]
    print 'reports frequencies of torsion bin strings up to <substring length>.   default substring length is 5 \n'
    print 'to output torsion bin strings corresponding to D/L/p/P sequence string use -output_seqbin \n'
    print 'to symmetrize bin strings (only lower alphabetical order isomer reported) use -symmetrize \n'
    print 'requires standard (non - binary) format silent file \n'
    print '\n'
    exit()

args=argv[1:]
infile=args[0]
substring_length=5
if args.count('-L'):
    pos=args.index('-L')
    substring_length=int(args[pos+1])
    print 'new substring length: ',substring_length


if args.count('-output_seqbin'):
    output_seqbin=1
else:
    output_seqbin=0

if args.count('-symmetrize'):
    symmetrize=1
else:
    symmetrize=0

decoy_num,info,name,seqs,score=input_silent(infile)
tor_bin_list=[]
ori_tor_bin_list=[]
#input torsion bin strings and sequences
for decoy in info.keys():
    bin_str=""
    for res in info[decoy]:
        bin_str=bin_str+info[decoy][res][0]
    ori_tor_bin_list.append(bin_str)
    a_bin_str,offset=alphabetize_ABBA(bin_str)
    a_bin_str_inv,offset=alphabetize_ABBA(invert_ABBA(bin_str))
    a_bin_str = min(a_bin_str,a_bin_str_inv)
    tor_bin_list.append(a_bin_str)

if output_seqbin:
 #convert sequences to seq_bin_strings
 seq_torstr={}
 alph_binstr_data=[]
 alph_binstr=[]
 for decoy in info.keys():
    seq = parse_seq(seqs[decoy])
    tor_bin_str=tor_bin_list[decoy-1]
    bin_str=''   
    for a in seq:
        bin_str=bin_str+get_seqbin(a)
    a_binstr,offset=alphabetize_ABBA(bin_str)    
    alph_binstr_data.append( (a_binstr,tor_bin_str,seq ))
    alph_binstr.append(a_binstr)
 # track torsion bin strings seen for each bin string 
    if a_binstr in seq_torstr.keys():
	if tor_bin_str not in seq_torstr[a_binstr]:
		seq_torstr[a_binstr].append(tor_bin_str)
    else:
	seq_torstr[a_binstr]=[]
	seq_torstr[a_binstr].append(tor_bin_str)

 #write out torsion bin strings observed for each sequence bin string
 alph_binstr_data.sort()
 for entry in alph_binstr_data:
    print 'SEQBIN',entry[0],alph_binstr.count(entry[0]),entry[1],string.join(entry[2]),seq_torstr[entry[0]]

#calculate frequencies of torsion bins
bin_freq_ori,ax,by,oz,tot=get_tor_bin_freqs(ori_tor_bin_list)
bin_freq,ax,by,oz,tot=get_tor_bin_freqs(tor_bin_list)
print ' ************* \n'
print ' total counts all bins: ',tot
print 'total counts of torsion bins before and after string alphabetization: '
for key in bin_freq.keys():
    print key, bin_freq_ori[key],bin_freq[key]

#output per position frequencies
get_tor_bin_freqs_per_position(ori_tor_bin_list)
#print 'bin_freq after alphabetize',bin_freq,' symmetrized: ',ax,by,oz

if symmetrize:
 bin_list=tor_bin_list
 invert=1
 bin_freq['A']=ax
 bin_freq['B']=by
 bin_freq['X']=ax
 bin_freq['Y']=by
 bin_freq['O']=oz
 bin_freq['Z']=oz

else:
 bin_list=ori_tor_bin_list
 invert=0
 for key in bin_freq_ori.keys():
     bin_freq[key]=bin_freq_ori[key]/float(tot)
#     print bin_freq
nres=len(bin_str)
# write out torsion substrings by frequency
for subL in range(1,substring_length):
 print ' '
 print '****  substring: ',subL,' *****'
 sub_strings=substring_counter(bin_list,nres,subL,invert)
 sorted_subs=sorted(sub_strings, key=sub_strings.__getitem__)
 for sub in  sorted_subs:
    frac=float(tot)
    for t in list(sub):
 #       print frac
	frac=frac*bin_freq[t]
  #      print t, frac
    if not symmetrize or sub < invert_ABBA(sub):
      print '%s %s %.2f    %smers'%(sub,sub_strings[sub],(float(sub_strings[sub]))/frac,subL)
     
