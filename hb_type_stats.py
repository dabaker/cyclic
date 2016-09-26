#!/usr/bin/python
from  collections import Counter
from cyclic_utilities import *

#############################

if len(argv) <=1:
    print '\n'
    print 'USAGE: %s <silent file> <pdb_list> <distance_cutoff>  '%argv[0]
    print '\n output counts of turns classified by ABXY string and hb pattern \n'
    print ' if <silent file> = NONE get torsions from idealized pdbs \n'
    print ' to suppress individual pdb annotation use "-no_output" \n'
    print '\n\n'
    exit()

args = argv[1:]
output =1
if args.count('-no_output'):
   output=0

input_file=args[0]
pdb_list=args[1]
cutoff=float(args[2])
cutoff2=cutoff*cutoff

contents={}
list_file=open(pdb_list,'r').readlines()
res_column=5
pathology={}
hbond_list={}
sc_hbond_list={}

abba=[]
tor_bin_list=[]
for line in list_file:
    pdb = string.split(line)[0]
#    print map(string.split,popen('grep CA %s | wc'%pdb).readlines())
    nres=int(map(string.split,popen('grep CA %s | wc'%pdb).readlines())[0][0])
    if nres > 150: 
        print nres, pdb, 'OVER 200 residues \n'
        continue
    contacts, N_hbonds, max_hb_to_O, sc_O_contacts, sc_N_contacts=compute_contacts(pdb,res_column,cutoff2)
    value='pass'
    if max_hb_to_O > 2:
	value='fail'
    pdb1=string.split(pdb,'/')[-1]
    if '_0001' in pdb1:
        index= pdb1.find('_0001')
        pdb_tag= pdb1[0:index]
    else:
        index= pdb1.find('.pdb')
        pdb_tag= pdb1[0:index]

    pathology[pdb_tag]=value
    hbond_list[pdb_tag]=contacts
    sc_hbond_list[pdb_tag]=sc_O_contacts + sc_N_contacts
    if input_file=='NONE':   #pdbs are idealized, get torsions from them
        torsions=get_torsions_from_ideal_pdb(pdb)
        res=torsions.keys()
        max_res=max(res)
        bin_str=''
        for n in range(1,max_res+1):
            if (n) in torsions.keys():
                r=torsions[n]
            else:
                r=' '
            bin_str=bin_str+r
        
#        print bin_str
        abba.append((bin_str,bin_str,bin_str,bin_str,bin_str,contacts,value,pdb_tag))
        contents[pdb_tag]=[]
        tor_bin_list.append(bin_str)
if input_file != 'NONE':
 decoy_num,info,names,seqs,score=input_silent(input_file)
 #tor_bin_list=[]
 #input torsion bin strings and sequences
 for decoy in info.keys():
    bin_str=""
    for res in info[decoy]:
        bin_str=bin_str+info[decoy][res][0]
    a_bin_str,offset=alphabetize_ABBA(bin_str)
    a_bin_str_inv,offset=alphabetize_ABBA(invert_ABBA(bin_str))
    a_bin_str = min(a_bin_str,a_bin_str_inv)
    tor_bin_list.append(bin_str)
    index= names[decoy].find('_0001')
    pdb_tag= names[decoy][0:index]
    if pdb_tag in pathology.keys():
        abba.append( (bin_str,a_bin_str,offset,parse_seq(seqs[decoy]),score[decoy],hbond_list[pdb_tag],pathology[pdb_tag],pdb_tag) )
    else:
        print 'Missing PDB: %s'%pdb_tag
        continue
    contents[pdb_tag]=[]
nres=len(bin_str)

hb_3=[]
hb_4=[]
hb_2=[]
hb_4c=[]
# three cases: i-i+3 only, i-i+3 AND (i-1)-i+3, i-i+3 AND (i+1)-(i+4)
t1=[]
t2=[]
t3=[]
t4=[]   # if none of above, consider i,i+4 hbonds
t5=[]  #case 2 AND 3
t6=[]  # i,i+4 and i+1,i+5
t7=[] # i,i+4 and i,i+5
t0=[] # all 2 residue turns

for entry in abba:
    exclude=[]
    if entry[6]=='fail': continue
    hbonds=entry[5]
    bin_str=entry[0]
    pdb_tag=entry[7]
    for hbond in hbonds:
        if hbond in exclude: continue   # only count hbond in one configuration
        res1=hbond[0]-1
        res2=hbond[1]-1
    
        if ( (res2-res1) == 3 or (res1-res2) == (nres-3)):
# have i, i+3 hbond from CO of res1 to NH of res2
# count all cases
            if res2 > res1:
                    t0.append(bin_str[res1+1:res2])
            else:
                    t0.append(bin_str[(res1+1):nres] + bin_str[0:res2])  
# check for case 2 and 3 together
            if ( ((res1+1)%nres)+1,((res2+1)%nres)+1 ) in hbonds and  (  res1 + 1, (res2+1)%nres+1) in hbonds:
                if res2 > res1:
                    t5.append(bin_str[res1+1:res2+1])
                else:
                    t5.append(bin_str[(res1+1):nres] + bin_str[0:res2+1])
                exclude.append(((res1+1)%nres+1,((res2+1)%nres)+1 ))
                exclude.append ((  res1 + 1, (res2+1)%nres+1))
                contents[pdb_tag].append( (bin_str,res1+1,res2+2,hbonds,'i,i+3 and i+1,i+4 and i,i+4') )
            elif (  res1+ 1, (res2+1)%nres+1) in hbonds:
#case 2
                exclude.append( ( ( (res1)%nres) + 1, (res2+1)%nres+1))
                contents[pdb_tag].append( (bin_str,res1+1,res2+2,hbonds,'i,i+3 and i,i+4') )
                if res2 > res1:
                    t2.append(bin_str[ res1+1:res2+1])
                else:
                    t2.append(bin_str[(res1+1):nres]+bin_str[0:res2+1])
            elif ( ((res1+1)%nres)+1,((res2+1)%nres)+1 ) in hbonds:
#case 3
                contents[pdb_tag].append( (bin_str,res1+1,res2+2,hbonds,'i,i+3 and i+1,i+4') )
                exclude.append( ((res1+1)%nres+1,((res2+1)%nres)+1) )
                if res2 > res1:
                    t3.append(bin_str[res1+1:res2+1])
                else:
                    t3.append(bin_str[(res1+1):nres] + bin_str[0:res2+1])
            else:
#case 1
                contents[pdb_tag].append( (bin_str,res1+1,res2+1,hbonds,'i,i+3') )
                if res2 > res1:
                    t1.append(bin_str[res1+1:res2])
                else:
                    t1.append(bin_str[(res1+1):nres] + bin_str[0:res2])           
        else:
# 3 resideue turn
             if ( (res2-res1) == 4 or (res1-res2) == (nres-4)):
                if res2==0: 
                    r2=nres
                else:
                    r2=res2
                if ( (res1+1), r2) in hbonds: continue

                if ( (res1+1)%nres +1, (res2+1)%nres +1) in hbonds:
                  exclude.append(( (res1+1)%nres +1, (res2+1)%nres +1))
                  contents[pdb_tag].append( (bin_str,res1+1,res2+2,hbonds,'i,i+4 and i+1,i+5') )
                  if res2 > res1:
                    t6.append(bin_str[res1+1:res2+1])
                  else:
                    t6.append(bin_str[(res1+1):nres] + bin_str[0:res2+1])

                elif ( (res1+1),(res2+1)%nres+1) in hbonds:
                    contents[pdb_tag].append( (bin_str,res1+1,res2+2,hbonds,'i,i+4 and i,i+5'))
                    exclude.append( ( (res1+1),(res2+1)%nres+1))
                    if res2 > res1:
                        t7.append(bin_str[res1+1:res2+1])
                    else:
                        t7.append(bin_str[res1+1:nres]+bin_str[0:res2+1])
                        
                elif res2 > res1:
                    contents[pdb_tag].append( (bin_str,res1+1,res2+1,hbonds,'i,i+4') )
                  #  if bin_str[res1+1:res2]=='BXX': print 'BXX ',bin_str,hbonds,entry[7]
                    t4.append(bin_str[res1+1:res2])

                else:
                    contents[pdb_tag].append( (bin_str,res1+1,res2+1,hbonds,'i,i+4') )
                 #   if bin_str[res1+1:nres]+bin_str[0:res2] == 'BXX':  print 'wBXX ',bin_str,hbonds,entry[7]
                  
                    t4.append(bin_str[res1+1:nres]+bin_str[0:res2])

t0_ct= Counter(t0).most_common()
t2_ct= Counter(t2).most_common()
t3_ct= Counter(t3).most_common()
t4_ct= Counter(t4).most_common()
t1_ct= Counter(t1).most_common()
t5_ct= Counter(t5).most_common()
t6_ct= Counter(t6).most_common()
t7_ct= Counter(t7).most_common()
print '# total counts of substrings in two central residues of i, i+3 hbond'
for entry in t0_ct:
    if len(entry[0])< 6:    print entry[0],entry[1]
print
print '#counts of substrings in two central residues of i, i+3 hbond excluding cases below'
for entry in t1_ct:
 if len(entry[0])< 6:    print entry[0],entry[1]
print
print '#counts of substrings in three central residues of i, i+4 hbond excludign cases below'
for entry in t4_ct:
 if len(entry[0])< 6:    print entry[0],entry[1]
print
print '#counts of substrings in i,i+3 AND i+1,i+4 hbonds'
for entry in t3_ct:
 if len(entry[0])< 6:    print entry[0],entry[1]
print
print '#counts of substrings in i,i+3 AND i,i+4 hbonds'
for entry in t2_ct:
 if len(entry[0])< 6:    print entry[0],entry[1]
print
print '#counts of substrings in i,i+3 AND i,i+4 AND i+1,i+4 hbonds'
for entry in t5_ct:
 if len(entry[0])< 6:    print entry[0],entry[1]
print
print '#counts of substrings in i,i+4 AND i +1, i+5 hbonds'
for entry in t6_ct:
 if len(entry[0])< 6:    print entry[0],entry[1]
print
print '#counts of substrings in i,i+4 AND i,i+5 hbonds'
for entry in t7_ct:
 if len(entry[0])< 6:    print entry[0],entry[1]
print
print '************ '
print 'pdb turn types'
if output:
 for pdb in contents.keys():
  print
  for hbond in contents[pdb]:
    if hbond[2] > hbond[1]: 
        str=hbond[0][hbond[1]:hbond[2]-1]
    else:
        str=hbond[0][hbond[1]:nres]+ hbond[0][0:hbond[2]-1]
    print pdb, hbond[0],hbond[1],hbond[2],hbond[3],hbond[4],str
exit()
#calculate frequencies of torsion bins
bin_freq,ax,by,tot=get_tor_bin_freqs(tor_bin_list)

for a in bin_freq.keys():
    bin_freq[a]=bin_freq[a]/float(tot)
#bin_freq['A']=bin_freq['A']/float(tot)
#bin_freq['B']=bin_freq['B']/float(tot)
#bin_freq['X']=bin_freq['X']/float(tot)
#bin_freq['Y']=bin_freq['Y']/float(tot)
print 
print
print 'bin_frequencies',bin_freq,' symmetrized: ',ax,by
# write out torsion substrings by frequency
for subL in range(2,5):
 print ' '
 print '****  substring: ',subL,' *****'
 sub_strings=substring_counter(tor_bin_list,nres,subL,0)
 sorted_subs=sorted(sub_strings, key=sub_strings.__getitem__)
 sorted_subs.reverse()
 for sub in  sorted_subs:
    frac=tot
    for t in list(sub):
	frac=frac*bin_freq[t]
#    if sub < invert_ABBA(sub):
#    print sub,sub_strings[sub]
    print '%s %s   %.2f  %.2f %.3f '%(sub,sub_strings[sub],frac,tot,(float(sub_strings[sub]))/frac)



#hb_3_ct=collections.Counter(hb_3)
#print hb_3_ct
#hb_4_ct=collections.Counter(hb_4)

#for entry in hb_3:

 #   print entry, hb_3.count(entry)

#for entry in hb_4:
 #   print entry, hb_4.count(entry)
