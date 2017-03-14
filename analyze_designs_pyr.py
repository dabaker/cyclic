#!/usr/bin/python
from cyclic_utilities import *
from pdb_utils_noclass import *
from  collections import Counter
#############################

if len(argv) <=1:
    print '\n'
    print 'USAGE: %s <silent file>  '%argv[0]
    print '\n output torsion bin strings and hb maps of non pathological designs '
    print '\n to set number of structures written to cluster files use num_to_output \n'
    print 'designs with hbond pathology and aro-pro motifs flagged \n'
    print 'SubGroup.dat file lists lowest energy member of each group'
    exit()
args = argv[1:]
input_file=args[0]
num_to_output=20
if args.count('num_to_output'):
    pos=args.index('num_to_output')
    num_to_output=int(args[pos+1])
    print 'number of structures written to cluster files: ',num_to_output
init_pyrosetta()
t,h,score_list_ori,seq_list_ori=input_silent_score_seq(input_file)
#print seq_list_ori
tor_string_list={}
hbond_list={}
score_list={}
seq_list={}
pathology={}
pdbs = t.keys()
for pdb in pdbs:
    if '_0001' in pdb:
        index= pdb.find('_0001')
        pdb_tag= pdb[0:index]
    else:
        index= pdb.find('.pdb')
        pdb_tag= pdb[0:index]
    score_list[pdb_tag]=score_list_ori[pdb]

    seq=seq_list_ori[pdb]
    seq_list[pdb_tag]=parse_seq_pose(seq)
    hb_to_O=[]
    contacts=[]
    for hbond in h[pdb]:
        contacts.append( (hbond[0],hbond[1]))
        hb_to_O.append(hbond[0])

    count=Counter(hb_to_O)
#    print count.most_common()
    try:
        if count.most_common()[0][1] > 2:
            value='fail'
            print 'pathology: ',count.most_common(),pdb
        else:
            value='pass'
    except:
        value='pass'   
    pathology[pdb_tag]=value

    hbond_list[pdb_tag]=contacts
    nres=len(t[pdb])
    bin_str=''
    for tor in t[pdb]:
            bin_str=bin_str+get_abba(tor[0],tor[1],tor[2])
    tor_string_list[pdb_tag]=bin_str
#    abba.append((bin_str,contacts,value,pdb_tag))
#    tor_bin_list.append(bin_str)

#print pathology.keys()
#print hbond_list.keys()


contact_count={}
abba_hb_info={}
abba=[]
abba_strings={}
for pdb_tag in hbond_list.keys():
    if pdb_tag in pathology.keys():
        if pathology[pdb_tag]=="fail": continue
    else:
        print 'tag not in pathology.keys: ',pdb_tag
    score=score_list[pdb_tag]
    bin_str=tor_string_list[pdb_tag]
    a_bin_str,offset=alphabetize_ABBA(bin_str)
#    a_bin_str_inv,offset=alphabetize_ABBA(invert_ABBA(bin_str))
#    a_bin_str = min(a_bin_str,a_bin_str_inv)
#    tor_bin_list.append(a_bin_str)
#NEED TO GET SCORE AND SEQ FROM SILENT
    num = len(bin_str)
    offset_hbonds=[]
    for hbond in hbond_list[pdb_tag]:
       if (hbond[1]-hbond[0])%num > 2: 
         offset_hbonds.append( ( (hbond[0] -1 -offset)%num,(hbond[1]-1 -offset)%num))
    offset_hbonds.sort()
    hbonds_str='  bb:'
    for hb in offset_hbonds:
        bond='-'.join(map(str,hb))
        hbonds_str=hbonds_str+' '+bond
    seq_l=list(seq_list[pdb_tag])
    offset_seq=[]    
    for i in range(num):
        offset_seq.append( (seq_l[(i+offset)%num]) )
    seq_str=string.join(offset_seq)

    if (a_bin_str < alphabetize_ABBA(invert_ABBA(bin_str))[0]):
     if (a_bin_str,hbonds_str) in abba_hb_info.keys():
        abba_hb_info[ (a_bin_str,hbonds_str) ].append((score,seq_str+" "+pdb_tag+hbonds_str,offset_hbonds))
     else:
        abba_hb_info[ (a_bin_str,hbonds_str)]=[]
        abba_hb_info[ (a_bin_str,hbonds_str) ].append((score,seq_str+" "+pdb_tag+hbonds_str,offset_hbonds))

    abba.append( (bin_str,pdb_tag,a_bin_str,offset,seq_str,hbonds_str,score_list[pdb_tag]))
    if a_bin_str in abba_strings.keys():
        abba_strings[a_bin_str]=abba_strings[a_bin_str]+1
    else:
        abba_strings[a_bin_str]=1

#now do same for mirror seq/struct
    bin_str_inv=invert_ABBA(bin_str)
    a_bin_str_inv,offset=alphabetize_ABBA(invert_ABBA(bin_str))
    seq_l_inv=mirror_seq(seq_l)
    offset_seq=[]
    for i in range(num):
        offset_seq.append( (seq_l_inv[(i+offset)%num]) )
    seq_str=string.join(offset_seq)

    offset_hbonds=[]
    for hbond in hbond_list[pdb_tag]:
         if (hbond[1]-hbond[0])%num > 2:
             offset_hbonds.append( ( (hbond[0] -1 -offset)%num,(hbond[1]-1 -offset)%num))
    offset_hbonds.sort()
    hbonds_str='  bb:'
    for hb in offset_hbonds:
        bond='-'.join(map(str,hb))
        hbonds_str=hbonds_str+' '+bond
    if  a_bin_str_inv < a_bin_str:
     if (a_bin_str_inv,hbonds_str) in abba_hb_info.keys():
        abba_hb_info[ (a_bin_str_inv,hbonds_str) ].append((score,seq_str+" "+'%s_inv'%pdb_tag+hbonds_str,offset_hbonds))
     else:
        abba_hb_info[ (a_bin_str_inv,hbonds_str)]=[]
        abba_hb_info[ (a_bin_str_inv,hbonds_str) ].append((score,seq_str+" "+'%s_inv'%pdb_tag+hbonds_str,offset_hbonds))
    abba.append( (bin_str_inv,pdb_tag+"_inv",a_bin_str_inv,offset,seq_str,hbonds_str,score_list[pdb_tag]) )
    if a_bin_str_inv in abba_strings.keys():
        abba_strings[a_bin_str_inv]=abba_strings[a_bin_str_inv]+1
    else:
        abba_strings[a_bin_str_inv]=1

# torsion_strings={}
# torsion_strings_inv={}
# seqs={}
# scores={}
# offset_list_inv={}
# offset_list={}
# for entry in abba:
# #   tot_pro=entry[0].count('P') + entry[0].count('Q')
# #   print "ABBAPQ", entry[3]," ",entry[0]," ", entry[2]," ", tot_pro
# #   print "PRO ", tot_pro,res_file_cutoff,entry[3]
# #   if tot_pro >= res_file_cutoff:
# #      output_resfile(entry[3],entry[0])
 
#    index= entry[3].find('_0001')
#    pdb_tag= entry[3][0:index]
# #   print entry[3],entry[0],entry[2],pdb_tag 
#    if 'inv' in entry[3]:
#  #       print "inverse"
#         torsion_strings_inv[pdb_tag]=(entry[0],entry[2])
#         offset_list_inv[pdb_tag]=entry[4]
#    else:
#         torsion_strings[pdb_tag]=(entry[0],entry[2])
#         offset_list[pdb_tag]=entry[4]
#         seqs[pdb_tag]=entry[5]
#         scores[pdb_tag]=entry[6]
# #print 'pdb   pathology  seq   score bb_hbonds  tor_bins  alphabet_tor_bins  inv_tor_bins  alphabet_inv_ter_bins'
# #for pdb_tag in torsion_strings.keys():

#  #  print pdb_tag,pathology[pdb_tag],seqs[pdb_tag],scores[pdb_tag],hbond_list[pdb_tag],torsion_strings[pdb_tag][0],torsion_strings[pdb_tag][1],torsion_strings_inv[pdb_tag][0],torsion_strings_inv[pdb_tag][1]
# abba_strings={}
# abba_names={}
# cluster_strings={}
# cluster_size={}
# hb_in_abba={}
# seq_in_abba={}
# sc_in_abba={}
# abba_hb_info={}
# for entry in abba:
# #    clust=int(string.split(entry[3],'.')[1])
#     score=entry[6]
#     index= entry[3].find('_0001')
#     inv=0
#     if "inv" in entry[3]:
#         inv=1
#     pdb_tag= entry[3][0:index]
#     if pdb_tag not in pathology.keys(): 
#         print 'Skipping missing pdb: ',pdb_tag
#         continue
#     if pathology[pdb_tag]=='fail': 
#         print 'Removed 3 hb pathology: ', pdb_tag
#         continue

#     sort_str=entry[2]
#     if (sort_str > alphabetize_ABBA(invert_ABBA(sort_str))[0]): 
# #        print 'gt',sort_str,  alphabetize_ABBA(invert_ABBA(sort_str))[0]
#         continue
# #    print 'lt',sort_str,  alphabetize_ABBA(invert_ABBA(sort_str))[0]  
#     offset_seq=[]
#     offset=int(entry[4])
#     seq_list=parse_seq(seqs[pdb_tag])
    
#     seq_l=seq_list
# #    print seq_l,seqs[pdb_tag],'parse'

#     if inv:
#         seq_l=mirror_seq(seq_list)

#     if len(seq_l) < num:
#         print 'WARNING: seq_l,seqs[pdb_tag]'
#     for i in range(num):
#         offset_seq.append( (seq_l[(i+offset)%num]) )
    
#     seq_str=string.join(offset_seq)
#   #  print inv,entry[0],sort_str,offset,seqs[pdb_tag],seq_l,seq_str,'seq'
#     if flag_aro_pro(seq_str) == 0: 
#         print 'Warning aro-pro dipeptide: ',entry

#     if inv:  
#         clust=pdb_tag+'_inv'
#     else:
#         clust=pdb_tag
 
#     offset_hbonds=[]
#     sc_offset_hbonds=[]
    
#     for sc_hbond in sc_hbond_list[pdb_tag]:
# #        print sc_hbond
#         offset0=(sc_hbond[0]-1 -offset)%num
#         offset1=(sc_hbond[1]-1 -offset)%num
# #        if string(sc_hbond[0])[0]=='s':
# #            offset0='s%s'%(int(sc_hbond[0][1:])-offset)%num
# #        else:
# #             offset0=(sc_hbond[0]-offset)%num
# #       if sc_hbond[1][0]=='s':
# #            offset1='s%s'%(int(sc_hbond[1][1:])-offset)%num
        
#         sc_offset_hbonds.append( ( offset0,offset1))
#     sc_offset_hbonds.sort()
#     sc_hbonds_str='  sc:'
#     for hb in sc_offset_hbonds:
#         bond='-'.join(map(str,hb))
#         sc_hbonds_str=sc_hbonds_str+' '+bond

#     for hbond in hbond_list[pdb_tag]:
#         offset_hbonds.append( ( (hbond[0] -1 -offset)%num,(hbond[1]-1 -offset)%num))
#   #  if pdb_tag=='c.11.96': print offset,hbond_list[pdb_tag],offset_hbonds,num, bin_str
#     offset_hbonds.sort()
#     hbonds_str='  bb:'
#     for hb in offset_hbonds:
#         bond='-'.join(map(str,hb))
#         hbonds_str=hbonds_str+' '+bond
   
#     if clust in cluster_strings.keys():
#        if sort_str not in cluster_strings[clust]:
#           cluster_strings[clust].append( sort_str +' '  +  hbonds_str )
#        cluster_size[clust]=cluster_size[clust]+1
#     else:
#        cluster_strings[clust]=[]
#        cluster_strings[clust].append( ( sort_str + ' ' + hbonds_str)  )
#        cluster_size[clust]=1

#     if sort_str in abba_strings.keys():
#         abba_strings[sort_str]=abba_strings[sort_str]+1
# #        index= entry[3].find('_0001')
# #        pdb_tag= entry[3][0:index]
#         abba_names[sort_str].append( clust + hbonds_str)
#         if hbonds_str in hb_in_abba[sort_str].keys():
#             hb_in_abba[sort_str][hbonds_str]=hb_in_abba[sort_str][hbonds_str]+1
#         else:
#             hb_in_abba[sort_str][hbonds_str]=1
#         seq_in_abba[sort_str].append(seq_str + "  " + clust + hbonds_str +sc_hbonds_str)
#         if (sort_str,hbonds_str) in abba_hb_info.keys():
#             abba_hb_info[(sort_str,hbonds_str)].append((score,seq_str + "  " + clust + hbonds_str +sc_hbonds_str,offset_hbonds))
#         else:
#             abba_hb_info[(sort_str,hbonds_str)]=[]
#             abba_hb_info[(sort_str,hbonds_str)].append((score,seq_str + "  " + clust + hbonds_str +sc_hbonds_str,offset_hbonds))
#     else:
#         abba_strings[sort_str]=1
#         hb_in_abba[sort_str]={}
#         hb_in_abba[sort_str][hbonds_str]=1
#         seq_in_abba[sort_str]=[]
#         abba_hb_info[(sort_str,hbonds_str)]=[]
#         abba_hb_info[(sort_str,hbonds_str)].append( (score, seq_str + "  " + clust + hbonds_str +sc_hbonds_str,offset_hbonds))
#         seq_in_abba[sort_str].append(seq_str + "  " +  clust + hbonds_str + sc_hbonds_str)
# #        index= entry[3].find('_0001')
# #        pdb_tag= entry[3][0:index]
#         abba_names[sort_str]=[]
#         abba_names[sort_str].append(  clust + hbonds_str)
        
# hbs_abba={}
# for sort_str in abba_strings.keys():
#     hbs=[]
#     for hb in hb_in_abba[sort_str].keys():
#         hbs.append( (hb_in_abba[sort_str][hb],hb))
#     hbs.sort()
#     hbs.reverse()
#     hb_str=''
#     for hb in hbs:
#         hb_str=hb_str + str(hb[0])+':'+hb[1] + ' /  '
#     hbs_abba[sort_str]=hb_str

# str_list=[]
# for sort_str in abba_strings.keys():
# #    PQ=sort_str.count('P') + sort_str.count('Q')
#     str_list.append( (abba_strings[sort_str],sort_str,' | '.join(abba_names[sort_str] ),hbs_abba[sort_str],alphabetize_ABBA(invert_ABBA(sort_str))[0] ))



# str_list.sort()

# for entry in str_list:
#     print 'cluster members: ', entry[1]," ",entry[4]," ", entry[0],"  ",entry[2]

# for entry in str_list:
#     print 'hb patterns: ',entry[1],"  ",entry[4]," ",entry[0],"  ",entry[3]

# # to reduce number of groups, for each ABBA string, assign groups with hbond patters
# # completely contained within a second group to that second group


inf=abba_hb_info.keys()
inf.sort()
prev_ABBA=''

hb_list_in_ABBA={}
## ****not currently using next two blocks which are for collapsing hb_patterns ***
for entry in inf:
    bin_str=entry[0]
    hb=entry[1]
    if bin_str != prev_ABBA:
        if prev_ABBA != '': hb_list_in_ABBA[prev_ABBA]=hb_patterns
        hb_patterns=[]
        prev_ABBA=bin_str
#    print bin_str,hb,abba_hb_info[(bin_str,hb)]
    hb_patterns.append( (abba_hb_info[(bin_str,hb)][0][2],hb) )

for bin_str in hb_list_in_ABBA.keys():
    for hb_l1 in hb_list_in_ABBA[bin_str]:
        hb_list1=hb_l1[0]
        hb1=hb_l1[1]
        for hb_l2 in hb_list_in_ABBA[bin_str]:
            hb_list2=hb_l2[0]
            hb2=hb_l2[1]
            if hb_list1 != hb_list2:
                if set(hb_list1) < set(hb_list2):
                    hb_list_in_ABBA[(bin_str,hb1)]=hb_list2
#                    print bin_str,hb_list1,hb_list2
                if set(hb_list2) < set(hb_list1):
                    hb_list_in_ABBA[(bin_str,hb2)]=hb_list1
 #                   print bin_str,hb_list1,hb_list2
SubGroup_list=[]
for entry in inf:
    if entry[0] != prev_ABBA:
        print
        print '*****Group: ', entry[0],'****'
        prev_ABBA=entry[0]
        if abba_strings[entry[0]]>6:
            output=1
            out_file=open('%s_%s'%(abba_strings[entry[0]],entry[0]),'w')
            out_file.write('%s \n'%entry[0])
            
        else:
            output=0

    sort_hb_info=abba_hb_info[entry]
    sort_hb_info.sort()
    hb_list=sort_hb_info[0][2]
    bin_str=entry[0]
    score=sort_hb_info[0][0]
# following two lines collapse hb_lists
# disable for now as sometimes fewer hbonds give better energies
#    if (bin_str,entry[1]) in hb_list_in_ABBA.keys():
#        hb_list=hb_list_in_ABBA[(bin_str,entry[1])]
    turns= find_turn_types(bin_str,hb_list)    
    turn_str=''
    for turn in turns:
        turn_str=turn_str+'%s%s '%(turn[2],turn[3])
    print 'SubGroup size: %s '%len(sort_hb_info),' ',bin_str,' ',sort_hb_info[0][1],' turns: ',turn_str,' score: %s'%score
    SubGroup_list.append( (score,bin_str,sort_hb_info,turn_str) )
#    print 'Group size: %s '%len(sort_hb_info),' ',entry[0],' ',sort_hb_info[0][1],hb_list,find_turn_types(entry[0],sort_hb_info[0][2])

    counter=0
    for x in sort_hb_info:
        print '%s %s energy: %7.3f '%(entry[0],x[1],x[0])
        if output:
            counter=counter+1
            if counter <= num_to_output: 
                 out_file.write('%s %7.3f \n'%(x[1],x[0]))

SubGroup_list.sort()
outf=open('SubGroups.dat','w')
for entry in SubGroup_list:
    if entry[0] < -3.:
        outf.write('SubGroup size: %5d %s %25s  turns: %s energy: %7.3f \n'%(len(entry[2]),entry[1],entry[2][0][1],entry[3],entry[0]))



# for sort_str in seq_in_abba.keys():
# #    print sort_str
    
#     for seq in seq_in_abba[sort_str]:
#         mismatch= seq_bin_mismatch(sort_str,seq) 
# #        if mismatch < 0.5: 
# #            print seq
        
# print
# #clus=cluster_strings.keys()
# #clus.sort()
# #for entry in clus:
# #   print entry, cluster_size[entry], cluster_strings[entry]
