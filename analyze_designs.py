#!/usr/bin/python
from cyclic_utilities import *

#############################

if len(argv) <=1:
    print '\n'
    print 'USAGE: %s <silent file> <pdb_list> <distance_cutoff>  '%argv[0]
    print '\n output torsion bin strings and hb maps of non pathological designs '
    print '\n\n'

args = argv[1:]
input_file=args[0]
pdb_list=args[1]
cutoff=float(args[2])
cutoff2=cutoff*cutoff



list_file=open(pdb_list,'r').readlines()
contact_count={}
res_column=5
pathology={}
hbond_list={}
sc_hbond_list={}
for line in list_file:
    pdb = string.split(line)[0]
    contacts, N_hbonds, max_hb_to_O, sc_O_contacts, sc_N_contacts=compute_contacts(pdb,res_column,cutoff2)
    value='pass'
    if max_hb_to_O > 2:
	value='fail'
    pdb1=string.split(pdb,'/')[-1]
    index= pdb1.find('_0001')
    pdb_tag= pdb1[0:index]
#    pdb_tag=string.split(string.split(pdb,'/')[-1],'.')[0]
#    print pdb_tag
    pathology[pdb_tag]=value
    hbond_list[pdb_tag]=contacts
    sc_hbond_list[pdb_tag]=sc_O_contacts + sc_N_contacts
#    print pdb_tag, max_hb_to_O,value
#    print 'backbone hbonds: ', pdb, contacts
#    print 'sc to O hbonds: ',sc_O_contacts
#    print 'sc to N hbonds: ',sc_N_contacts



decoy_num,info,names,seq=input_silent(input_file)
#print info
abba=[]
for decoy in info.keys():
    bin_str=""
   # print names[decoy]
    for res in info[decoy]:
#        print info[decoy][res]
        bin_str=bin_str+info[decoy][res][0]
    permute=[]
    num=len(bin_str)
#    print "BIN",bin_str,num
    if num>0:
        for aa in range(num):
            sort_str=""
            for i in range(num):
                j=(i+aa)%num
#                print "j",i+aa,num,j
                sort_str=sort_str+bin_str[j:j+1]
            permute.append(sort_str)
        ori=[]
        for i in range(len(permute)):
            ori.append(permute[i]) 
        permute.sort()
        offset=ori.index(permute[0])
 #       print offset, bin_str,permute[0],ori,permute
 #       print permute
        abba.append( (bin_str,decoy,permute[0],names[decoy],offset,seq[decoy]) )
 
ndecoy=len(abba)
for i in range(ndecoy):
    entry=abba[i]
    bin_str=invert_ABBA(entry[0])
    permute=[]
    num=len(bin_str)
#    print "BIN",bin_str,num
    if num>0:
        for aa in range(num):
            sort_str=""
            for i in range(num):
                j=(i+aa)%num
#                print "j",i+aa,num,j
                sort_str=sort_str+bin_str[j:j+1]
            permute.append(sort_str)
        ori=[]
        for i in range(len(permute)):
            ori.append(permute[i])  
        permute.sort()
        offset=ori.index(permute[0])

  #      print offset, bin_str,permute[0]
 #       print permute
        abba.append( (bin_str,entry[1],permute[0],entry[3]+"_inv",offset,entry[5]+"_inv") )

torsion_strings={}
torsion_strings_inv={}
seqs={}

offset_list_inv={}
offset_list={}
for entry in abba:
#   tot_pro=entry[0].count('P') + entry[0].count('Q')
#   print "ABBAPQ", entry[3]," ",entry[0]," ", entry[2]," ", tot_pro
#   print "PRO ", tot_pro,res_file_cutoff,entry[3]
#   if tot_pro >= res_file_cutoff:
#      output_resfile(entry[3],entry[0])
 
   index= entry[3].find('_0001')
   pdb_tag= entry[3][0:index]
#   print entry[3],entry[0],entry[2],pdb_tag 
   if 'inv' in entry[3]:
 #       print "inverse"
        torsion_strings_inv[pdb_tag]=(entry[0],entry[2])
        offset_list_inv[pdb_tag]=entry[4]
   else:
        torsion_strings[pdb_tag]=(entry[0],entry[2])
        offset_list[pdb_tag]=entry[4]
        seqs[pdb_tag]=entry[5]
print 'pdb   pathology  seq   bb_hbonds  tor_bins  alphabet_tor_bins  inv_tor_bins  alphabet_inv_ter_bins'
for pdb_tag in torsion_strings.keys():
    
    print pdb_tag,pathology[pdb_tag],seqs[pdb_tag],hbond_list[pdb_tag],torsion_strings[pdb_tag][0],torsion_strings[pdb_tag][1],torsion_strings_inv[pdb_tag][0],torsion_strings_inv[pdb_tag][1]
abba_strings={}
abba_names={}
cluster_strings={}
cluster_size={}
hb_in_abba={}
seq_in_abba={}
sc_in_abba={}
abba_hb_info={}
for entry in abba:
#    clust=int(string.split(entry[3],'.')[1])
    index= entry[3].find('_0001')
    inv=0
    if "inv" in entry[3]:
        inv=1
    pdb_tag= entry[3][0:index]
    if pathology[pdb_tag]=='fail': continue
    sort_str=entry[2]
    if (sort_str > alphabetize_ABBA(invert_ABBA(sort_str))[0]): 
#        print 'gt',sort_str,  alphabetize_ABBA(invert_ABBA(sort_str))[0]
        continue
#    print 'lt',sort_str,  alphabetize_ABBA(invert_ABBA(sort_str))[0]  
    offset_seq=[]
    offset=int(entry[4])
    seq_list=parse_seq(seqs[pdb_tag])
    
    seq_l=seq_list
#    print seq_l,seqs[pdb_tag],'parse'

    if inv:
        seq_l=mirror_seq(seq_list)

    if len(seq_l) < num:
        print 'WARNING: seq_l,seqs[pdb_tag]'
    for i in range(num):
        offset_seq.append( (seq_l[(i+offset)%num]) )
    
    seq_str=string.join(offset_seq)
    print inv,entry[0],sort_str,offset,seqs[pdb_tag],seq_l,seq_str,'seq'
        

    if inv:  
        clust=pdb_tag+'_inv'
    else:
        clust=pdb_tag
 
    offset_hbonds=[]
    sc_offset_hbonds=[]
    
    for sc_hbond in sc_hbond_list[pdb_tag]:
#        print sc_hbond
        offset0=(sc_hbond[0]-1 -offset)%num
        offset1=(sc_hbond[1]-1 -offset)%num
#        if string(sc_hbond[0])[0]=='s':
#            offset0='s%s'%(int(sc_hbond[0][1:])-offset)%num
#        else:
#             offset0=(sc_hbond[0]-offset)%num
#       if sc_hbond[1][0]=='s':
#            offset1='s%s'%(int(sc_hbond[1][1:])-offset)%num
        
        sc_offset_hbonds.append( ( offset0,offset1))
    sc_offset_hbonds.sort()
    sc_hbonds_str='  sc:'
    for hb in sc_offset_hbonds:
        bond='-'.join(map(str,hb))
        sc_hbonds_str=sc_hbonds_str+' '+bond

    for hbond in hbond_list[pdb_tag]:
        offset_hbonds.append( ( (hbond[0] -1 -offset)%num,(hbond[1]-1 -offset)%num))
    offset_hbonds.sort()
    hbonds_str='  bb:'
    for hb in offset_hbonds:
        bond='-'.join(map(str,hb))
        hbonds_str=hbonds_str+' '+bond
   
    if clust in cluster_strings.keys():
       if sort_str not in cluster_strings[clust]:
          cluster_strings[clust].append( sort_str +' '  +  hbonds_str )
       cluster_size[clust]=cluster_size[clust]+1
    else:
       cluster_strings[clust]=[]
       cluster_strings[clust].append( ( sort_str + ' ' + hbonds_str)  )
       cluster_size[clust]=1

    if sort_str in abba_strings.keys():
        abba_strings[sort_str]=abba_strings[sort_str]+1
#        index= entry[3].find('_0001')
#        pdb_tag= entry[3][0:index]
        abba_names[sort_str].append( clust + hbonds_str)
        if hbonds_str in hb_in_abba[sort_str].keys():
            hb_in_abba[sort_str][hbonds_str]=hb_in_abba[sort_str][hbonds_str]+1
        else:
            hb_in_abba[sort_str][hbonds_str]=1
        seq_in_abba[sort_str].append(seq_str + "  " + clust + hbonds_str +sc_hbonds_str)
        if (sort_str,hbonds_str) in abba_hb_info.keys():
            abba_hb_info[(sort_str,hbonds_str)].append(seq_str + "  " + clust + hbonds_str +sc_hbonds_str)
        else:
            abba_hb_info[(sort_str,hbonds_str)]=[]
            abba_hb_info[(sort_str,hbonds_str)].append(seq_str + "  " + clust + hbonds_str +sc_hbonds_str)
    else:
        abba_strings[sort_str]=1
        hb_in_abba[sort_str]={}
        hb_in_abba[sort_str][hbonds_str]=1
        seq_in_abba[sort_str]=[]
        abba_hb_info[(sort_str,hbonds_str)]=[]
        abba_hb_info[(sort_str,hbonds_str)].append(seq_str + "  " + clust + hbonds_str +sc_hbonds_str)
        seq_in_abba[sort_str].append(seq_str + "  " +  clust + hbonds_str + sc_hbonds_str)
#        index= entry[3].find('_0001')
#        pdb_tag= entry[3][0:index]
        abba_names[sort_str]=[]
        abba_names[sort_str].append(  clust + hbonds_str)
        
hbs_abba={}
for sort_str in abba_strings.keys():
    hbs=[]
    for hb in hb_in_abba[sort_str].keys():
        hbs.append( (hb_in_abba[sort_str][hb],hb))
    hbs.sort()
    hbs.reverse()
    hb_str=''
    for hb in hbs:
        hb_str=hb_str + str(hb[0])+':'+hb[1] + ' /  '
    hbs_abba[sort_str]=hb_str

str_list=[]
for sort_str in abba_strings.keys():
#    PQ=sort_str.count('P') + sort_str.count('Q')
    str_list.append( (abba_strings[sort_str],sort_str,' | '.join(abba_names[sort_str] ),hbs_abba[sort_str],alphabetize_ABBA(invert_ABBA(sort_str))[0] ))



str_list.sort()

for entry in str_list:
    print 'cluster members: ', entry[1]," ",entry[4]," ", entry[0],"  ",entry[2]

for entry in str_list:
    print 'hb patterns: ',entry[1],"  ",entry[4]," ",entry[0],"  ",entry[3]

inf=abba_hb_info.keys()
inf.sort()
prev_ABBA=''
for entry in inf:
    if entry[0] != prev_ABBA:
        print
        print entry[0]
        prev_ABBA=entry[0]
        if abba_strings[entry[0]]>6:
            output=1
            out_file=open('%s_%s'%(abba_strings[entry[0]],entry[0]),'w')
            out_file.write('%s \n'%entry[0])
            
        else:
            output=0
    for x in abba_hb_info[entry]:
        print x
        if output:
            out_file.write('%s \n'%x)

for sort_str in seq_in_abba.keys():
#    print sort_str
    
    for seq in seq_in_abba[sort_str]:
        mismatch= seq_bin_mismatch(sort_str,seq) 
#        if mismatch < 0.5: 
#            print seq
        
print
#clus=cluster_strings.keys()
#clus.sort()
#for entry in clus:
#   print entry, cluster_size[entry], cluster_strings[entry]
