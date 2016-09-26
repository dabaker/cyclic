#!/usr/bin/python
from os import system,popen,path
import string
from sys import argv,stderr
from os.path import exists


def get_Pbin(phi,psi,omega,psi_prev):
   Pbin="N"
#   if (psi_prev <  -25 or psi_prev > 40):
#      if (phi < -50 and phi > -60 and psi  > -40 and psi < -20) or (phi < -55 and phi > -75 and psi > 135 and psi < 160): Pbin="P"

#   if (psi_prev <  -40 or psi_prev > 25):
#      if (phi > 50 and phi < 60 and psi  < 40 and psi > 20) or (phi > 55 and phi < 75 and psi < -135 and psi > - 160): Pbin="Q"

   return Pbin

def get_abba(phi,psi,omega):
   split_B_bin=0
   if (omega==None and phi==None and psi==None):
       return "I"
   elif (omega<90.0 and omega>-90.0):
       if (phi < 0.0):
           return "O"
       else:
           return "Z"
   else:
       if phi < 0.0:
           if (psi <50.0 and psi > -80.0):
               return "A"
           else:
             if  split_B_bin:
               if (psi >0):
                  return "B"
               else:
                  return "C"
             else:
                return "B"
       else:
           if (psi <80.0 and psi > -50.0):
               return "X"  	### Aprime
           else:
             if split_B_bin:
              if (psi < 0) :
                 return "Y"
              else:
                 return "W"  ### Bprime
             else:
                return "Y"

def get_tor_bin_freqs(tor_bin_list):
#calculate frequencies of torsion bins
 bin_freq={}
 tot=0
 for tor_bin in tor_bin_list:
	for tor in list(tor_bin):
           tot=tot+1
           if tor in bin_freq.keys():
		bin_freq[tor]=bin_freq[tor]+1
           else:
                bin_freq[tor]=1

 ax=(bin_freq['A']+bin_freq['X'])/float(tot)
 by=(bin_freq['B']+bin_freq['Y'])/float(tot)
 oz=(bin_freq['O']+bin_freq['Z'])/float(tot)
# print 'bin_freq',bin_freq,' symmetrized: ',ax,by
 return bin_freq,ax,by,oz,tot
 
def get_tor_bin_freqs_per_position(tor_bin_list):
   bin_count=[]
   nres=len(tor_bin_list[0])
   ndecoy=len(tor_bin_list)
   for i in range(nres):
      bin_count.append({})

   for i in range(ndecoy):
      bin_str=list(tor_bin_list[i])
   
      for j in range(nres):
       r=bin_str[j]
       if r in bin_count[j].keys():
         bin_count[j][r]=bin_count[j][r]+1
       else:
         bin_count[j][r]=1

   print '########## per position counts  '
   print 'position  A     X      B     Y  '
   for i in range(nres):
      print '%4d %6d %6d %6d %6d'%(i+1,bin_count[i]['A'],bin_count[i]['X'],bin_count[i]['B'],bin_count[i]['Y'])



def get_seqbin(seq):
  aa_1=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','H']
  d_aa_1=['a','c','d','e','f','g','h','i','k','l','m','n','p','q','r','s','t','v','w','y','h']
  if seq in aa_1:
     bin='L'
  if seq in d_aa_1:
     bin='D'
  if seq == 'p':
     bin='p'
  if seq == 'P':
     bin='P'
  return bin

def substring_counter(strings,length,n_sub,invert):
   sub_strings={}
   for string in strings:
      if len(string) != length:
         print 'ERROR-wrong length string'
      for i in range(length):
         sub=''
         for j in range(n_sub):
            sub=sub+string[(i+j)%length]
         if sub in sub_strings.keys():
            sub_strings[sub]=sub_strings[sub]+1
         else:
            sub_strings[sub]=1
	 if invert:
          inv_sub=invert_ABBA(sub)
          if inv_sub in sub_strings.keys():
            sub_strings[inv_sub]=sub_strings[inv_sub]+1
          else:
            sub_strings[inv_sub]=1
         
   return sub_strings


def invert_ABBA(instring):
   bin_list=list(instring)
   outstring=""
   for bin in bin_list:
      if bin == 'A': outstring=outstring+'X'
      if bin == 'B': outstring=outstring+'Y'
      if bin == 'X': outstring=outstring+'A'
      if bin == 'Y': outstring=outstring+'B'
      if bin == 'O': outstring=outstring+'Z'
      if bin == 'Z': outstring=outstring+'O'
      if bin == 'P': outstring=outstring+'Q'
      if bin == 'Q': outstring=outstring+'P'
      if bin == 'C': outstring=outstring+'W'
      if bin == 'W': outstring=outstring+'C'
#   print instring, outstring
   return outstring

def alphabetize_ABBA(bin_str):
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
    return permute[0],offset

def output_resfile(name,res_string):
    #print name, res_string
    resfile=open('%s.resfile'%name,'w')
    resfile.write('ALLAAxc \nEX 1 EX 2 \nUSE_INPUT_SC \n\nstart\n\n')
    reslist=list(res_string)
    for i in range(len(reslist)):

       if reslist[i]=='P' :
          resfile.write('%s A PIKAA P\n'%(i+1))
       elif reslist[i]=='Q':
              resfile.write('%s A EMPTY NC DPR\n'%(i+1))
       elif reslist[i]=='X' or reslist[i]=='Y':
              resfile.write('%s A EMPTY NC DAL NC DAS NC DGU NC DPH NC DHI NC DIL NC DLY NC DLE NC DME NC DAN NC DPR NC DGN NC DAR NC DSE NC DTH NC DVA NC DTR NC DTY\n'%(i+1))
       else:
          resfile.write("\n")
    resfile.close()
    return

def mirror_seq(input_seq):
   d_aa=['DALA','DCYS','DASP','DGLU','DPHE','DGLY','DHIS','DILE','DLYS','DLEU','DMET','DASN','DPRO','DGLN','DARG','DSER','DTHR','DVAL','DTRP','DTYR','DHIS_D','DHIS_']

   aa=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','HIS_D','HIS_']
   aa_1=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','H','H']
   d_aa_1=['a','c','d','e','f','g','h','i','k','l','m','n','p','q','r','s','t','v','w','y','h','h']
   output_seq=[]
   for s in input_seq:
      
      if s in d_aa_1:
         out=aa_1[d_aa_1.index(s)]
      else:
         out=d_aa_1[aa.index(s)]
      output_seq.append(out)
   return output_seq

def input_silent(file):
    info = {}
    bb={}
    seq={}
    name={}
    score={}
    decoy_num=0
    silent_file =  open(file,'r').readlines()
    for line in silent_file :
        l = string.split(line)
    #    print l[0],' ',l[1]
#        if l[0]=='SEQUENCE:' :
#             sequence=l[1]
#             for i in range(len(sequence)):
#                    seq[i+1]=sequence[i]



        if l[0]=='SCORE:' and l[1] != 'score':
	    scorel=float(l[1])
            decoy_num=decoy_num+1
            psi_prev=[-100]  # for first residue, don't have to worry about pre-PRO restrictions
            if decoy_num >= 1:
    #            print decoy_num,bb
 #               info[decoy_num]= bb

                bb={}

        if l[0]=='ANNOTATED_SEQUENCE:' :
            decoy_name=l[2]
            
            if decoy_num >= 1:
                name[decoy_num]=decoy_name
                seq[decoy_num]=l[1]
                score[decoy_num]=scorel
        if l[1] in ['E','H','L']  :
            tors=map(float, l[2:5] )
            tbin=get_abba(tors[0],tors[1],tors[2])
            Pbin=get_Pbin(tors[0],tors[1],tors[2],psi_prev)
            if Pbin=="N":
               bin=tbin
            else:
               bin=Pbin
            bb[int(l[0])] = (bin,tors)
            info[decoy_num]=bb   #this is really dumb,but works (only need to update after last residue)
            psi_prev=tors[1]
        
 #           name[decoy_num]=decoy_name
 #   info[decoy_num]=bb
 #   name[decoy_num]=decoy_name
 #   print 'decoy_num', decoy_num
#    print seq
    return decoy_num,info,name,seq,score

def seq_bin_mismatch(sort_str,seq):
    d_aa_1=['a','c','d','e','f','g','h','i','k','l','m','n','p','q','r','s','t','v','w','y','h']
    aa_1=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','H']
    num=len(sort_str)
    mismatch=0
    for i in range(num):
       if sort_str[i]=='A' or sort_str[i]=='B':
          if seq[(2*i)] not in aa_1: 
             mismatch=mismatch+1
    #         print i,seq,sort_str,seq[i],sort_str[i]
       if sort_str[i]=='X' or sort_str[i]=='Y':
          if seq[(2*i)] not in d_aa_1: mismatch=mismatch+1
    return mismatch
             
def flag_aro_pro(s):
   aro=['W','w','Y','y','F','f']
   pro=['P','p']
   n=len(s)
   a_p=1
   for i in range(n):
      j=(i+1)%n
      if (s[i] in aro and s[j] in pro) or (s[j] in aro and s[i] in pro):
         a_p=0
   return a_p

def parse_seq(s):
    d_aa=['DALA','DCYS','DASP','DGLU','DPHE','DGLY','DHIS','DILE','DLYS','DLEU','DMET','DASN','DPRO','DGLN','DARG','DSER','DTHR','DVAL','DTRP','DTYR','DHIS_D']
    aa=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','HIS_D']
    d_aa_1=['a','c','d','e','f','g','h','i','k','l','m','n','p','q','r','s','t','v','w','y','h']
    aa_1=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','H']
    seq_n=[]
    d=0
    res=0
    for c in range(len(s)):
#       print c,s[c],res,seq_n
       if s[c] == '[':
          d=1
          str=seq_n[res-1]=s[c+1:c+5]
          if str in d_aa:
             str=d_aa_1[d_aa.index(str)]
          if str=='HIS_':
             str='H'
          seq_n[res-1]=str         
       elif  s[c]== ']':
          d=0
       else: 
         if d == 0:
            str=s[c]
            if s[c] in aa:
               str=aa_1[aa.index(s[c])]
            seq_n.append(str)
            res=res+1
    return seq_n   

def get_torsions_from_ideal_pdb(pdb):
   torsions={}
   
   lines=map(string.split,popen('grep REMARK %s | grep -v torsions'%pdb).readlines())
   for line in lines:
      torsions[int(line[2])]=(get_abba(float(line[6]),float(line[7]),float(line[8])))
   return torsions

def compute_contacts(pdb,res_column,cutoff2):
# note changed from "N" to "H"; list names don't reflect this
    min_sep=2   ## don't count contacts between adjacent residues
    max_hb_to_O=0
    contacts= []
    N_hbonds=[]
    coords_N = map(string.split,popen('grep " H   "  %s | grep ATOM'%pdb).readlines())
    coords_O = map(string.split,popen('grep " O   " %s | grep ATOM'%pdb).readlines())



    coords_sc= map(string.split,popen('grep -v " H   " %s  | grep -v "HB" | grep -v " HA "| grep -v " O   " | grep -v " N  " | grep -v "C" | grep ATOM'%pdb).readlines())

    nres = len(coords_O)

    for i in range(len(coords_O)):
        N_contacts=0
        for j  in range(len(coords_N)):
          
          res1=int(coords_O[i][res_column])
          try:
             res2=int(coords_N[j][res_column])
          except:
             print 'PDB read error:',pdb,res2
             exit()
          if res1 != res2:
            xyz1=map(float,[coords_O[i][res_column +1],coords_O[i][res_column+2],coords_O[i][res_column+3]])
            xyz2=map(float,[coords_N[j][res_column+1],coords_N[j][res_column+2],coords_N[j][res_column+3]])
#            print coords[i][5]

            dist=0
            for k in range(3):
                     dist=dist+(xyz1[k]-xyz2[k])**2
#            print res1,res2,dist
            if dist < cutoff2:
                if abs(res1-res2) > min_sep and abs(res1-res2) < (nres-min_sep): 
                 N_contacts=N_contacts+1
                 contacts.append( (res1,res2) )
        N_hbonds.append( (res1,N_contacts) )
        if N_contacts> max_hb_to_O: max_hb_to_O = N_contacts

    sc_O_contacts=[]
    sc_N_contacts=[]
    for i in range(len(coords_O)):
        for j in range(len(coords_sc)):

                xyz1=map(float,[coords_O[i][res_column +1],coords_O[i][res_column+2],coords_O[i][res_column+3]])
                xyz2=map(float,[coords_sc[j][res_column+1],coords_sc[j][res_column+2],coords_sc[j][res_column+3]])
                res1=int(coords_O[i][res_column])
                res2=int(coords_sc[j][res_column])
                if res1 != res2:
                    dist=0
                    for k in range(3):
                        dist=dist+(xyz1[k]-xyz2[k])**2
                    if dist < cutoff2:

                        sc_O_contacts.append( (res2,res1) )
  #                      print 'sc hb to O',coords_O[i],coords_sc[j]
    for i in range(len(coords_N)):
        for j in range(len(coords_sc)):

                xyz1=map(float,[coords_N[i][res_column +1],coords_N[i][res_column+2],coords_N[i][res_column+3]])
                xyz2=map(float,[coords_sc[j][res_column+1],coords_sc[j][res_column+2],coords_sc[j][res_column+3]])
                res1=int(coords_N[i][res_column])
                res2=int(coords_sc[j][res_column])

                if res1 != res2 and coords_sc[j][11] !='H':
                    dist=0
                    for k in range(3):
                        dist=dist+(xyz1[k]-xyz2[k])**2
                    if dist < cutoff2:

                        sc_N_contacts.append( (res2,res1) )
#                        print 'sc to N', coords_N[i],coords_sc[j]
              ## print res1,res2,dist
    return contacts , N_hbonds, max_hb_to_O, sc_O_contacts, sc_N_contacts

def find_turn_types(bin_str,hbonds):
  names=['i:i+3','i:i+3 | i:i+4','i:i+3 | i+1:i+4','i:i+3 | i+1:i+4 | i:i+4','i:i+4 ','i:i+4 | i:i+5','i:i+4 | i+1:i+5']
  names=['','a','b','ab','','a','b']
  exclude=[]
  turn_list=[]
  nres=len(bin_str)
  for hbond in hbonds:
        if hbond in exclude: continue   # only count hbond in one configuration
        res1=hbond[0]
        res2=hbond[1]
        if ( (res2-res1) == 3 or (res1-res2) == (nres-3)):
#  i, i+3 hbond from CO of res1 to NH of res2

# check for case 2 and 3 together
            if ( ((res1+1)%nres),((res2+1)%nres) ) in hbonds and  (  res1 , (res2+1)%nres) in hbonds:
                if res2 > res1:
                    str=bin_str[res1+1:res2+1]
                else:
                    str=bin_str[(res1+1):nres] + bin_str[0:res2+1]
                exclude.append(((res1+1)%nres+1,((res2+1)%nres)+1 ))
                exclude.append ((  res1 + 1, (res2+1)%nres+1))
                turn_list.append( (res1,res2,str,names[3]) )

            elif (  res1, (res2+1)%nres) in hbonds:
#case 2
                exclude.append( ( ( (res1)%nres) , (res2+1)%nres))

                if res2 > res1:
                    str=bin_str[ res1+1:res2+1]
                else:
                    str=bin_str[(res1+1):nres]+bin_str[0:res2+1]
                turn_list.append( (res1,res2,str,names[1] ) )

            elif ( ((res1+1)%nres),((res2+1)%nres) ) in hbonds:
#case 3
                exclude.append( ((res1+1)%nres,((res2+1)%nres)) )
                if res2 > res1:
                    str=bin_str[res1+1:res2+1]
                else:
                    str=bin_str[(res1+1):nres] + bin_str[0:res2+1]
                turn_list.append( (res1,res2,str,names[2] ) )
            else:
#case 1

                if res2 > res1:
                    str=bin_str[res1+1:res2]
                else:
                    str=bin_str[(res1+1):nres] + bin_str[0:res2]
                turn_list.append( (res1,res2,str,names[0]) )
        else:
# 3 resideue turn
             if ( (res2-res1) == 4 or (res1-res2) == (nres-4)):
                if res2==0: 
                    r2=nres
                else:
                    r2=res2
                if ( (res1+1), r2) in hbonds: continue

                if ( (res1+1)%nres , (res2+1)%nres ) in hbonds:
                  exclude.append(( (res1+1)%nres , (res2+1)%nres ))
                  if res2 > res1:
                    str=bin_str[res1+1:res2+1]
                  else:
                    str=bin_str[(res1+1):nres] + bin_str[0:res2+1]
                  turn_list.append( (res1,res2,str,names[5]))
                elif ( (res1),(res2+1)%nres) in hbonds:
                    exclude.append( ( (res1),(res2+1)%nres))
                    if res2 > res1:
                        str=bin_str[res1+1:res2+1]
                    else:
                        str=bin_str[res1+1:nres]+bin_str[0:res2+1]
                    turn_list.append( (res1,res2,str,names[6]) )    
                elif res2 > res1:
                    str=bin_str[res1+1:res2]
                    turn_list.append( (res1,res2,str,names[4]) )
                else:
                    str=bin_str[res1+1:nres]+bin_str[0:res2]
                    turn_list.append( (res1,res2,str,names[4]))

  return turn_list 
