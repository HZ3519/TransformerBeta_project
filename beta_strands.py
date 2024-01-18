import os
import gzip
import sys
from traceback import print_exc #to manage exceptions
#from itertools import product
#import pickle
import pickle as pickle
import numpy
import scipy.stats
from collections import OrderedDict
import subprocess
import threading



#import csv_dict
#import misc

try :
    import mybio
except ImportError :
    sys.stderr.write("WARNING mybio.py module not found - excecution may fail\n")
    class mybio :
        a=1
        def composition(self,*args,**kwargs):
            raise Exception("required mybio.py MODULE NOT AVAILABLE\n")
try :
    import structs
except ImportError :
    class structs :
        a=1
        def Stats(self,*args,**kwargs):
            raise Exception("required structs.py MODULE NOT AVAILABLE\n")

RUN_CAMSOL_PUBLISHED_VERSION=False
try :
    from camsol_intrinsic_linear import camsol_intrinsic
    HAVE_camsol=True
except Exception :
    HAVE_camsol=False
    sys.stderr.write("WARNING zyggregator not found!!!\n")
    def camsol_intrinsic(*args,**kwargs):
        return

# check ids map for !





amino_list1=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

"""
################  DSSP DATABASE CREATION AND UPTADE ################
basically is all down to run the function:
Have a backup of the existing db before moving on..

import beta_strands
dssp_database_folder_path='DSSP'
chain_list_file_name=beta_strands.update_dssp_database(dssp_database_folder_path, update_repr_chain_list=True, chain_list_identity=90,chain_list_name=None, overwrite_dssp=False)


provided that the rsync service is working (a good reason to leave overwrite_dssp=False)
it also download a chain list from the rcsb database with the required degree of identity.
it then converts it into a representative chains list (the format that can be used afterwards) convert_to_representative_chain_list=True
it does nothing if the folder has been updated in the last 10 days (only_if_older_than_days=10)

CAREFUL!!!! The above procedure sometimes duplicate files if some are .gz and some are not... Double check afterwards!!  <<<<<<------------------******

################ CREATION OF THE BS and BSn DATABASES FROM THE DSSP FILES ################

python beta_strands.py --all DSSP pdb90_chain_list.list

# WHERE the list file pdb90_chain_list.list must contain the keyword 'list' in its filename or extension.
#       DSSP is the name of the folder with the DSSP files




##  LOOK INTO STRAND DICT TO LOCATE A source PDB ID 
## DISUSED by function get_pdb_from_complementary_strands()

import cPickle as pickle
from beta_strands import Beta_strands
strand_dict= pickle.load(open('/Users/francesco/Dropbox/CAM/Grafting_project/DB_DSSP/strand_dict.pkl','rb'))
strange_strand_dict= pickle.load(open('/Users/francesco/Dropbox/CAM/Grafting_project/DB_DSSP/strange_strand_dict.pkl','rb'))

bit_on_protein='TAV' # this must be N to C oriented
bit_complementary='QFT' # this must be oriented so that first is complementary to first of bit_on_protein and so on, thus for antiparalle it is C to N

for key in strand_dict :
   if bit_on_protein in key :
      for j,comp in enumerate(strand_dict[key].complementary) : # comp will be like (('VQVTVYP', '-'), 1)
         if bit_complementary in comp[0][0] and comp[0][0].find(bit_complementary)==key.find(bit_on_protein):
             print key, comp, strand_dict[key].pdb_ids[comp[0]], strand_dict[key].hb_patterns[j]

for key in strange_strand_dict :
   if bit_on_protein in key :
      for j,comp in enumerate(strange_strand_dict[key].complementary) : # comp will be like (('VQVTVYP', '-'), 1)
         if bit_complementary in comp[0][0] and comp[0][0].find(bit_complementary)==key.find(bit_on_protein) :
             print key, comp, strange_strand_dict[key].pdb_ids[comp[0]], strange_strand_dict[key].hb_patterns[j]
"""


def get_pdb_from_complementary_strands(bit_on_protein,bit_complementary,is_anti=None, first_has_hb=None,strand_dict='/Users/francesco/Dropbox/CAM/Grafting_project/DB_DSSP/strand_dict.pkl',strange_strand_dict='/Users/francesco/Dropbox/CAM/Grafting_project/DB_DSSP/strange_strand_dict.pkl'):
    '''
    bit_on_protein must be N to C oriented
    bit_complementary must be oriented so that first residue is complementary to first residue of bit_on_protein and so on (hence N to C for parallel complementary and C to N for antiparallel ones)
    strand_dict,strange_strand_dict= beta_strands.get_pdb_from_complementary_strands('NLYFQ','LEEIK',strand_dict,strange_strand_dict)
    '''
    if type(strand_dict) is str :
        if os.path.isfile(strand_dict) :
            #from beta_strands import Beta_strands # does not work
            print('Extracting strand_dict from pickle, might take a while..')
            strand_dict= pickle.load(open(strand_dict,'rb'))
            print('Extraction done')
        else :
            sys.stderr.write("ERROR cant find strand_dict %s, USING only strange_strand_dict\n" % (strand_dict))
            strand_dict=None
    if type(strange_strand_dict) is str :  
        if os.path.isfile(strange_strand_dict) :  
            print('Extracting strange_strand_dict from pickle, might take a while..')
            strange_strand_dict= pickle.load(open(strange_strand_dict,'rb'))
            print('Extraction done')
        else :
            sys.stderr.write("ERROR cant find strange_strand_dict %s, USING only strand_dict\n" % (strange_strand_dict))
            strange_strand_dict=None
    if strange_strand_dict is None and strand_dict is None :
        raise Exception("\n\n*** get_pdb_from_complementary_strands() can't run without strand_dict and strange_strand_dict\n")
    if strand_dict is not None :
        print('doing strand_dict...')
        for key in strand_dict :
            if bit_on_protein in key :
                for j,comp in enumerate(strand_dict[key].complementary) : # comp will be like (('VQVTVYP', '-'), 1)
                    if bit_complementary in comp[0][0] :
                        if is_anti is not None :
                            if is_anti==True and comp[0][1]=='+' : continue
                            if is_anti==False and comp[0][1]=='-': continue
                        if first_has_hb is not None :
                            if first_has_hb==True and (strand_dict[key].hb_patterns[j][0]==[None,None] or strand_dict[key].hb_patterns[j][0]==(None,None) ): continue
                            if first_has_hb==False and strand_dict[key].hb_patterns[j][0]!=[None,None] and strand_dict[key].hb_patterns[j][0]!=(None,None) : continue
                        print(key, comp, strand_dict[key].pdb_ids[comp[0]], strand_dict[key].hb_patterns[j])
    if strange_strand_dict is not None :
        print('doing strange_strand_dict...')
        for key in strange_strand_dict :
            if bit_on_protein in key :
                for j,comp in enumerate(strange_strand_dict[key].complementary) : # comp will be like (('VQVTVYP', '-'), 1)
                    if bit_complementary in comp[0][0] :
                        if is_anti is not None :
                            if is_anti==True and comp[0][1]=='+' : continue
                            if is_anti==False and comp[0][1]=='-': continue
                        if first_has_hb is not None :
                            if first_has_hb==True and (strange_strand_dict[key].hb_patterns[j][0]==[None,None] or strange_strand_dict[key].hb_patterns[j][0]==(None,None) ): continue
                            if first_has_hb==False and strange_strand_dict[key].hb_patterns[j][0]!=[None,None] and strange_strand_dict[key].hb_patterns[j][0]!=(None,None) : continue
                        print(key, comp, strange_strand_dict[key].pdb_ids[comp[0]], strange_strand_dict[key].hb_patterns[j])
    return strand_dict, strange_strand_dict      



def make_dssp_database_from_pdb_folder(pdb_folder,dssp_folder, check_for_already_present=True, only_these_files=None, pipe=True):
    '''
    requires dssp installed and it is a mess to run, we strongly advice to use update_dssp_database() instead, it will try to download pre-computed dssp files.
    this function tries to make the dssp database starting from a folder of pdb files
    '''
    if pdb_folder[-1]!='/' :  pdb_folder+='/'
    if dssp_folder[-1]!='/' :  dssp_folder+='/'
    os.system("mkdir "+dssp_folder)
    os.system("mkdir "+dssp_folder[:-1]+"_failed")
    if only_these_files is None :
        to_do=os.listdir(pdb_folder)
    else :
        to_do=only_these_files
    sys.stdout.write("%d files found in folder %s\n" % (len(to_do),pdb_folder))
        
    if check_for_already_present :
        already_done=os.listdir(dssp_folder)
        already_done=[ f.split('.')[0].upper() for f in already_done]
        to_do=[f for f in to_do if f[:4].upper() not in already_done ]
        
    to_do=[ pdb for pdb in to_do if pdb[0]!='.' and '.pdb' in pdb ]
    N=len(to_do)  
    sys.stdout.write("%d files to be analyzed\n" % (N))
        
    i=0
    j=0
    problematic=[]
    for  pdb in to_do :
        if pipe :    
            p=subprocess.Popen("dssp -i "+pdb_folder+pdb+" -o "+dssp_folder+pdb[:4]+'.dssp > /dev/null',shell=True, preexec_fn=os.setsid,stdout=subprocess.PIPE,stderr=subprocess.PIPE) # run dssp and kill its stdout
            t = threading.Timer( 60.0, timeout, [p, pdb] ) # allow at most 60 seconds before terminating it (some pdb stuck dssp)
            t.start()
            out=p.communicate()[1]
            del p
            t.cancel()
            del t
            if out!='' :
                problematic+=[pdb]
                sys.stderr.write("Problem %s    with %s\n" % (out,pdb))
                os.system('cp '+pdb_folder+pdb+' '+dssp_folder[:-1]+'_failed/.')
                j+=1
        else :
            os.system("dssp -i "+pdb_folder+pdb+" -o "+dssp_folder+pdb[:4]+'.dssp > /dev/null')
        i+=1
        if i%20==0 :
            sys.stdout.write('%4d (problem %d) -> %5.2lf %% |' % (i,j,100.*i/(1.*N)) )
            sys.stdout.flush()
    sys.stdout.write('\nRUN FINISHED %d files processed, %d with problems\n\t%s\n\n' % (i,j,str(problematic)))
#callback, process is process = Popen( ... )
# then use t = threading.Timer( 10.0, timeout, [process] )
#    t.start()
# process.join()
def timeout( process, tag=None , use_violence=True):
    if process.poll() is None:
        try:
            if use_violence : os.killpg(process.pid,9)
            else :process.kill()
            sys.stderr.write('\n**Error**: process %d taking too long to complete. Terminating...\n   %s\n\n' % (process.pid,str(tag)))
        except:
            pass

def update_dssp_database(dssp_database_folder_path, just_chain_list=False, update_repr_chain_list=True, chain_list_identity=90,chain_list_name=None, overwrite_dssp=False):
    '''
    this uses rsync to update database
    it also download a chain list from the rcsb database with the required degree of identity.
    it then converts it into a representative chains list (the format that can be used afterwards) convert_to_representative_chain_list=True
    it does nothing if the folder has been updated in the last 10 days (only_if_older_than_days=10)
    '''
    if just_chain_list or dssp_database_folder_path is None  :
        _,chain_list_filename, repr_dict = mybio.download_pdb_chain_list(output_filename=chain_list_name, only_if_older_than_days=0,identity_percentage=chain_list_identity, convert_to_representative_chain_list=True)
        return chain_list_filename, repr_dict
    if dssp_database_folder_path[-1]!='/' : dssp_database_folder_path+='/'
    _,chain_list_filename, repr_dict = mybio.download_pdb_chain_list(output_filename=chain_list_name, only_if_older_than_days=10,identity_percentage=chain_list_identity, convert_to_representative_chain_list=True)
    if repr_dict==False : return chain_list_filename, repr_dict # means that it has been modified in the last only_if_older_than_days days
    
    if overwrite_dssp : os.system("rsync -avz --delete rsync://rsync.cmbi.ru.nl/dssp/ "+dssp_database_folder_path)
    else : os.system("rsync -avz rsync://rsync.cmbi.ru.nl/dssp/ "+dssp_database_folder_path)
    print('rsync done')
    
    return chain_list_filename, repr_dict






SECONDARY_STRUCTURES=['H','E','B','G','I','T','S']


def parse_dssp_line(line,advanced=False,get_hb_patterns=False,energy_threshold=-0.35,warn_file=sys.stderr, warn_tag=''):
    '''
    parse a line from the last field (the actual output) of a dssp file.
    recall that resname=='!' is a missing residue, and !* is a chain break
    it returns resnumber,resname,res_id,chain_id,ss_kind,beta_partner_1,beta_partner_2
    if advanced it also returns beta_bridge_label,beta_sheet_label,chirality,solvent_accessibility
    beta_partner_1 is the resnumber of the 'standard' partner in a beta strand
    beta_partner_2 is non 0 only if the strand actually is part of a beta sheet, forming bridges with two chains (one per side)
    ss_kind_legend:
        H = alpha helix
        B = residue in isolated beta-bridge
        E = extended strand, participates in beta ladder
        G = 3-helix (3/10 helix)
        I = 5 helix (pi helix)
        T = hydrogen bonded turn
        S = bend
    '''
    try :
        chain_id=line[11]
        resname=line[13]
        if resname.islower() : resname='C' # when CYS form Sulfur bond they are labelled as a,b,c,... according to the number of sulfur bonds in the protein.
        if warn_file is not None and resname!='!' and resname not in amino_list1 : warn_file.write('   warn %s non standard amino acid %s in line |%s|\n' % (warn_tag,resname,line.replace('   ',' ')))
        res_id=line[5:11].strip() # remove white spaces
        resnumber=int(line[0:5])
        ss_kind=line[16]
        beta_partner_1=int(line[25:29])
        beta_partner_2=int(line[29:33])
        hb_pattern={}
        if advanced:
            chirality=line[22]
            beta_bridge_label=line[23:25]
            beta_sheet_label=line[33]
            solvent_accessibility=int(line[34:38])
        if get_hb_patterns :
            NtoO_p1,NtoO_e1=line[39:50].split(',')
            NtoO_p1=int(NtoO_p1)
            NtoO_e1=float(NtoO_e1)
            NtoO_p2,NtoO_e2=line[61:72].split(',')
            NtoO_p2=int(NtoO_p2)
            NtoO_e2=float(NtoO_e2)
            OtoN_p1,OtoN_e1=line[50:61].split(',')
            OtoN_p1=int(OtoN_p1)
            OtoN_e1=float(OtoN_e1)
            OtoN_p2,OtoN_e2=line[72:83].split(',')
            OtoN_p2=int(OtoN_p2)
            OtoN_e2=float(OtoN_e2)
            if -1<=NtoO_p1+resnumber-beta_partner_1<=1 or NtoO_e1 < energy_threshold :
                hb_pattern[(NtoO_p1+resnumber,'NtoO')]=NtoO_e1
            if -1<=NtoO_p2+resnumber-beta_partner_2<=1 or NtoO_e2 < energy_threshold :
                hb_pattern[(NtoO_p2+resnumber,'NtoO')]=NtoO_e2
            if -1<=OtoN_p1+resnumber-beta_partner_1<=1 or OtoN_e1 < energy_threshold :
                hb_pattern[(OtoN_p1+resnumber,'OtoN')]=OtoN_e1
            if -1<=OtoN_p2+resnumber-beta_partner_2<=1 or OtoN_e2 < energy_threshold :
                hb_pattern[(OtoN_p2+resnumber,'OtoN')]=OtoN_e2       
    except Exception :
        sys.stderr.write("ERROR in parse_dssp_line %s while processing line:\n%s\n" % (warn_tag,line) )
        if warn_file is not None : print_exc(warn_file)
        raise
    if advanced: return resnumber,resname,res_id,chain_id,ss_kind,beta_partner_1,beta_partner_2,beta_bridge_label,beta_sheet_label,chirality,solvent_accessibility,hb_pattern
    return resnumber,resname,res_id,chain_id,ss_kind,beta_partner_1,beta_partner_2,hb_pattern


def simplify_ss(ss, return_different_in_lower_case=False,unkown_symbol=' '):
    '''
    convert observed secondary structure to either alpha helix, beta strand or coil
    '''
    
    if ss=='G' :
        if return_different_in_lower_case : return 'h'
        return 'H'
        #if myres[k][1][j]=='H' :  number_of_match+=1
        #elif myres[k][1][j]=='C' :  number_of_match+=0.5
    elif ss=='B' :
        if return_different_in_lower_case : return 'c'
        return 'C'
        #if myres[k][1][j]=='E' : number_of_match+=0.5
        #if myres[k][1][j]=='C' : number_of_match+=1
    elif ss=='T' :
        if return_different_in_lower_case : return 'c'
        return 'C'
        #if myres[k][1][j]=='C' : number_of_match+=1
    elif ss=='I' :
        if return_different_in_lower_case : return 'h'
        return 'H'
        #if myres[k][1][j]=='H' : number_of_match+=1
    elif ss=='U' :
        return unkown_symbol
    elif ss not in ['H','E','C' ] : 
        print(ss)
        return unkown_symbol
    return ss


def read_dssp_file(filename, simplify=False, chain_to_grab=None, ingnore_gaps=False, unkown_symbol=' '):
    '''
    reads ss info from a dssp file
    unkown_symbol is always space, add an this if statement if you want to change it..
    note that gaps are reported but always of length 1 (and are reported only if dssp inserts a ! there, which is not always the case).
    '''
    if '.gz' in filename : f= gzip.open(filename,'r')
    else : f=open(filename,'r')
    wholefile= f.read().splitlines() # save everything in wholefile
    f.close()
    chain_to_avoid=[] # we need this as the dssp format is wrong for pdbs with more than 10000 (10^4) residues. We need to exclude those chains with this label
    resnumber,_,_,chain_id,_,_,_,_ = parse_dssp_line(wholefile[-1], warn_file=sys.stderr,get_hb_patterns=True, warn_tag=filename) # process last line
    i=0
    if resnumber > 9999 :
        while resnumber>9999 :
            wholefile.pop()
            i+=1
            resnumber,_,_,chain_id,_,_,_,_ = parse_dssp_line(wholefile[-1], warn_file=sys.stderr,get_hb_patterns=True, warn_tag=filename)
        chain_to_avoid+=[chain_id] # we exclude the chain between the numbering 10000 NOTE THAT EXTRA CHAIN BETA STRANDS COULD MAKE A MESS IN THIS CASE (hopefully the noise they introduce is negligible)
        sys.stderr.write('     file '+filename+' excluding %d lines from the end as resnumber exceeded 10000, also excluding chain %s\n' % (i,chain_id) )
   
    read=False
    ss_kind_str={} # dict of chain ids, values are string with ss results
    sequence={} # dict of chain ids, values is the sequence.
    oldresnumber=None
    # read the file, save the sequence and the strand map
    for line in wholefile : #we read for specific position,as split length can vary in some lines
        if read and line!='' :
            resnumber,resname,_,chain_id,ss_kind,_,_,_ = parse_dssp_line(line, warn_file=None,get_hb_patterns=False, warn_tag=filename)
            if resname=='!' : continue # TER end of chain
            if chain_to_avoid is not None and chain_id in chain_to_avoid : continue
            if chain_to_grab is not None and chain_id !=chain_to_grab : continue
            if chain_id not in ss_kind_str : 
                oldresnumber=None
                ss_kind_str[chain_id]=''
                sequence[chain_id]=''
            if not ingnore_gaps and oldresnumber is not None and oldresnumber!=resnumber-1 :
                if oldresnumber>resnumber : print('Error in read_dssp_file for file %s chain %s oldresnumber>resnumber %d %d' % (filename,chain_id,oldresnumber,resnumber))
                else : 
                    sequence[chain_id]+= '-'*(resnumber-oldresnumber-1)
                    ss_kind_str[chain_id]+=unkown_symbol*(resnumber-oldresnumber-1)
            sequence[chain_id]+=resname
            if simplify : ss_kind_str[chain_id] += simplify_ss(ss_kind, return_different_in_lower_case=True)
            else : ss_kind_str[chain_id] += ss_kind
            oldresnumber=resnumber
        if line[:12]=='  #  RESIDUE' : read=True #line of the header, after the actual entries begin
    del line,wholefile
    return sequence,ss_kind_str

def population_from_dssp(filename, chain_to_grab=None, header=None,amino_acids_pool=amino_list1):
    '''
    Used by population_from_dssp_folder()
    get secondary structure populations from dssp file
    returns an ordered dictionary of numpy arrays of int
    '''
    if '.gz' in filename : f= gzip.open(filename,'r')
    else : f=open(filename,'r')
    wholefile= f.read().splitlines() # save everything in wholefile
    f.close()
    chain_to_avoid=[] # we need this as the dssp format is wrong for pdbs with more than 10000 (10^4) residues. We need to exclude those chains with this label
    resnumber,_,_,chain_id,_,_,_,_ = parse_dssp_line(wholefile[-1], warn_file=sys.stderr,get_hb_patterns=True, warn_tag=filename) # process last line
    i=0
    if resnumber > 9999 :
        while resnumber>9999 :
            wholefile.pop()
            i+=1
            resnumber,_,_,chain_id,_,_,_,_ = parse_dssp_line(wholefile[-1], warn_file=sys.stderr,get_hb_patterns=True, warn_tag=filename)
        chain_to_avoid+=[chain_id] # we exclude the chain between the numbering 10000 NOTE THAT EXTRA CHAIN BETA STRANDS COULD MAKE A MESS IN THIS CASE (hopefully the noise they introduce is negligible)
        sys.stderr.write('     file '+filename+' excluding %d lines from the end as resnumber exceeded 10000, also excluding chain %s\n' % (i,chain_id) )
   
    read=False
    if header is None : hd={'total_count':0,'B': 3, 'E': 2, 'G': 4, 'I': 5, 'H': 1, 'S': 7, 'T': 6}
    else : hd=header
    aa_ss={}
    for aa in amino_acids_pool :
        aa_ss[aa]=numpy.zeros(len(hd),'int') # use numpy as we can easily sum these when looping over multiple files
    # read the file, save the sequence and the strand map
    for line in wholefile : #we read for specific position,as split length can vary in some lines
        if read and line!='' :
            resnumber,resname,_,chain_id,ss_kind,_,_,_ = parse_dssp_line(line, warn_file=None,get_hb_patterns=False, warn_tag=filename)
            if chain_to_avoid is not None and chain_id in chain_to_avoid : continue
            if chain_to_grab is not None and chain_id !=chain_to_grab : continue
            if resname in aa_ss : 
                aa_ss[resname][hd['total_count']] += 1
                if ss_kind in hd :
                    aa_ss[resname][hd[ss_kind]] += 1
                #if ss_kind=='E' :
        if line[:12]=='  #  RESIDUE' : read=True #line of the header, after the actual entries begin
    del line,wholefile
    return aa_ss,hd


def population_from_dssp_folder(input_dir, chain_list_file=None,  out_to_pickle=True,amino_acids_pool=amino_list1, warnfile_name='stderr'):
    '''
    determines secondary structures population statistics, not useful to generate complementary peptides but used to reparametrize zyggregator and other staff.
    it returns data,hd and saves results in ss_populations.tsv
    it pickles output using protocol  in population_dict.pkl
    '''
    if chain_list_file!='' and chain_list_file is not None : use_chain_list=True
    else : use_chain_list= False
    if input_dir[-1]!='/' : input_dir +='/'
    file_list=os.listdir(input_dir)
    N=len(file_list)
    sys.stdout.write('%d files found in folder %s\n\n' % (N,input_dir) )
    if warnfile_name=='stderr' : warnings=sys.stderr
    elif warnfile_name=='stdout' : warnings=sys.stdout
    else :
        warnings=open(warnfile_name,'w')
    warnings.write('warning file from analyze_dssp_folder(%s)\n\n' %(input_dir))
    if use_chain_list : # use list with pdbname and the chains we are interested in
        chain_list={}
        if type(chain_list_file) is str and os.path.isfile(chain_list_file): pass
        elif os.path.isfile(input_dir+chain_list_file) : chain_list_file=input_dir+chain_list_file
        chain_file=file(chain_list_file).read().splitlines()
        for line in chain_file :
            if line[0]!='#' and line[0]!=' ' and line[0]!='\n' :
                try :
                    chain_list[ line.split()[0].split('.')[0] ] = line.split()[1].split(';')[:-1] #save so that pdb_id is the key and the entry is a list of chain ids
                except IndexError :
                    sys.stderr.write('**Index error in chain list %s while processing line :\n\t|%s|\n' % (chain_list_file,line))
                    pass
        sys.stderr.write('%d files in chain_list\n\n' % len(chain_list) )
        N=len(chain_list)
    #INITIALIZE OUTPUT
    hd={'total_count':0,'B': 3, 'E': 2, 'G': 4, 'I': 5, 'H': 1, 'S': 7, 'T': 6}
    aa_ss={}
    for aa in amino_acids_pool :
        aa_ss[aa]=numpy.zeros(len(hd),'int') # use numpy as we can easily sum these when looping over multiple files
   
    run_on_chain=[ None ] # in this case it will run on everything
    processed_pdb_id=[]
    i=0
    skip_exception=0
    just_another_file_format=0
    for dssp_file in file_list :
        if '.swp' in dssp_file or '.dssp' not in dssp_file or dssp_file[0]=='.': 
            just_another_file_format+=1
            continue
        _,pdb_id,_=misc.get_file_path_and_extension(dssp_file)
        if use_chain_list :
            if pdb_id.upper() not in chain_list : continue
            run_on_chain=chain_list[pdb_id.upper()]
        if i%40==0 :
            sys.stdout.write('%4d -> %5.2lf %% |' % (i,100.*i/(1.*N)) )
            sys.stdout.flush()
        i+=1
        try :
            for ch_id in run_on_chain : # note that if ch_id=None population_from_dssp() will run on all the chains
                aa_ss_file,_=population_from_dssp(input_dir+dssp_file, chain_to_grab=ch_id, header=hd,amino_acids_pool=amino_acids_pool)
                for k in aa_ss_file :
                    aa_ss[k]+=aa_ss_file[k] # sum up numpy lists
                    
            processed_pdb_id+=[ pdb_id.upper() ]
        except Exception :
            warnings.write( '\n**ERROR** while processing %s ch=%s\n File skipped\n' % (dssp_file,str(ch_id)))
            print_exc(file=warnings)
            warnings.flush()
            skip_exception+=1
            continue
    
 
    #pickle output using protocol 0.
    if out_to_pickle :
        output=open('ss_populations.pkl','wb')
        pickle.dump(aa_ss, output)
        output.close()
    # convert output to Data() class
    data=csv_dict.Data()
    data.hd=hd.copy()
    for aa in aa_ss :
        data[aa]=list(aa_ss[aa]) # convert in standard list format
    data.Print('ss_populations.tsv') # save in excel readable format
    if len(processed_pdb_id)!=N :
        sys.stdout.write('\n\n%d files correctly analyzed among %d expected (%d in folder but %d of different format)...\n' % (len(processed_pdb_id),N-just_another_file_format,N,just_another_file_format))
        if use_chain_list :
            for pdb_id in chain_list :
                if pdb_id not in processed_pdb_id :
                    sys.stdout.write('%s ' % pdb_id)
            sys.stdout.write('\n\n')
    sys.stdout.write('\n\n%d entries in pickle, corresponding to %d file processed and %d amino acids\n' %(len(aa_ss),i,sum(data.column('total_count'))))
    if warnfile_name!='stderr' and warnfile_name!='stdout' : warnings.close()
    return data,hd      
                
                
                
                
                
                
 
def assess_direction(comp_ids_i) :
    '''
    auxiliary function to determine the direction of a beta strand (parallel or antiparallel). Use only within beta_from_dssp()
    '''
    a=0
    b=-1
    while comp_ids_i[a]==0 : a+=1
    while comp_ids_i[b]==0 : b-=1
    if len(comp_ids_i)+b>a :   # this could be wrong in case the ids are not consecutive, later we fix it (see below and resnumber 115-123 in 1kt4.dssp)
        if comp_ids_i[a]<comp_ids_i[b] : return '+' 
        elif comp_ids_i[a]>comp_ids_i[b]: return '-'


def beta_from_dssp(filename, chain_to_grab=None, only_longer_than=2,warn_file=sys.stderr,debug=False, get_aa_abundance=False):
    '''
    The MAIN function that parses a dssp file, extracts beta-strand information and tries to handle whole exceptions!
    
    # NOTE: some dssp file contains lower letter aa such as 'i','f','j',... e.g. 2dsf.dssp
    # NOTE1: some strands are between different chains
    # NOTE2: some residues (e.g. 1aon.dssp resnumber 3667) are in a strand but don't have any partner.
    #        in such a case we give it as a partner '-', if it is in the middle of a beta sheet (hopefully unluckily) we give it two '-', '-'.
    # NOTE3 UNADDRESSED (or maybe now addressed): in some file two residues complement on a single one 1xr0.dssp 78,79 both on 92, 1blj.dssp  76,77 both on 85, 1ax3.dssp 58,59 both on 73 in 2lr4.dssp 13,14 both on 46
    
    this function parses into a dssp file (secondary structure info about a protein, as returned by the dssp program)
    optimized for .gz files as well
    it returns (check, not updated) two lists of str and strand_map. One with the sequences that participate in beta strands and the other one with the complementary bit.
      always read info from one of the two list otherwise the information is redundant. (if a protein has a single beta strand between segment A and B
      the function will return [A , B], [B , A])
      if a strand is part of a beta sheet (e.g. it pairs with one segment per side) it will appear twice. E.g. [A,B,C,C],[C,C,A,B] (this means that segment
      C is between A and B in a beta strand) 
    strand_map is a list of lists of lists which contain the resnumber of the strand in triplets [ [[10,24,38],..],[next strand in triplets],...] 
      means that resnumber 10 is in a beta sheet between resnumber 24 and 39.
    get_aa_abundance can either be a bool or a dict, if it is a dict the keys must be aa name (1 letter) and the entries are just updated (new occurrences are added to previous values).
      if True or dict a dictionary whose key are resname (1 letter) and values are the numbers they have been found in the file is returned
    '''
    if '.gz' in filename : f= gzip.open(filename,'r')
    else : f=open(filename,'r')
    if debug :
        if warn_file!=sys.stdout and type(warn_file) is file :
            debug_out=warn_file
        else :
            debug_out=sys.stderr
    wholefile= f.read().splitlines() # save everything in wholefile
    f.close()
    chain_to_avoid=[] # we need this as the dssp format is wrong for pdbs with more than 10000 (10^4) residues. We need to exclude those chains with this label
    resnumber,_,_,chain_id,_,_,_,_ = parse_dssp_line(wholefile[-1], warn_file=warn_file,get_hb_patterns=True, warn_tag=filename)
    i=0
    if resnumber > 9999 :
        while resnumber>9999 :
            wholefile.pop()
            i+=1
            resnumber,_,_,chain_id,_,_,_,_ = parse_dssp_line(wholefile[-1], warn_file=warn_file,get_hb_patterns=True, warn_tag=filename)
        chain_to_avoid+=[chain_id] # we exclude the chain between the numbering 10000 NOTE THAT EXTRA CHAIN BETA STRANDS COULD MAKE A MESS IN THIS CASE (hopefully the noise they introduce is negligible)
        warn_file.write('     file '+filename+' excluding %d lines from the end as resnumber exceeded 10000, also excluding chain %s\n' % (i,chain_id) )
    strand_seq=[]
    strand_ids=[] #list of list with resunmber(s) corresponding to the amino acids in stran_seq
    direction=[]
    complementary_seq=[]
    complementary_ids=[]#list of list with resunmber(s) corresponding to the amino acids in complementary_seq
    hb_patterns={}
    strand_map=[] # a list of lists of lists which contain the resnumber of the strand in triplets [ [[10,24,38],..],[next strand in triplets],...] means that resnumber 10 is in a beta shett between resnumber 24 and 39.
    read=False
    is_strand=False
    sequence=''
    processed_sequence=''
    strand=''
    first_res_num=None
    compute_abundance=False
    if get_aa_abundance==True or get_aa_abundance=={} : 
        get_aa_abundance={}
        for res in amino_list1 :
            get_aa_abundance[res]=0
        compute_abundance=True
    elif type(get_aa_abundance) is dict : compute_abundance=True
    # read the file, save the sequence and the strand map
    for line in wholefile : #we read for specific position,as split length can vary in some lines
        if read and line!='' :
            resnumber,resname,_,chain_id,ss_kind,beta_partner_1,beta_partner_2,hb_pattern = parse_dssp_line(line, warn_file=None,get_hb_patterns=True, warn_tag=filename)
            if sequence=='' : first_res_num=resnumber # save the fist resnumber to map the strand_map on the sequence (particularly useful if we are reading only one chain)
            sequence+= resname
            if chain_to_avoid is not None and chain_id in chain_to_avoid : continue
            if chain_to_grab is not None and chain_id not in chain_to_grab : continue #MOVE THIS 2 ROWS ABOVE IF YOU DON'T WANT EXTRA-CHAIN PARTNERS!!
            if compute_abundance and resname in get_aa_abundance : get_aa_abundance[resname]+=1
            processed_sequence+=resname
            if ss_kind=='E' :
                hb_patterns[resnumber]=hb_pattern
                #if resname not in amino_list1 : continue #avoid saving non-standard amino acids
                if not is_strand : # begin a new strand
                    strand_map+= [ [ [resnumber, beta_partner_1, beta_partner_2 ] ] ]
                    is_strand=True
                else : #add info to current strand
                    strand_map[-1]+=[ [resnumber, beta_partner_1, beta_partner_2 ] ]
            else : is_strand=False            
        if line[:12]=='  #  RESIDUE' : read=True #line of the header, after the actual entries begin
    del line,wholefile
    #now match resumber stored in strand_map to the sequence
    for strand in strand_map :
        strand_seq1=[''] # sequence (read form stand_map[][0]) that pairs with residues in beta_partner1
        comp_seq1=['']  # pairing sequence in  beta_partner1
        comp_ids1=[[]]  # ids (resnumber) of pairing (complementary) sequence
        strand_seq2=['']  # sequence (read form stand_map[][0]) that pairs with residues in beta_partner2
        comp_seq2=['']  # complementary pairing sequence in  beta_partner2
        comp_ids2=[[]] # ids (resnumber) of pairing (complementary) sequence
        direction2=[] # +/- to distinguish parallel from antiparallel
        for i,residue in enumerate(strand) :
            #check if it's part of a sheet or if it's just a strand or if NOTE2 in the header
            if residue[1]!=0 and (residue[1]- first_res_num)>=0 and (residue[1]- first_res_num)<=len(sequence) :
                strand_seq1[-1]+=sequence[ residue[0]- first_res_num ]
                comp_seq1[-1]+=sequence[ residue[1]- first_res_num ]
                comp_ids1[-1]+=[ residue[1] ]
            elif residue[1]==0 and residue[2]!=0 and comp_seq1[-1]!='' : #this means that the fragment forms strands on one side and then on the other
                strand_seq1+=['']
                comp_seq1+=['']
                comp_ids1+=[[]]
            if residue[2]!=0 and  (residue[2]- first_res_num)>=0 and (residue[2]- first_res_num)<=len(sequence) : 
                strand_seq2[-1]+=sequence[ residue[0]- first_res_num ]
                comp_seq2[-1]+=sequence[ residue[2]- first_res_num ]
                comp_ids2[-1]+=[ residue[2] ]
            elif residue[2]==0 and residue[1]!=0 and comp_seq2[-1]!='' : #this means that the fragment forms strands on one side and then on the other
                strand_seq2+=['']
                comp_seq2+=['']
                comp_ids2+=[[]]
            if residue[1] == 0 and residue[2]==0 : #this is a rare case where a residue is part of a beta strand but doesn't have a partner (NOTE2 in function header)
                if i>0 and i<len(strand)-1 : #if we are not looking at the extrema of the strand
                    if strand[i-1][1]!=0 and strand[i+1][1]!=0 : # both next and previous residues have partners in 1
                        strand_seq1[-1]+=sequence[ residue[0]- first_res_num ]
                        comp_seq1[-1]+='-'
                        comp_ids1[-1]+=[ residue[1] ] # that is 0
                    if strand[i-1][2]!=0 and strand[i+1][2]!=0 : # both next and previous residues have partners in 2
                        strand_seq2[-1]+=sequence[ residue[0]- first_res_num ]
                        comp_seq2[-1]+='-'
                        comp_ids2[-1]+=[ residue[2] ] # that is 0
                    elif strand[i+1][2]==0 and strand[i+1][1]==0: # this can happen (see 2VVZ res 840 and other position. also 3blx.dssp and many more). We just split the sequence
                        if strand_seq1[-1]!='':
                            strand_seq1+=['']
                            comp_seq1+=['']
                            comp_ids1+=[[]]
                        if strand_seq2[-1]!='':    
                            strand_seq2+=['']
                            comp_seq2+=['']
                            comp_ids2+=[[]]

                        #warn_file.write('WARNING in beta_from_dssp() file '+filename+' strange map for resnumber and pairs: '+str(residue)+' Check it out!\n')
                        #warn_file.flush()
                #else :
                    #warn_file.write('WARNING in beta_from_dssp() file '+filename+' strange map for resnumber and pairs: '+str(residue)+' Check it out!\n')
                    #warn_file.flush()
        # Merge the entries to strand_seq and complementary_seq (rather than keeping comp_seq1 and comp_seq2
        # also check that complementary ids are consecutive (the case in NOTE2 can be symmetrical)
        # and also assing parallel or antiparallel label to the strand 
        direction1=['?' for i in range(0,len(strand_seq1)) ] # +/- to distinguish parallel from antiparallel
        direction2=['?' for i in range(0,len(strand_seq2)) ]           
        for i in range(0,len(strand_seq1)) :
            if  strand_seq1[i]!='' and comp_seq1[i]!='' and len(strand_seq1[i])>only_longer_than :
                
                direction1[i]=assess_direction(comp_ids1[i])
                    #if comp_ids1[i][0]<comp_ids1[i][-1] : direction1[i]='+'
                    #elif comp_ids1[i][0]>comp_ids1[i][-1]: direction1[i]='-'
###
                if debug : debug_out.write( '1a:'+strand_seq1[i]+'<->'+comp_seq1[i]+str(comp_ids1[i])+'  %d  %s\n' % (i,direction1[i]) )
                #check that complementary ids are consecutive (the case in NOTE2 can be symmetrical)
                tmp1=[ strand_seq1[i] ]
                tmpC=[ comp_seq1[i] ]
                tmp_dir=[ direction1[i]]
                old=-1
                off=0
                for j,rid in enumerate(comp_ids1[i]) :
                    if rid==0 : # this means that an extra residue is present in this sequence ('-' in the comp, NOTE2)
                        if j>0 and j<len(comp_ids1[i])-1 :
                            if abs(comp_ids1[i][j-1]-comp_ids1[i][j+1])==2 : #sometimes e.g. 3995-4008 in 1aon.dssp) they are part of a strand but they do not form H bond.. hence are not reported as pairing
                                tmp1[-1] =  tmp1[-1][:j+off]+tmp1[-1][j+off].lower()+tmp1[-1][j+off+1:] 
                                if tmpC[-1][j+off]=='-':
                                    missing_res=sequence[(comp_ids1[i][j-1]+comp_ids1[i][j+1])/2- first_res_num]
                                    tmpC[-1] = tmpC[-1][:j+off]+missing_res.lower()+tmpC[-1][j+off+1:]
#
                                warn_file.write("\tfile %s changing to pair %s %s comp_ids %s\n" % (filename,tmp1[-1],tmpC[-1],str(comp_ids1[i])))
                                old=-1
                            elif abs(comp_ids1[i][j-1]-comp_ids1[i][j+1])>2 : # we need to split the sequence as in 2yi3 at resnumber 231
                                old= comp_ids1[i][j-1] # we do it with the framework written below at the next cycle
                        continue
                    if old!=-1 and abs(rid-old)>1 : # ids are not consecutive
                        if abs(rid-old)==2 :  # there must be a single residue in the complementary sequence that is squeezed so that it doesn't pair with these but is still part of the strand, insert a dash
                            tmp1[-1] =  tmp1[-1][:j+off]+'-'+tmp1[-1][j+off:] 
                            missing_res=sequence[(rid+old)/2- first_res_num]
                            tmpC[-1] = tmpC[-1][:j+off]+missing_res+tmpC[-1][j+off:] 
                            off+=1
                        else : # rare case (e.g. resnumber 296-301 of 3dfq.dssp) where the beta sheet topology is so complicated that strands that actually are at the two sides have to be recorded both as beta_partner_1 (or 2)
                            # split the sequence and write the warning...
                            if tmp1[-1]!=strand_seq1[i] :
                                #print '--a %d %s <-> %s   %d' % (off,str([ tmp1[-1][:j+off], tmp1[-1][j+off:] ] ),str([ tmpC[-1][:j+off], tmpC[-1][j+off:] ]),j)
                                tmp1[-1:]=[ tmp1[-1][:j+off], tmp1[-1][j+off:] ] 
                                tmpC[-1:]=[ tmpC[-1][:j+off], tmpC[-1][j+off:] ]
                                off=-1*(j) #we reset eventual +1 from above as they shouldn't be accounted for in the new sequence

                            else  :
                                #print '--- %d %s <-> %s   %d' % (off,str([ strand_seq1[i][:j],strand_seq1[i][j:] ] ),str([ comp_seq1[i][:j],comp_seq1[i][j:] ]),j)
                                tmp1[-1:] = [ strand_seq1[i][:j],strand_seq1[i][j:] ] 
                                tmpC[-1:] = [ comp_seq1[i][:j],comp_seq1[i][j:] ]
                                off=-1*(j)
                            st_d=j-len(tmpC[-2])
                            if st_d<0:st_d=0
                            tmp_dir[-1:]= [ assess_direction(comp_ids1[i][st_d:j]) , assess_direction(comp_ids1[i][j:]) ]
                            #warn_file.write('   split1 in beta_from_dssp() file '+filename+' complementary ids 1 not consecutive in strange way:\n    %s\n    now %s \n        %s  dir:%s\n' % (str(strand),str(tmp1),str(tmpC),str(tmp_dir)) )
                            #warn_file.flush()
                        
                    old=rid
                    
                if len(tmp1)!=len(tmpC) : 
                    warn_file.write("WARNING 1 in file %s ch=%s while splitting strand %s, len(tmp1)!=len(tmpC) %d %d  %s %s\n" % (filename,str(chain_to_grab),str(strand),len(tmp1),len(tmpC),str(tmp1),str(tmpC)) )
                strand_seq += [seq for seq in tmp1 if len(seq)>only_longer_than ]
                del seq
                complementary_seq += [seq for seq in tmpC if len(seq)>only_longer_than ]
                if len(tmp_dir)!=len(tmpC) : 
                    warn_file.write("WARNING 1 in file %s ch=%s while splitting strand %s, len(tmp_dir)!=len(tmpC) %d %d  %s %s\n" % (filename,str(chain_to_grab),str(strand),len(tmp_dir),len(tmpC),str(tmp_dir),str(tmpC)) )
                direction+= [ tmp_dir[k] for k,seq in enumerate(tmp1) if len(seq)>only_longer_than ]
###
                if debug : debug_out.write( '1b:'+str(len(strand_seq))+str(strand_seq[-len(tmp1):])+'<->'+str(complementary_seq[-len(tmp1):])+'  %d  %s\n' % (i,str(direction[-len(tmp1):])) )
                if debug and len(tmp1)>1 :
                    debug_out.write( str(tmp1)+'  '+str(tmpC)+'!!!!!!!!\n')

        del i
        for i in range(0,len(strand_seq2)) :
            if  strand_seq2[i]!='' and comp_seq2[i]!='' and len(strand_seq2[i])>only_longer_than :
                direction2[i]=assess_direction(comp_ids2[i])
###
                if debug : debug_out.write( '2a:'+strand_seq2[i]+'<->'+comp_seq2[i]+str(comp_ids2[i])+'  %d  %s\n' % (i,direction2[i]) )
                tmp1=[ strand_seq2[i] ]
                tmpC=[ comp_seq2[i] ]
                tmp_dir=[ direction2[i]]
                old=-1
                off=0
                for j,rid in enumerate(comp_ids2[i]) :
                    if rid==0 : # this means that an extra residue is present in this sequence ('-' in the comp, NOTE2)
                        if j>0 and j<len(comp_ids2[i])-1 :
                            if abs(comp_ids2[i][j-1]-comp_ids2[i][j+1])==2 : #sometimes e.g. 3995-4008 in 1aon.dssp) they are part of a strand but they do not form H bond.. hence are not reported as pairing
                                tmp1[-1] =  tmp1[-1][:j+off]+tmp1[-1][j+off].lower()+tmp1[-1][j+off+1:]
                                if tmpC[-1][j+off]=='-':
                                    missing_res=sequence[(comp_ids2[i][j-1]+comp_ids2[i][j+1])/2- first_res_num]
                                    tmpC[-1] = tmpC[-1][:j+off]+missing_res.lower()+tmpC[-1][j+off+1:]
#
                                warn_file.write("\tfile %s changing to pair %s %s comp_ids %s\n" % (filename,tmp1[-1],tmpC[-1],str(comp_ids2[i])))
                                old=-1
                            elif abs(comp_ids2[i][j-1]-comp_ids2[i][j+1])>2 : # we need to split the sequence as in 2yi3 at resnumber 231
                                old= comp_ids2[i][j-1] # we do it with the framework written below at the next cycle
                        continue
                    if old!=-1 and abs(rid-old)>1 : # ids are not consecutive
                        if abs(rid-old)==2 :  
                            tmp1[-1] =  tmp1[-1][:j+off]+'-'+tmp1[-1][j+off:] 
                            missing_res=sequence[(rid+old)/2- first_res_num]
                            tmpC[-1] = tmpC[-1][:j+off]+missing_res+tmpC[-1][j+off:]
                            off+=1
                        else : # rare case (e.g. resnumber 296-301 of 3dfq.dssp) where the beta sheet topology is so comlicated that strands that actually are at the two sides have to be recorde both as beta_partner_1 (or 2)
                            # split the sequence and write the warning...
                            if tmp1[-1]!=strand_seq2[i] :
                                tmp1[-1:]=[ tmp1[-1][:j+off], tmp1[-1][j+off:] ] 
                                tmpC[-1:]=[ tmpC[-1][:j+off], tmpC[-1][j+off:] ]
                                off=-1*(j) #we reset eventual +1 from above as they shouldn't be accounted for in the new sequence
                            else  :
                                tmp1[-1:] = [ strand_seq2[i][:j],strand_seq2[i][j:] ] 
                                tmpC[-1:] = [ comp_seq2[i][:j],comp_seq2[i][j:] ]
                                off=-1*(j)
                            st_d=j-len(tmpC[-2])
                            if st_d<0:st_d=0
                            tmp_dir[-1:]= [ assess_direction(comp_ids2[i][st_d:j]) , assess_direction(comp_ids2[i][j:]) ]
                            #warn_file.write('   split2 in beta_from_dssp() file '+filename+' complementary ids 2 not consecutive in strange way:\n    %s\n    now %s \n        %s\n' % (str(strand),str(tmp1),str(tmpC)) )
                            #warn_file.flush()
                    old=rid
                if len(tmp1)!=len(tmpC) : 
                    warn_file.write("WARNING 2 in file %s ch=%s while splitting strand %s, len(tmp1)!=len(tmpC) %d %d  %s %s\n" % (filename,str(chain_to_grab),str(strand),len(tmp1),len(tmpC),str(tmp1),str(tmpC)) )
                strand_seq += [seq for seq in tmp1 if len(seq)>only_longer_than ]
                complementary_seq += [seq for seq in tmpC if len(seq)>only_longer_than ]
                if len(tmp_dir)!=len(tmpC) : 
                    warn_file.write("WARNING 2 in file %s ch=%s while splitting strand %s, len(tmp_dir)!=len(tmpC) %d %d  %s %s\n" % (filename,str(chain_to_grab),str(strand),len(tmp_dir),len(tmpC),str(tmp_dir),str(tmpC)) )
                direction+= [ tmp_dir[k] for k,seq in enumerate(tmp1) if len(seq)>only_longer_than ]
###
                if debug : debug_out.write('2b:'+str(len(strand_seq))+str(strand_seq[-len(tmp1):])+'<->'+str(complementary_seq[-len(tmp1):])+'  %d  %s\n' % (i,str(direction[-len(tmp1):])))
                if debug and len(tmp1)>1 :
                    debug_out.write( str(tmp1)+'  '+str(tmpC)+'!!!!!!!!\n') 
                    
    # clean possible entries with X or other unknown amino aicds. THIS ARE ALSO MODRES in .dssp files
    tmp1=[]
    tmpC=[]
    tmpD=[]
    keep_list=[a.lower() for a in amino_list1 ]
    keep_list+= amino_list1+['-']
    for i,seq in enumerate(strand_seq) :
        bits= mybio.split_string_at(seq,keep_list,complementary_seq[i],only_longer_than=only_longer_than)
        tmp1+=bits[0]
        tmpC+=bits[1]
        tmpD+=[ direction[i] ]*len(bits[0]) 
    strand_seq=tmp1
    complementary_seq=tmpC
    direction=tmpD
#
    #print tmpD,tmp1,tmpC
    if debug : debug_out.write( " length of new strand,comp, direction: %d %d %d\n" % (len(tmp1),len(tmpC),len(tmpD)) ) 
    del tmpD,tmp1,tmpC
    # now check that there are no pairs with a '-' in both. sometimes e.g. 3995-4008 in 1aon.dssp) they are part of a strand but they do not form H bond.. hence are not reported as pairing
    strange_strands=[]    
    for i,seq in enumerate(strand_seq) :
        if '-' in seq and '-' in complementary_seq[i]: # it shouldn't get in here any more as I added this correction before, I leave this bit just in case (this is why I added the last endif with.upper())
            j=seq.index('-')
            k=complementary_seq[i].index('-')
            if abs(j-k)<=1 : #if the dash represents complementary amino acids
#
#                print 'STRAND:',strand_seq[i]
#                print 'COMPLE:',complementary_seq[i]
               
                warn_file.write( '\tchanging: %s  %s  %s %s\n' % (str(strand_seq[i]), str(complementary_seq[i]),str(i),filename) )
#               print i,j,k
                strand_seq[i] = strand_seq[i][:j]+strand_seq[i][j+1].lower()
                if j+2<len(strand_seq[i]) : strand_seq[i] += strand_seq[i][j+2:] #add next bit
                complementary_seq[i] = complementary_seq[i][:k-1]+complementary_seq[i][k-1].lower()
                if k+1<len(complementary_seq[i]): complementary_seq[i]+= complementary_seq[i][k+1:] #add next bit
#                
            warn_file.write( '\tstrange   '+str(strand_seq[i])+'   '+str(complementary_seq[i])+' '+filename+' '+str(i)+'\n' ) 
            strange_strands+=[i]
        elif '-' in seq or '-' in complementary_seq[i]:
            strange_strands+=[i]
        elif seq!=seq.upper() or  complementary_seq[i]!=complementary_seq[i].upper(): # if it contains lower case letters
            strange_strands+=[i]
            
            
            
    # we now map our strands on the resnumbers to incorporate the hydrogen bond patterns. It was clearly smarter to this before but this is a latter addition!
    offset=0  # this is useful only when reading specific chains of very large complexes (e.g. 4diw.dssp) so that the first strand_seq is more readily found
    try :
        if strand_map[0][0][0]>500 :
            offset=strand_map[0][0][0]-500
    except :
        pass
    template_seq=sequence[offset:]
    ids_pool=list(range(first_res_num,first_res_num+len(sequence)))
    template_pool=ids_pool[offset:]
    if debug : debug_out.write( 'now matching hb_patterns, we have a hb_pattern dict of len %d and a sequence of len %d and %d strands. using offset= %d hence first id is %d ->actual seq_len=%d\n' % (len(hb_patterns),len(sequence),len(strand_seq),offset,ids_pool[0],len(template_seq)) )
    for i,seq in enumerate(strand_seq) :
        if debug : debug_out.write( '     %d %s <--> %s\n' % (i,strand_seq[i],complementary_seq[i]) )
        notfound=True
        oldj=-5
        oldk=-5
        dash_pos=None
        comp_dash_pos=None
        comp_seq = complementary_seq[i]
        jump_count=0
        debug_count=0
        last_considered_pos=0
        while notfound :
            debug_count+=1
            if debug_count>500:
                warn_file.write('ERROR while matching hb patterns, too many void iterations in %d,%s from sequence in file %s ch=%s (failed %d len_seq %d len_template %d)\n' % (i,seq,filename,str(chain_to_grab),jump_count,len(sequence),len(template_seq)))
                strand_ids+=[ ['?']]
                complementary_ids+=[ ['?']]
                break
            if i in strange_strands :
                if '-' in seq : # recall that - can be also at the beginning/end in the rare case where the complementary has an X next to it (see 2vne)
                    dash_pos=[k for k in range(0,len(seq)) if seq[k]=='-']
                    seq=seq.replace('-','')
                if '-' in comp_seq :
                    comp_dash_pos=[k for k in range(0,len(comp_seq)) if comp_seq[k]=='-']
                    comp_seq=comp_seq.replace('-','')
                if not seq.isupper() :
                    if seq[0].isupper() : # if it's not at the begining than it's not a problem
                        seq=seq.upper()
                    #think about it...
                    else :
                        if dash_pos is None : dash_pos=[0]
                        elif 0 not in dash_pos : dash_pos[:0]=[0]
                        if comp_dash_pos is None : comp_dash_pos=[0] # we are not going to find the compelementary id anyway
                        elif 0 not in comp_dash_pos : comp_dash_pos[:0]=[0]
                        seq=seq[1:].upper()
                        comp_seq=comp_seq[1:].upper()
                        #strand_ids+=[ ['?']]
                        #complementary_ids+=[ ['?']]
                        #break
                if not comp_seq.isupper() :
                    if comp_seq[0].isupper() :
                        comp_seq=comp_seq.upper()
                    #think about it...
                    else :
                        if dash_pos is None : dash_pos=[0]
                        elif 0 not in dash_pos : dash_pos[:0]=[0]
                        if comp_dash_pos is None : comp_dash_pos=[0] # we are not going to find the compelementary id anyway
                        elif 0 not in comp_dash_pos : comp_dash_pos[:0]=[0]
                        seq=seq[1:].upper()
                        comp_seq=comp_seq[1:].upper()
                        #strand_ids+=[ ['?']]
                        #complementary_ids+=[ ['?']]
                        #break
            j=template_seq.find(seq,last_considered_pos) # this will look in the slice starting at last_considered_pos which is zero unless we have already looped once
            if debug and last_considered_pos>0:
                try : 
                    debug_out.write( 'new_pos %s\n' % (str((j,'id=',template_pool[j],oldj,last_considered_pos,template_seq[last_considered_pos:],seq,i,jump_count))) )
                except IndexError : 
                    debug_out.write( 'new_pos indexError in ids_pool | %s\n' % (str((j,oldj,last_considered_pos,seq,i,jump_count))) )
                    pass
            if j==oldj :
                if debug :
                    try : 
                        debug_out.write( 'double_found %s\n' % (str((j,'id=',template_pool[j],oldj,template_seq[j:],seq,i,jump_count,len(sequence),len(template_seq)))))
                    except IndexError : 
                        debug_out.write( 'new_pos indexError in ids_pool | %s\n' % (str((j,oldj,seq,i,jump_count,len(sequence),len(template_seq)))))
                        pass
                if oldk==-5: 
                    oldk=sequence[offset:].find(seq)
                    if oldk!=-1 and 0<oldj<len(template_pool) and ids_pool[offset:][oldk]<template_pool[oldj] :
                        if debug : debug_out.write( '  in double_found going back to %d id=%d\n' % (oldk,ids_pool[offset:][oldk]))
                        template_seq=sequence[offset:][oldk:]
                        template_pool=ids_pool[offset:][oldk:]
                        continue
                if jump_count>100 : 
                    warn_file.write('ERROR while matching hb patterns, cant retrieve resnumbers of %d,%s from sequence in file %s ch=%s (failed %d len_seq %d len_template %d)\n' % (i,seq,filename,str(chain_to_grab),jump_count,len(sequence),len(template_seq)))
                    strand_ids+=[ ['?']]
                    complementary_ids+=[ ['?']]
                    break
                jump_count+=1 
                last_considered_pos= j+1
                oldj=j
                continue
                
            if j<0:
                try : # try to go 150 residues backwards (just in case)
                    last_found=strand_ids[-1][0]
                    if type(last_found) is int :
                        template_seq=sequence[offset:][last_found-250:]
                        template_pool=ids_pool[offset:][last_found-250:]
                        if debug : debug_out.write( '   sliding 250 res backward when looking for %s len_template=%d\n' % (seq,len(template_seq)))
                        oldj=-5
                        continue
                except :
                    pass
                warn_file.write('ERROR while matching hb patterns, cant retrieve resnumbers of %d,%s %s from sequence in file %s ch=%s (failed %d len_seq %d len_template %d)\n' % (i,seq,direction[i],filename,str(chain_to_grab),jump_count,len(sequence),len(template_seq)))
                strand_ids+=[ ['?']]
                complementary_ids+=[ ['?']]
                break # if nothing found break from while
            if len([1 for k in template_pool[j:j+len(seq)] if k in hb_patterns])==len(seq) : # since hb_patterns has been filled with E ss only this should prevent from considering other parts of the sequence
                strand_ids+=[ template_pool[j:j+len(seq)] ]
                if dash_pos is not None :
                    for dp in dash_pos :
                        strand_ids[-1]=strand_ids[-1][:dp]+[None]+strand_ids[-1][dp:]
                # get complementaries
                first_cmp_index=-1
                c_off=0
                if strand_ids[-1][0] is not None : 
                    target_id=strand_ids[-1][0]
                    if complementary_seq[i][0]=='-' or complementary_seq[i][0].islower() :
                        target_id+=1
                elif strand_ids[-1][1] is not None : 
                    target_id=strand_ids[-1][1]
                    if complementary_seq[i][0]!='-' and complementary_seq[i][0].isupper() :
                        if direction[i] =='-' :c_off=1
                        else : c_off=-1
                    if complementary_seq[i][1]=='-' or complementary_seq[i][1].islower() :
                        target_id+=1
                else :
                    warn_file.write('ERROR while matching hb patterns, too many None entries in strand_ids of %d,%s -> %s from sequence in file %s ch=%s\n' % (i,strand_seq[i],str(strand_ids[-1]),filename,str(chain_to_grab)))
                    complementary_ids+=[ ['?']]
                    break
                if debug : debug_out.write( '  hb_m: target: %s\n' % (str((target_id,c_off,strand_seq[i],str(strand_ids[-1]),seq,complementary_seq[i],comp_seq,direction[i],i))) )
                for c,strand in enumerate(strand_map) :
                    for sco,s in enumerate(strand) :
                        if s[0]==target_id :
##
                            if debug : 
                                debug_out.write('  hb_m: target_found: %s\n' % (str((str(s),s[1]-first_res_num+c_off,s[2]-first_res_num+c_off,'strand_map id=',c,' target_id',target_id))) )
                                if s[1]>0 : debug_out.write('     s1=%d --> tar_comp=%d->%s of=%d to compare with comp_seq[0]=%s\n' % (s[1],s[1]-first_res_num+c_off,sequence[s[1]-first_res_num+c_off],c_off,comp_seq))
                                if s[2]>0 : debug_out.write('     s2=%d --> tar_comp=%d->%s of=%d to compare with comp_seq[0]=%s\n' % (s[2],s[2]-first_res_num+c_off,sequence[s[2]-first_res_num+c_off],c_off,comp_seq))
                            if s[1]>0 and sequence[s[1]-first_res_num+c_off]==comp_seq[0]  and s[2]>0 and sequence[s[2]-first_res_num+c_off]==comp_seq[0] : #rare case where a residues is sandwiched between two identical residues
                                if direction[i]=='-': 
                                    if sequence[s[2]-first_res_num+c_off-len(comp_seq)+1:s[2]-first_res_num+c_off+1]==sequence[s[1]-first_res_num+c_off-len(comp_seq)+1:s[1]-first_res_num+c_off+1]==comp_seq[::-1] : #extremely rare case (see resnumber 86 in 2gow.dssp
                                        if debug : debug_out.write('       peculiar sandwitching - between s1->%d  %s s2->%d %s\n' % (s[1],sequence[s[1]-first_res_num+c_off-len(comp_seq)+1:s[1]-first_res_num+c_off+1],s[2],sequence[s[2]-first_res_num+c_off-len(comp_seq)+1:s[2]-first_res_num+c_off+1]))
                                        try :
                                            if s[1] not in complementary_ids[-1] : first_cmp_index=s[1]-first_res_num+c_off
                                            else : first_cmp_index=s[2]-first_res_num+c_off
                                        except IndexError:
                                            first_cmp_index=s[1]-first_res_num+c_off
                                            pass
                                    elif sequence[s[1]-first_res_num+c_off-len(comp_seq)+1:s[1]-first_res_num+c_off+1]==comp_seq[::-1] : first_cmp_index=s[1]-first_res_num+c_off
                                    elif sequence[s[2]-first_res_num+c_off-len(comp_seq)+1:s[2]-first_res_num+c_off+1]==comp_seq[::-1] : first_cmp_index=s[2]-first_res_num+c_off
                                elif direction[i]=='+' :
                                    if sequence[s[1]-first_res_num+c_off:s[1]-first_res_num+c_off+len(comp_seq)]==sequence[s[2]-first_res_num+c_off:s[2]-first_res_num+c_off+len(comp_seq)]==comp_seq : #extremely rare case (see resnumber 86 in 2gow.dssp
                                        if debug : debug_out.write('       peculiar sandwitching + between s1->%d  %s s2->%d %s\n' % (s[1],sequence[s[1]-first_res_num+c_off:s[1]-first_res_num+c_off+len(comp_seq)],s[2],sequence[s[2]-first_res_num+c_off:s[2]-first_res_num+c_off+len(comp_seq)]))
                                        try :
                                            if s[1] not in complementary_ids[-1] : first_cmp_index=s[1]-first_res_num+c_off
                                            else : first_cmp_index=s[2]-first_res_num+c_off
                                        except IndexError:
                                            first_cmp_index=s[1]-first_res_num+c_off
                                            pass
                                    elif sequence[s[1]-first_res_num+c_off:s[1]-first_res_num+c_off+len(comp_seq)]==comp_seq : first_cmp_index=s[1]-first_res_num+c_off
                                    elif sequence[s[2]-first_res_num+c_off:s[2]-first_res_num+c_off+len(comp_seq)]==comp_seq : first_cmp_index=s[2]-first_res_num+c_off
                            elif s[1]>0 and sequence[s[1]-first_res_num+c_off]==comp_seq[0] :
                                first_cmp_index=s[1]-first_res_num+c_off
                            elif s[2]>0 and sequence[s[2]-first_res_num+c_off]==comp_seq[0] :
                                first_cmp_index=s[2]-first_res_num+c_off
                if first_cmp_index!=-1 :
                    if direction[i]=='+' :
                        if sequence[first_cmp_index:first_cmp_index+len(comp_seq)]==comp_seq :
                            complementary_ids+= [ ids_pool[first_cmp_index:first_cmp_index+len(comp_seq)]]
                            if comp_dash_pos is not None :
                                for dp in comp_dash_pos :
                                    complementary_ids[-1]=complementary_ids[-1][:dp]+[None]+complementary_ids[-1][dp:]
                        else :
                            tmp = template_seq.find(seq,j+1)
                            if  tmp>0 :
                                last_considered_pos= j+1
                                warn_file.write('   sliding over at %d, while matching hb patterns, found complementary start resnumbers (%d) + of %d,%s %s (now removed) complementary: %s but whole complment doesnt map on sequence in file %s ch=%s\n' % (tmp,first_cmp_index,i,strand_seq[i],str(strand_ids[-1]),complementary_seq[i],filename,str(chain_to_grab)))
                                strand_ids.pop()
                                continue
                            warn_file.write('ERROR while matching hb patterns, found complementary start resnumbers (%d) + of %d,%s %s complementary: %s but whole complment doesnt map on sequence in file %s ch=%s\n' % (first_cmp_index,i,strand_seq[i],str(strand_ids[-1]),complementary_seq[i],filename,str(chain_to_grab)))
                            complementary_ids+=[ ['?']]
                    elif direction[i]=='-' :
                        if sequence[first_cmp_index-len(comp_seq)+1:first_cmp_index+1]==comp_seq[::-1] :
                            complementary_ids+= [ ids_pool[first_cmp_index-len(comp_seq)+1:first_cmp_index+1][::-1] ]
                            if comp_dash_pos is not None :
                                for dp in comp_dash_pos :
                                    complementary_ids[-1]=complementary_ids[-1][:dp]+[None]+complementary_ids[-1][dp:]
                        else :
                            tmp = template_seq.find(seq,j+1)
                            if  tmp>0 :
                                last_considered_pos= j+1
                                strand_ids.pop()
                                warn_file.write('   sliding over at %d, while matching hb patterns, found complementary start resnumbers (%d) - of %d,%s %s (now removed) complementary: %s but whole complment doesnt map on sequence in file %s ch=%s\n' % (tmp,first_cmp_index,i,strand_seq[i],str(strand_ids[-1]),complementary_seq[i],filename,str(chain_to_grab)))
                                continue
                            warn_file.write('ERROR while matching hb patterns, found complementary start resnumbers (%d) - of %d,%s %s complementary: %s but whole complment doesnt map on sequence in file %s ch=%s\n' % (first_cmp_index,i,strand_seq[i],str(strand_ids[-1]),complementary_seq[i],filename,str(chain_to_grab)))
                            complementary_ids+=[ ['?']]
                else :
                    # in this case check that the seq we want is not afterwards along the sequence
                    tmp = template_seq.find(seq,j+1)
                    if  tmp>0 :
                        last_considered_pos= j+1
                        warn_file.write('   sliding over at %d, while matching hb patterns, cant find complementary resnumbers of %d,%s %s %s (now removed) complementary: %s from sequence in file %s ch=%s\n' % (tmp,i,strand_seq[i],str(strand_ids[-1]),direction[i],complementary_seq[i],filename,str(chain_to_grab)))
                        strand_ids.pop()
                        continue
                    warn_file.write('ERROR while matching hb patterns, cant find complementary resnumbers of %d,%s %s %s complementary: %s from sequence in file %s ch=%s\n' % (i,strand_seq[i],str(strand_ids[-1]),direction[i],complementary_seq[i],filename,str(chain_to_grab)))
                    complementary_ids+=[ ['?']]
                if debug : debug_out.write('   hb_m: comple_found: first_cmp_index=%d ids=%s\n' % (first_cmp_index,str(complementary_ids[-1])) )
                notfound=False
            oldj=j
            if j-30>0 :
                template_seq=template_seq[j-30:]  # so to avoid doulbe counting
                template_pool=template_pool[j-30:]
                oldj=30
                last_considered_pos=0

                    
    if not ( len(strand_seq)==len(complementary_seq)==len(direction)==len(strand_ids)==len(complementary_ids) ) :
        warn_file.write('ERROR unequal list length len(strand_seq)=%d len(complementary_seq)=%d len(direction)=%d len(strand_ids)=%d len(complementary_ids)=%d in file %s ch=%s\n' % (len(strand_seq),len(complementary_seq),len(direction),len(strand_ids),len(complementary_ids),filename,str(chain_to_grab)))
    warn_file.flush()
    if compute_abundance  : return strand_seq,complementary_seq,direction,strand_map,sequence,strange_strands,strand_ids,complementary_ids,hb_patterns, get_aa_abundance,processed_sequence
    return strand_seq,complementary_seq,direction,strand_map,sequence,strange_strands,strand_ids,complementary_ids,hb_patterns,processed_sequence


def debug_print_comp(strand_seq,complementary_seq,direction,strand_ids,complementary_ids,sequence,hb_trace=None,hb_energies=None, hb_patterns=None):
    '''
    auxiliary function that prints out the complemntaries in a debug friendly way
    '''
    for i,si in enumerate(complementary_ids) :
        if si!=['?'] :
            hb='\t'
            if hb_patterns is not None :
                for k in strand_ids[i] :
                    if k is not None : hb+=str(hb_patterns[k])
            print(i)
            k=0
            le=None
            while le is None : 
                le=strand_ids[i][k]
                k+=1
            k=-1
            re=None
            while re is None : 
                re=strand_ids[i][k]
                k-=1
            print(len(si),len(strand_seq[i]),strand_ids[i],'\t'+strand_seq[i],sequence[le-1:re],direction[i])
            if hb_trace is not None : print('\t\t\t\t'+str(hb_trace[i]))
            if hb_energies is not None : print('\t\t\t\t'+str(hb_energies[i]))
            if hb_patterns is not None : print(hb)
            k=0
            le=None
            while le is None : 
                le=si[k]
                k+=1
            k=-1
            re=None
            while re is None : 
                re=si[k]
                k-=1
            if le>re :
                print(len(si),len(complementary_seq[i]),si,'\t'+complementary_seq[i],sequence[(re-1):(le)])
            else :
                print(len(si),len(complementary_seq[i]),si,'\t'+complementary_seq[i],sequence[(le-1):(re)])
        else : print(str(['?']))
    return





def process_hb_pattern(strand_seq,sequence,strand_ids,complementary_ids,hb_patterns,direction, doublecheck=True,strange_strands=None,complement_seq=None, clear_threshold=-1.):
    '''
    CONSIDER ADDING EXTREMA FOR PARALLEL STRANDS
    returns hb_trace, hb_energies. The first is a  list of lists, outer list indeces are posistion in each_strand_seq,
      each inner list is [NtoO, OtoN] both can be 0->residue in front on the facing strand (mainly antiparallel), 1-> next residue, -1->residue before
      in corresponding complementary strand,  None-> no hbond
     hb_energies is the same but with corresponding energy values. This list structure is chosen as it suits the fragmentation process
     Energies are as read from the DSSP file
    '''
    problematic_patterns=[]
    hb_trace=[[] for k in range(0,len(strand_seq)) ] # list of lists, outer list indeces are posistion in each_strand_seq, each inner list is [NtoO, OtoN] both can be 0->residue in front (mainly antiparallel), 1-> next residue, -1->residue before in corresponding complementary strand,None-> no hbond
    hb_energies=[[] for k in range(0,len(strand_seq)) ] # same with energies
    warnings=''
    for i,st in enumerate(strand_seq) : # we use strand_seq and not strand_id as the latter can be ['?'] in some cases
        hb_trace[i]=[ [None,None] for k in range(0,len(st)) ]
        hb_energies[i]=[ [None,None]for k in range(0,len(st)) ]
        if strand_ids[i]!=['?'] and complementary_ids[i]!=['?'] :
            for j in range(0,len(st)) :
                if strand_ids[i][j] is not None and strand_ids[i][j] in hb_patterns :
                    for hbp in hb_patterns[strand_ids[i][j]] :
                        #comp_ids=complementary_ids[i][:]
                        jj=j
                        if direction[i]=='-' : comp_ids=complementary_ids[i][:]
                        else :
                            comp_ids=complementary_ids[i][:] # add one extra so that we know whether first/last res of parallel strand was engaged with HB outside the facing fragment or not 
                            if comp_ids[0] is not None : 
                                comp_ids=[comp_ids[0]-1]+comp_ids
                                jj=j+1
                            if comp_ids[-1] is not None : comp_ids= comp_ids+[comp_ids[-1]+1]
                        if hbp[0] in comp_ids : # loop on all the keys (each key is like (8231, 'NtoO')->(complementary_id, kind))
                            partner=None
                            if comp_ids[jj] is not None and hbp[0]== comp_ids[jj] : # None could be when analyzing strange_strand_dict
                                partner=0
                            elif jj+1<len(comp_ids) and comp_ids[jj+1] is not None and hbp[0]== comp_ids[jj+1] : 
                                partner=1
                            elif jj-1>=0 and comp_ids[jj-1] is not None and hbp[0]== comp_ids[jj-1] :
                                partner=-1
                            if partner is not None:
                                save_it=True
                                if hbp[1]=='NtoO' :
                                    if doublecheck  : # given the looping you may have already save this
                                        if direction[i]=='+' and hb_trace[i][j][0]==-1 and partner!=-1 :  # we would overwrite something that looks ok.
                                            save_it=False
                                        elif direction[i]=='-' and hb_trace[i][j][0]==0 and partner!=0 :
                                            save_it=False
                                    if save_it :
                                        hb_trace[i][j][0]=partner
                                        hb_energies[i][j][0]=hb_patterns[strand_ids[i][j]][hbp]
                                else : #save OtoN
                                    if doublecheck  :
                                        if direction[i]=='+' and hb_trace[i][j][1]==1 and partner!=1 : 
                                            save_it=False
                                        elif direction[i]=='-' and hb_trace[i][j][1]==0 and partner!=0 :
                                            save_it=False
                                    if save_it :
                                        hb_trace[i][j][1]=partner
                                        hb_energies[i][j][1]=hb_patterns[strand_ids[i][j]][hbp]
        else :
            warnings+='   in strand %d %s %s %s pattern are ? ? \n' % (i,strand_seq[i],direction[i],str(strand_ids[i]))
            problematic_patterns+=[i]
            hb_trace[i]=[ ['?','?'] ]
            hb_energies[i]=[ ['?','?'] ]
    # double check what was done if asked
    if doublecheck :
        possible_patterns=create_possible_hb_patterns( max([len(seq) for seq in strand_seq]) )
        for i,pattern in enumerate(hb_trace) :
            if '?' in pattern : 
                continue
            if strange_strands is not None and complement_seq is not None and i in strange_strands :
                s, c,bits = mybio.split_string_at(strand_seq[i], amino_list1, second_string=complement_seq[i], only_longer_than=2)
                hb_p_bits, map_bit_pos_hb= mybio.split_at_bits(bits, pattern, only_longer_than=2)
                hb_energy_bits, _= mybio.split_at_bits(bits, hb_energies[i][:], only_longer_than=2)
                #print bits
                #print pattern,' --> ',hb_p_bits,'map',map_bit_pos_hb,'\n seq ',s,c,strand_seq[i],complement_seq[i]
                for k, pattern_bit in enumerate(hb_p_bits) :
                    p_id, scores,warn_str = classify_pattern( pattern_bit, possible_patterns ,direction=direction[i], hb_energy=hb_energy_bits[k])
                    if warn_str!='' : warnings+=warn_str+' (energy_bit %d) strand %i %s dssp file ids %s\n' % (k,i,strand_seq[i],str(strand_ids[i]))
                    for j,pair in enumerate(pattern_bit) :
                        if pair[0]!= possible_patterns[p_id][j][0] or pair[1]!= possible_patterns[p_id][j][1] : # most of the times one is lists and one tuple so its easier to look at the individual elements 
                            # if there shouldn't be hb bond delete it if it is weak
                            if possible_patterns[p_id][j][0] is None and hb_energies[i][map_bit_pos_hb[k][j]][0] > clear_threshold :
                                hb_trace[i][map_bit_pos_hb[k][j]][0]=None
                                hb_energies[i][map_bit_pos_hb[k][j]][0]=None
                            if possible_patterns[p_id][j][1] is None and hb_energies[i][map_bit_pos_hb[k][j]][1] > clear_threshold :
                                hb_trace[i][map_bit_pos_hb[k][j]][1]=None
                                hb_energies[i][map_bit_pos_hb[k][j]][1]=None
                    if (direction[i]=='+' and p_id < 2) or (direction[i]=='-' and p_id >= 2) :
                        warnings+= '  in strand %d %s %s pattern %s bit %d %s maps on opposite direction! scores %s  ids %s\n' % (i,strand_seq[i],direction[i],pattern,k, str(pattern_bit),str(scores),str(strand_ids[i]))
                        problematic_patterns+=[i]
                 
            else :
                p_id, scores,warn_str = classify_pattern( pattern, possible_patterns,direction=direction[i] , hb_energy=hb_energies[i])
                if warn_str!='' : warnings+=warn_str+' strand %i %s dssp file ids %s\n' % (i,strand_seq[i],str(strand_ids[i]))
                for j,pair in enumerate(pattern) :
                    if pair[0]!= possible_patterns[p_id][j][0] or pair[1]!= possible_patterns[p_id][j][1] : 
                        # if there shouldn't be hb bond delete it if it is weak
                        if possible_patterns[p_id][j][0] is None and hb_energies[i][j][0] > clear_threshold :
                            hb_trace[i][j][0]=None
                            hb_energies[i][j][0]=None
                        if possible_patterns[p_id][j][1] is None and hb_energies[i][j][1] > clear_threshold :
                            hb_trace[i][j][1]=None
                            hb_energies[i][j][1]=None
                if (direction[i]=='+' and p_id < 2) or (direction[i]=='-' and p_id >= 2) :
                    warnings+= '  in strand %d %s %s pattern %s maps on opposite direction! scores %s  ids %s\n' % (i,strand_seq[i],direction[i], str(pattern),str(scores),str(strand_ids[i]))
                    problematic_patterns+=[i]
    if doublecheck and complement_seq is not None and strange_strands is not None :
        hb_trace, hb_energies, warn=check_hb_trace_symmetry(hb_trace, hb_energies,direction,strand_ids,complementary_ids, strange_strands, strand_seq,complement_seq,hb_patterns,problematic_strands=problematic_patterns)
        warnings+=warn
    return hb_trace, hb_energies,warnings,problematic_patterns

def check_hb_trace_symmetry(hb_trace, hb_energies,direction,strand_ids,complementary_ids, strange_strands, strand_seq,complementary_seq,hb_patterns,problematic_strands=[],debug=True):
    '''
    elegant function most likely full of bugs.
     this fixes for instance  TFSLIIEAL facing HRQTALRSI - from file ['4XBM'] (it makes hb symmetric again)
    '''
    warn_str=''
    for j,stra_ids in enumerate(strand_ids) :
        if j in strange_strands or j in problematic_strands: 
            continue
# there should be a better way
        if direction[j]=='+' :off=1
        else : off=-1
        if stra_ids[::off] not in complementary_ids :
            continue
# SYMMETRY BREAKING (maybe solved by using the above criterion based on the ids in check_consistency_of_extracted_strand rather than that based on the sequence.
# However we may start to double count too much..
#if complementary_seq[j][::off]==strand_seq[j] :  (this would mean double counting.. i'm not sure whether the two following ones would mean the same)
# if we included the below we would have for example an homodimer with strand A and B inter chain so in the first chain you have -A--B- that map on -B--A- on the second
# if we use this complementary_ids mapping we would get that A,B counts twice, but in truth we should not consider two identical chains.
# On the other hand however we may have repeated protein, and we would count all occurrences as long as they are on one chain, and moreover strand that holds protein complexes
# together will more likely deliver optimal affinity.
# can happen if for dimers or such, in this case the sequence complementary to itself would not have been duplicated
# there are also multimer with weird symmetries such as 1a4e.dssp chain ['A'] (N and C terminal strands), better avoid raising the warning
#warn_str+=' *WARN* in check_hb_trace_symmetry the symmetrical sequence to strand at index %d cannot be found\n' % (j)
# another breaker is 'DSSP/1bi9.dssp', chain_to_grab=['A'] with last and third last parallel strands
        else :
            sym_ind= complementary_ids.index(stra_ids[::off])
            #if debug : print'  check_hb_trace_symmetry mathced',j,sym_ind
            while complementary_ids[j][::off]!=strand_ids[sym_ind] :
                try : 
                    sym_ind=complementary_ids.index(stra_ids[::off],sym_ind+1)
                    #if debug : print'  check_hb_trace_symmetry mathced',j,sym_ind
                except Exception :
                    sym_ind=None
                    break
            if sym_ind is None or sym_ind in problematic_strands :
                if debug : print('    check_hb_trace_symmetry skipping' ,j,sym_ind,complementary_ids.index(stra_ids[::off]))
                continue
            if direction[j]=='+' :
                for i,pattern in enumerate(hb_trace[j]) : # look for conflicting HB note that given the symmetry even if we don't update pattern here it will be updated with its symmetrical at some point (or before)
                    if pattern[0] is None and hb_trace[sym_ind][i][0]!=-1 and not (i==0 and hb_trace[sym_ind][i][0] is None): # at i=0 None would be a tolerated discrepancy
                        if i>0 and  hb_trace[j][i-1][1]==1 and (hb_energies[sym_ind][i][0] is None or hb_energies[sym_ind][i][0]>= hb_energies[j][i-1][1] ) : 
                            # provided that the other bond (hb_trace[j][i-1][1]) seem to confirm that pattern is right, if the energy of the potentially wrong bond is smaller than that of the potentially right one the wrong one is changed.  
                            if debug : print(' DEBcheck_hb_trace_symmetry + NtoO j=',j,'sym_ind=',sym_ind,'changing bond from',hb_trace[sym_ind][i][0],'e',hb_energies[sym_ind][i][0],'to -1 e',hb_energies[j][i-1][1])
                            hb_trace[sym_ind][i][0]=-1
                            hb_energies[sym_ind][i][0]=hb_energies[j][i-1][1]
                    elif pattern[0]==-1  and hb_trace[sym_ind][i][0] is not None :
                        if pattern[1]==1 and hb_trace[sym_ind][i][0] in [0,1] and (hb_energies[sym_ind][i][0] is None or hb_energies[sym_ind][i][0]>=(hb_energies[j][i][0]+hb_energies[j][i][1])/2. ) : # [0,1] is most likely wrong in a parallel strand
                            if debug : print(' DEBcheck_hb_trace_symmetry + None j=',j,'sym_ind=',sym_ind,'changing bond from',hb_trace[sym_ind][i][0],'e',hb_energies[sym_ind][i][0],'to None e None')
                            hb_trace[sym_ind][i][0]=None
                            hb_energies[sym_ind][i][0]=None
                            
                    if pattern[1] is None and hb_trace[sym_ind][i][1]!=1 and not (i==len(hb_trace[j])-1 and hb_trace[sym_ind][i][1] is None): # at i=end None would be a tolerated discrepancy
                        if i<len(hb_trace[j])-1 and  hb_trace[j][i+1][0]==-1 and (hb_energies[sym_ind][i][1] is None or hb_energies[sym_ind][i][1]>= hb_energies[j][i+1][0] ) : 
                            if debug : print(' DEBcheck_hb_trace_symmetry + OtoN j=',j,'sym_ind=',sym_ind,'changing bond from',hb_trace[sym_ind][i][1],'e',hb_energies[sym_ind][i][1],'to 1 e',hb_energies[j][i+1][0])
                            hb_trace[sym_ind][i][1]=1
                            hb_energies[sym_ind][i][1]=hb_energies[j][i+1][0]
                    elif pattern[1]==1  and hb_trace[sym_ind][i][1] is not None :
                        if pattern[0]==-1 and hb_trace[sym_ind][i][1] in [0,-1] and (hb_energies[sym_ind][i][1] is None or hb_energies[sym_ind][i][1]>=(hb_energies[j][i][0]+hb_energies[j][i][1])/2. ) : # [0,1] is most likely wrong in a parallel strand
                            if debug : print(' DEBcheck_hb_trace_symmetry + None j=',j,'sym_ind=',sym_ind,'changing bond from',hb_trace[sym_ind][i][1],'e',hb_energies[sym_ind][i][1],'to None e None')
                            hb_trace[sym_ind][i][1]=None
                            hb_energies[sym_ind][i][1]=None
            else : # -
                lo=len(hb_trace[j])-1
                p=['NtoO','OtoN']
                for i,pattern in enumerate(hb_trace[j]) : # look for conflicting HB note that given the symmetry even if we don't update pattern here it will be updated with its symmetrical at some point (or before)
                    for jj in [0,1] :
                        if pattern[jj] is None and hb_trace[sym_ind][lo-i][1-jj] is not None :
                            hbound_id=None
                            aa_id=strand_ids[sym_ind][lo-i]
                            if hb_trace[sym_ind][lo-i][1-jj]==-1 and lo-i>0 : # [-1,1] is probably wrong for a - strand
                                hbound_id= (complementary_ids[sym_ind][lo-i - 1],p[1-jj])
                            elif hb_trace[sym_ind][lo-i][1-jj]==1 and i>0 :
                                hbound_id= (complementary_ids[sym_ind][lo-i + 1],p[1-jj])
                            if hbound_id is not None and hbound_id[0] in hb_patterns : # if so check if that atom is actually forming another stronger HB
                                e=0
                                for pa_k in  hb_patterns[hbound_id[0]] :
                                    if pa_k[0]!=aa_id and pa_k[1]==p[jj] : e=hb_patterns[hbound_id[0]][pa_k]
                                
                                if hb_patterns[aa_id][hbound_id]>=e : # energy comparison, hb_patterns[aa_id][hbound_id] shoudl be = hb_energies[sym_ind][lo-i][1-jj]
                                    if debug : print(' DEBcheck_hb_trace_symmetry - None j=',j,'sym_ind=',sym_ind,'changing bond from',hb_trace[sym_ind][lo-i][1-jj],'e',hb_energies[sym_ind][lo-i][1-jj],'to None e None')
                                    hb_trace[sym_ind][lo-i][1-jj]=None
                                    hb_energies[sym_ind][lo-i][1-jj]=None
                        #elif pattern[jj]==0 and hb_trace[sym_ind][lo-i][1-jj]!=0 : NOT IMPLEMENTED
                                
    return hb_trace,hb_energies,warn_str
                            
                        
                    

def check_consistency_of_extracted_strand(strand_seq,complementary_seq,direction,strange_strands,strand_ids,complementary_ids,hb_patterns, perform_homomer_check=True,max_lendiff_homocorrection=20,create_missing=True, warn_id='',delete_strange_strands_shorter_than=4,fix_easy_strange=True,symmetrise_hb_pattern=True,debug=False):
    '''
    this double checks the results of beta_from_dssp().
    firt it checks if all strands have the same length as their complementary and complementary ids, it checks that no '?' are present among the ids (could happen if errors are generated).
      in case any of this fails the strand is removed.
    Then it checks whether all complementary also exist as reference sequences.
    Note that if you extract using a pdb chain list this might not be the case, 
      especially in the case of heterodimers, which could form extra-chain beta strands.
    if create_missing is False the ones that do not match the criterion are deleted
    perform_homomer_check = True is meaningful when extracting information about only one chain in a homomer. In this case 
     some of the complementary sequences might belong to other chains, and if the symmetry is perfect the data contain no mistakes, 
     however if the symmetry in the crystal structure is not perfect (e.g. 2hkn.dssp) one ends up with a pair of complementary sequences 
     (e.g. RVEVI<-->TGRHG and HRGT<-->VEVR) which should be recipracal. In this case the perform_homomer_check would leave RVEVI<-->TGRHG and GHRGT<-->IVEVR if homomer_keep_longer==True otherwise RVEV<-->TGRH and HRGT<-->VEVR
    fix_easy_strange=True will try to fix those strands with a '-' at the beginning or the end of either strand_seq or complementary_seq by shortening both and removing the strange label. 
       see for instance 1hii.dssp where strand_seq contains ['VEIEV', 'KKVRAT'] and complementary_seq ['ARVKK', 'VEIEV-'] the latter is in stange ids but without the dash and the corresponding T it would be the opposite of the first.
    '''
    # first check that the extracted strands make sense.
    delete=[]
    for j,ref_seq in enumerate(strand_seq) :
        if not len(ref_seq)==len(complementary_seq[j])==len(strand_ids[j])==len(complementary_ids[j]) : # this does not happen very often but in the case we delete the strands.
            delete+=[j]
        elif '?' in strand_ids[j] or '?' in complementary_ids[j] :
            delete+=[j]
    if delete!=[] :
        delete.sort(reverse=True)
        if debug : print('First quality check deletes',delete)
        for j in delete :
            del strand_ids[j],strand_seq[j],direction[j],complementary_seq[j],complementary_ids[j]
            if j in strange_strands : strange_strands.remove(j)
            for i,jj in enumerate(strange_strands) : # rescale other indeces (necessary even if this strand was not strange).
                if jj>j : strange_strands[i]=jj-1  
    if delete_strange_strands_shorter_than is not None :
        delete=[]
        for jd in strange_strands :
            if len(strand_seq[jd])<delete_strange_strands_shorter_than : delete+=[jd]
        if delete!=[]:
            if debug : print('delete_strange_strands_shorter_than quality check deletes',delete)
            for jd in sorted(delete,reverse=True) :
                del strand_ids[jd],strand_seq[jd],direction[jd],complementary_seq[jd],complementary_ids[jd]
                if jd in strange_strands : strange_strands.remove(jd)
                for idd,jjd in enumerate(strange_strands) : # rescale other indeces (necessary even if this strand was not strange).
                    if jjd>jd : strange_strands[idd]=jjd-1  
    # try to fix some easy strange strand
    if fix_easy_strange  : # try to fix this on the fly
        strange_to_remove=[]
        strands_to_remove=[] # those that will become too short if fixed (this situation almost never arise as the shorter strands woul also need to break the symmetry, e.g. comp_seq[:-1][::off] must be in strand_seq
        for j in strange_strands :
            if direction[j]=='-' : off=-1
            else : off=1
            comp_seq=complementary_seq[j]
            if comp_seq[0]=='-' and comp_seq[1:][::off] in strand_seq :
                if '-' not in comp_seq[1:] and comp_seq[1:].isupper() and '-' not in strand_seq[j][1:] and strand_seq[j][1:].isupper() : # this checks that after fixing it may no longer considered a strange_strand, otherwise no fixing is carried out.
                    complementary_seq[j]=comp_seq[1:]
                    strand_seq[j]=strand_seq[j][1:]
                    complementary_ids[j]=complementary_ids[j][1:]
                    strand_ids[j]=strand_ids[j][1:]
                    # we don't touch hb_patterns as these are a dictionary at this stage that may contain more info..
                    if len(strand_seq[j])<3 : strands_to_remove+=[j]
                    else :
                        if debug: print('   fix_easy_strange, removing',j,'we had comp_seq[0]==- dir=',direction[j]) 
                        strange_to_remove+=[j]# no need to rescale everything as the index in strange_strands are only labels and we have not removed strand j it is simply no longer strange
            if comp_seq[-1]=='-' and comp_seq[:-1][::off] in strand_seq :
                if '-' not in comp_seq[:-1] and comp_seq[:-1].isupper() and '-' not in strand_seq[j][:-1] and strand_seq[j][:-1].isupper() : # this checks that after fixing it may no longer considered a strange_strand, otherwise no fixing is carried out.
                    complementary_seq[j]=comp_seq[:-1]
                    strand_seq[j]=strand_seq[j][:-1]
                    complementary_ids[j]=complementary_ids[j][:-1]
                    strand_ids[j]=strand_ids[j][:-1]
                    # we don't touch hb_patterns as these are a dictionary at this stage that may contain more info..
                    if len(strand_seq[j])<3 : strands_to_remove+=[j]
                    else :
                        if debug: print('   fix_easy_strange, removing',j,'we had comp_seq[-1]==- dir=',direction[j])
                        strange_to_remove+=[j] # no need to rescale everything as the index in strange_strands are only labels and we have not removed strand j
            if strand_seq[j][0]=='-' and comp_seq[1:][::off] in strand_seq : # i'm not fully sure this situation may arise
                if '-' not in comp_seq[1:] and comp_seq[1:].isupper() and '-' not in strand_seq[j][1:] and strand_seq[j][1:].isupper() : # this checks that after fixing it may no longer considered a strange_strand, otherwise no fixing is carried out.
                    complementary_seq[j]=comp_seq[1:]
                    strand_seq[j]=strand_seq[j][1:]
                    complementary_ids[j]=complementary_ids[j][1:]
                    strand_ids[j]=strand_ids[j][1:]
                    # we don't touch hb_patterns as these are a dictionary at this stage that may contain more info..
                    if len(strand_seq[j])<3 : strands_to_remove+=[j]
                    else :
                        if debug: print('   fix_easy_strange, removing',j,'we had strand_seq[j][0]==- dir=',direction[j])
                        strange_to_remove+=[j]
            if strand_seq[j][-1]=='-' and comp_seq[:-1][::off] in strand_seq : # i'm not fully sure this situation may arise
                if '-' not in comp_seq[:-1] and comp_seq[:-1].isupper() and '-' not in strand_seq[j][:-1] and strand_seq[j][:-1].isupper() : # this checks that after fixing it may no longer considered a strange_strand, otherwise no fixing is carried out.
                    complementary_seq[j]=comp_seq[:-1]
                    strand_seq[j]=strand_seq[j][:-1]
                    complementary_ids[j]=complementary_ids[j][:-1]
                    strand_ids[j]=strand_ids[j][:-1]
                    # we don't touch hb_patterns as these are a dictionary at this stage that may contain more info..
                    if len(strand_seq[j])<3 : strands_to_remove+=[j]
                    else :
                        if debug: print('   fix_easy_strange, removing',j,'we had strand_seq[j][-1]==- dir=',direction[j])
                        strange_to_remove+=[j]
        for j in strange_to_remove : strange_strands.remove(j) # remove these first as then we will do -1 to some of these ids if strands_to_remove is not empty.
        if strands_to_remove!=[] :
            if debug : print('fix_easy_strange deleting as became too short:',strands_to_remove)
            for j in sorted(strands_to_remove,reverse=True) :
                del strand_ids[j],strand_seq[j],direction[j],complementary_seq[j],complementary_ids[j]
                if j in strange_strands : strange_strands.remove(j)
                for i,jj in enumerate(strange_strands) : # rescale other indeces (necessary even if this strand was not strange).
                    if jj>j : strange_strands[i]=jj-1  
    # now check that reciprocal pairs of complementary sequences are represented.
    wrong_ids=[] 
    '''
    # this is in principle a better way but at the end the perform_homomer_check would revert everythin like the old way..
    for j,comp_ids in enumerate(complementary_ids) :
        if direction[j]=='-' : off=-1
        else : off=1
        if comp_ids[::off] not in strand_ids :
            if complementary_seq[j][::off]==strand_seq[j] :continue # sequence is complementary to itself, obviously this can happen.
            wrong_ids+=[j]
    ''' # I SHOULD USE COMPLEMENTARY_INDS IN THAT WAY I WOULD BE SURE TO START WITH..
    for j,comp_seq in enumerate(complementary_seq) :
        if direction[j]=='-' : off=-1
        else : off=1
    
        if comp_seq[::off] not in strand_seq :
            wrong_ids+=[j]
        else :
            found=False
            jj=-1
            while not found : # check that they are actually complementary to each other
                try : jj=strand_seq.index(comp_seq[::off],jj+1)
                except Exception : break
                if complementary_seq[jj][::off]==strand_seq[j] : # in principle we should also check hb pattern but the chances are low and the dict format of hb_pattern at this stage makes it quite hard.. 
                    found=True
                    # MAybe check also HB pattern above
            if not found : wrong_ids+=[j] # this means that there is another sequence equal to comp_seq in strand_seq, but this is complementary to something else and it is not the reciprocal
    
    # now try to gradually solve the problems identified. 
    if debug: print('wrong_ids',wrong_ids)
    if wrong_ids!=[] :
        if perform_homomer_check :
            solved=[]
            for j in wrong_ids :
                if direction[j]=='-' :off=-1
                else : off=1
                comp=complementary_seq[j][::off]
                for jj in wrong_ids :
                    if direction[jj]!=direction[j] or j==jj : continue
                    # if the complementary is contained in the reference sequence we shorten the reference sequence and its related quantity.
                    if comp in strand_seq[jj] and (j not in solved or jj not in solved) and strand_seq[j] in complementary_seq[jj][::off] and abs(len(comp)-len(strand_seq[jj]))<max_lendiff_homocorrection : 
                        # Note that while this might not be ideal it is the only way with the current format of hb_patterns and other extracted variables.
                        # keeping the shorter relies on the idea that the other h bonds are an artifact of crystal conditions.
                        start = strand_seq[jj].find(comp)# comp is already N to C, no need to check for direction again
                        end=start+len(comp)
                        if debug:print('1',comp,start,end,j,jj)
                        strand_seq[jj]=strand_seq[jj][start:end]     
                        strand_ids[jj]=strand_ids[jj][start:end]     
                        complementary_seq[jj]=complementary_seq[jj][start:end]
                        complementary_ids[jj]=complementary_ids[jj][start:end]
                        if jj in strange_strands and '-' not in complementary_seq[jj] and complementary_seq[jj].isupper() and '-' not in strand_seq[jj] and strand_seq[jj].isupper() :
                            strange_strands.remove(jj)
                        solved+=[j,jj]
                    # if part of this reference sequence is contained in the complementary
                    elif strand_seq[jj] in comp and (j not in solved or jj not in solved) and complementary_seq[jj][::off] in strand_seq[j] and abs(len(comp)-len(strand_seq[jj]))<max_lendiff_homocorrection :     
                        if direction[j]=='-' : start = comp[::-1].find(strand_seq[jj][::-1])
                        else : start = comp.find(strand_seq[jj])
                        end=start+len(strand_seq[jj])
                        if debug:print('2',comp,start,end,j,jj)
                        strand_seq[j]=strand_seq[j][start:end]
                        strand_ids[j]=strand_ids[j][start:end]
                        complementary_seq[j]=complementary_seq[j][start:end]
                        complementary_ids[j]=complementary_ids[j][start:end]
                        if j in strange_strands and '-' not in complementary_seq[j] and complementary_seq[j].isupper() and '-' not in strand_seq[j] and strand_seq[j].isupper() :
                            strange_strands.remove(j)
                        solved+=[j,jj]
            if solved != [] :
                if debug : print('perform_homomer_check solved',solved)
                solved=misc.uniq(solved )
                solved.sort(reverse=True)
                for j in solved :
                    wrong_ids.remove(j)
        if wrong_ids!=[] : # add the missing one or remove the wrong ones...
            delete=[]
            if not create_missing : delete=wrong_ids[:]
            else :
                n_added=0
                for j in wrong_ids :
                    if direction[j]=='-' : off=-1
                    else : off=1
                    # first add hb patterns copying from the other strand, in case something goes wrong this gets deleted rather than duplicated
                    skip=False
                    for i,jj in enumerate( complementary_ids[j][::off] ) : # jj will be the id-s of the strand_seq we will add
                        if direction[j]=='-' : 
                            if strand_ids[j][::off][i] is None : continue # this happens sometimes when in strange 
                            t=strand_ids[j][::off][i]
                            k1,new_k1=(jj,'NtoO'),(t,'OtoN')
                            k2,new_k2=(jj,'OtoN'),(t,'NtoO')
                            if jj not in hb_patterns : hb_patterns[jj]={}
                            if t in hb_patterns and k1 in hb_patterns[ t ] :
                                if new_k1 in hb_patterns[jj] and hb_patterns[jj][new_k1]>=-0.5 : del hb_patterns[jj][new_k1] # very little energy
                                if new_k1 in hb_patterns[jj] and hb_patterns[jj][new_k1]!=hb_patterns[ t ][k1] and any([t in c for c in complementary_ids] ) : # the second condition checks that we have not removed it from the complementary pattern before, if we have then we overwrite rather than deleting the fragments
                                    if j not in strange_strands or debug : sys.stderr.write('**ERROR** in check_consistency_of_extracted_strand() option create_missing -, at wrong id %d, new compl id %d key %s already present. DELETING complementary strand at index %d %s\n' % (j,jj,new_k1,j,warn_id) )
                                    delete+=[j]
                                    skip=True
                                    break
                                hb_patterns[jj][new_k1]= hb_patterns[ t ][k1]
                                if debug : print(' adding h-bond - for new strand at index',len(strand_seq),'to id ',jj,'key1 ',new_k1,' from h-bond of id ',t,'key',k1,' taken from strand at index',j)
                            #else : # sometimes for dimers of which we read one chain only the above is not there but the pattern should be reconstructed anyway (this was the case when reading B of 1aw8 for ids about 7) 
                            if t in hb_patterns and k2 in hb_patterns[ t ] :
                                if new_k2 in hb_patterns[jj] and hb_patterns[jj][new_k2]>=-0.5 : del hb_patterns[jj][new_k2] # very little energy
                                if new_k2 in hb_patterns[jj] and hb_patterns[jj][new_k2]!=hb_patterns[ t ][k2] and any([t in c for c in complementary_ids] ) : 
                                    if j not in strange_strands or debug : sys.stderr.write( '**ERROR** in check_consistency_of_extracted_strand() option create_missing -, at wrong id %d, new compl id %d key %s already present. DELETING complementary strand at index %d %s\n' % (j,jj,new_k2,j,warn_id) )
                                    delete+=[j]
                                    skip=True
                                    break
                                hb_patterns[jj][new_k2]= hb_patterns[ t ][k2]
                                if debug : print(' adding h-bond - for new strand at index',len(strand_seq),'to id ',jj,'key2 ',new_k2,' from h-bond of id ',t,'key',k2,' taken from strand at index',j)
                        else : # parallel strand
                            # first HB
                            if i>0 and strand_ids[j][::off][i-1] is None : k1,new_k1=None,None
                            elif i==0 :
                                if strand_ids[j][::off][i] is not None :
                                    t1=strand_ids[j][::off][i]-1 # previous one even if outside strand
                                    k1,new_k1=(jj,'OtoN'),(t1,'NtoO')
                                else : k1,new_k1=None,None
                            else :
                                t1=strand_ids[j][::off][i-1]
                                k1,new_k1=(jj,'OtoN'),(t1,'NtoO')
                            if k1 is not None :
                                if jj not in hb_patterns : hb_patterns[jj]={}
                                if t1 in hb_patterns and k1 in hb_patterns[ t1 ] :
                                    if new_k1 in hb_patterns[jj] and hb_patterns[jj][new_k1]>=-0.3 : del hb_patterns[jj][new_k1] # very little energy
                                    if new_k1 in hb_patterns[jj] and hb_patterns[jj][new_k1]!=hb_patterns[ t1 ][k1] and any([t1 in c for c in complementary_ids] ) : # the second condition checks that we have not removed it from the complementary pattern before, if we have then we overwrite rather than deleting the fragments
                                        if j not in strange_strands or debug : sys.stderr.write('**ERROR** in check_consistency_of_extracted_strand() option create_missing +, at wrong id %d, new compl id %d key %s already present. DELETING complementary strand at index %d %s\n' % (j,jj,new_k1,j,warn_id) )
                                        delete+=[j]
                                        skip=True
                                        break
                                    hb_patterns[jj][new_k1]= hb_patterns[ t1 ][k1]
                                    if debug : print(' adding h-bond + for new strand at index',len(strand_seq),'to id ',jj,'key1 ',new_k1,' from h-bond of id ',t1,'key',k1,' taken from strand at index',j)
                            # second HB
                            if i<len(strand_ids[j][::off])-1 and strand_ids[j][::off][i+1] is None : k2,new_k2=None,None
                            elif i==len(strand_ids[j][::off])-1 :
                                if strand_ids[j][::off][i] is not None :
                                    t2=strand_ids[j][::off][i]+1 # following one even if outside strand
                                    k2,new_k2=(jj,'NtoO'),(t2,'OtoN')
                                else :  k2,new_k2=None,None
                            else : 
                                t2=strand_ids[j][::off][i+1]
                                k2,new_k2=(jj,'NtoO'),(t2,'OtoN')
                            if k2 is not None :
                                if jj not in hb_patterns : hb_patterns[jj]={}
                                if t2 in hb_patterns and k2 in hb_patterns[ t2 ] :
                                    if new_k2 in hb_patterns[jj] and hb_patterns[jj][new_k2]>=-0.3 : del hb_patterns[jj][new_k2] # very little energy
                                    if new_k2 in hb_patterns[jj] and hb_patterns[jj][new_k2]!=hb_patterns[ t2 ][k2] and any([t2 in c for c in complementary_ids] ) : # the second condition checks that we have not removed it from the complementary pattern before, if we have then we overwrite rather than deleting the fragments
                                        if j not in strange_strands or debug : sys.stderr.write('**ERROR** in check_consistency_of_extracted_strand() option create_missing +, at wrong id %d, new compl id %d key %s already present. DELETING complementary strand at index %d %s\n' % (j,jj,new_k2,j,warn_id) )
                                        delete+=[j]
                                        skip=True
                                        break
                                    hb_patterns[jj][new_k2]= hb_patterns[ t2 ][k2]
                                    if debug : print(' adding h-bond + for new strand at index',len(strand_seq),' to id ',jj,'key2 ',new_k2,' from h-bond of id ',t2,'key',k2,' taken from strand at index',j)
                            #print jj,i,strand_ids[j][::off][i],t1,k1,new_k1,t2,k2,new_k2
                    if skip : continue
                    strand_seq+=[ complementary_seq[j][::off] ]
                    complementary_seq+=[ strand_seq[j][::off] ]
                    complementary_ids+=[ strand_ids[j][::off] ]
                    strand_ids+=[ complementary_ids[j][::off] ]
                    direction+=[ direction[j] ]
                    if j in strange_strands : strange_strands+=[ len(direction)-1 ]
                    n_added+=1
                if debug and n_added>0: 
                    print('in check_consistency_of_extracted_strand() added %d strands to make non-symmetrical ones symmetric.' % (n_added))
            if delete!=[] :
                delete.sort(reverse=True)
                if debug : 
                    print('delete',delete)
                for j in delete :
                    del strand_ids[j],strand_seq[j],direction[j],complementary_seq[j],complementary_ids[j]
                    if j in strange_strands : strange_strands.remove(j)
                    for i,jj in enumerate(strange_strands) :
                        if jj>j : strange_strands[i]=jj-1      
    if symmetrise_hb_pattern :
        ''' prevent cases such as this were the dssp file is erroneusly annotated (from 3v7a.dssp chain E)
        >>> hb_patterns[510]
        {(513, 'OtoN'): -3.1, (514, 'OtoN'): -1.2, (515, 'NtoO'): -3.9}
        >>> hb_patterns[515]
        {(510, 'NtoO'): -1.2, (510, 'OtoN'): -3.9}
        the idea is that if a bond exist (in this case 515 to 510) then the opposite must exist as well (510 to 515) above it does exist but it appears as a 514
        , which is probably a small bug of dssp or former inaccuarcies in the pdb file.
        '''
        #energy_threshold=-1. no longer use this as next function in workflow will remove the non-predicted ones that are weaker than -1..
        for j,stra_ids in enumerate(strand_ids) : 
            for i,aa_id in enumerate(stra_ids) :
                if aa_id is None : continue
                if aa_id not in hb_patterns : continue
                for bond in hb_patterns[aa_id] :
                    #if hb_patterns[aa_id][bond]<=energy_threshold : # we care only for relatively strong bonds (those that may actually exist)
                        aa_id_bound,bond_direction=bond
                        if direction[j]=='-' and complementary_ids[j][i] is not None and aa_id_bound!= complementary_ids[j][i] : continue # not worth adding potentially missing HB that would not be relevant (e.g. those with Turns that have not been read in the first place
                        if direction[j]=='+' and complementary_ids[j][i] is not None  and (aa_id_bound!=complementary_ids[j][i]+1 or aa_id_bound!= complementary_ids[j][i]-1) : continue # not worth adding potentially missing HB that would not be relevant (e.g. those with Turns that have not been read in the first place
                        if bond_direction=='NtoO' : op_direction='OtoN'
                        elif bond_direction=='OtoN' : op_direction='NtoO'
                        else :
                            sys.stderr.write('**ERROR** while symmetrise_hb_pattern bond_direction %s not recognised aa_id=%d in %s\n' % (bond_direction,aa_id,warn_id))
                            continue
                        if aa_id_bound not in hb_patterns :
                            if debug : print('WARN symmetrise_hb_pattern adding bond %s energy=%lf in aa_id %d, which was not even present in the pattern dictionary ' % (str((aa_id,op_direction)),hb_patterns[aa_id][bond], aa_id_bound))
                            hb_patterns[aa_id_bound]={}
                            hb_patterns[aa_id_bound][(aa_id,op_direction)]=hb_patterns[aa_id][bond]
                        else :
                            found=False 
                            for op_bond in hb_patterns[aa_id_bound] :
                                if op_bond==(aa_id,op_direction) :
                                    found=True
                                    if hb_patterns[aa_id_bound][op_bond]!=hb_patterns[aa_id][bond] : # check that the energy match
                                        sys.stderr.write('**Warn** while symmetrise_hb_pattern  aa_id=%d bound to %d has different energy depending on where we look, %lf %lf in %s\n' % (bond_direction,aa_id,aa_id_bound,hb_patterns[aa_id_bound][op_bond],hb_patterns[aa_id][bond],warn_id))
                                        hb_patterns[aa_id_bound][op_bond]=hb_patterns[aa_id][bond] # we don't really care about this but we set them equal just in case the op_bond energy is lower than the threshold (which would mean it will become None when creating the patterns).
                                    break
                            if not found :
                                if debug : print('symmetrise_hb_pattern adding bond %s  energy=%lf in aa_id %d ' % (str((aa_id,op_direction)), hb_patterns[aa_id][bond], aa_id_bound))
                                hb_patterns[aa_id_bound][(aa_id,op_direction)]=hb_patterns[aa_id][bond]
    return strand_seq,complementary_seq,direction,strange_strands,strand_ids,complementary_ids,hb_patterns


class Beta_strands():
    '''
    This is an auxiliary class used for the DB creation/usage.
    Will be the content of a dictionary of classes. 
    the key will be the strand, (ref_seq) and the
    list comp_seq=[] will contain all possible complementary sequences
    the list comp_scores will contain the number of occurrences of each complementary sequence
    when finalize() method is called the complementaries list will contain (seq,scores),... sorted with the higher scores on top
      consider using free_others to free the two initial lists (if memory problems)
    check out also add_comp() function
    '''
    ref_seq=''
    comp_seq=[]
    hb_patterns=[] # list of hb_pattrens relative to ref_seq <--> comp_seq
    hb_energies=[] # list of hb_energies relative to ref_seq <--> comp_seq, after finalize contains average values, before sums
    comp_scores=[]
    pdb_ids={} #ids where the ref_seq and the comp_seq are found, key are (comp_seq, direction) or comp_seq alone
    complementary=None
    def __init__(self):
        self.ref_seq=''
        self.comp_seq=[]
        self.hb_patterns=[] # list of hb_pattrens relative to ref_seq <--> comp_seq
        self.hb_energies=[] # list of hb_energies relative to ref_seq <--> comp_seq, after finalize contains average values, before sums
        self.comp_scores=[]
        self.pdb_ids={}
        self.complementary=None
    def __repr__(self): #used when class is represented (i.e. when in python you press enter!)
        tmp_str=''
        for attr, value in self.__dict__.items():
            if type(value) is list or type(value) is dict and len(value)<10 : tmp_str+=str(attr)+'='+repr(value)+' '
            elif type(value) is list or type(value) is dict : tmp_str+=str(attr)+'='+str(type(value))+' '
            else : tmp_str+=str(attr)+'='+repr(value)+' '
        return tmp_str[:-1]
    def __str__(self): #used when class is printed
        tmp_str=''
        for attr, value in self.__dict__.items():
            if type(value) is list or type(value) is dict and len(value)<10 : tmp_str+=str(attr)+'='+repr(value)+' '
            elif type(value) is list or type(value) is dict : tmp_str+=str(attr)+'='+str(type(value))+' '
            else : tmp_str+=str(attr)+'='+repr(value)+' '
        return tmp_str[:-1]
    
    def retrieve( self, variable_name_str ) :
        '''
        retrieve the value of a variable by calling it as a string. One can do residue.retrieve('variable')
        so one can use if 'variable' in dir(residues) : and staff like that...
        '''
        if variable_name_str not in dir(self) :
            sys.stderr.write('**ERROR** in method retrieve, %s not defined in class\n    %s\n' % (variable_name_str,str(dir(self))) )
            raise AttributeError 
        return self.__dict__[variable_name_str]
    
    def update( self, variable_name_str, new_variable_value, check=False  ) :
        '''
        used to change the value of a variable by calling it as a str (ranther than residue.variable=new_value we do residue.update('variable',new_value)) 
        so one can use if 'variable' in dir(residues) : and staff like that...
        '''
        if check and variable_name_str not in dir(self) :
            sys.stderr.write('**ERROR** in method retrieve, %s not defined in class\n    %s\n' % (variable_name_str,str(dir(self))) )
            raise AttributeError 
        self.__dict__[variable_name_str] = new_variable_value
        return
    
    def add_comp(self,new_comp_str, direction=None, pdb_id=None, hb_pattern=None, hb_energy=None):
        '''
        this adds another complementary strand,
         it updates the list comp_seq in the correct way, together with the comp_score field
        '''
        if hb_pattern== [ ['?','?'] ] : # don't add the complementary 
            sys.stderr.write(" **WARN** in add_comp reference %s comp strand %s %s from pdb %s rejected because pattern ??\n" % (self.ref_seq,new_comp_str,str(direction),str(pdb_id)))
            return
        add_it=True
        if direction is not None :
            if (new_comp_str,direction) in self.comp_seq : # update score, the complementary is already there
                if hb_pattern is not None :
                    add_it=False
                    index=self.comp_seq.index((new_comp_str,direction))
                    while True :
                        if hb_pattern== self.hb_patterns[ index ] :
                            self.comp_scores[ index ]+=1 # the pattern is the same, add 1 to the count score.
                            if hb_energy is not None :
                                for i in range(0,len(self.hb_energies[ index ])) :
                                    if self.hb_energies[index][i][0]  is not None : self.hb_energies[index][i][0]+=hb_energy[i][0] # sum the energy to the existing one, with finalize() we will take the average.
                                    if self.hb_energies[index][i][1]  is not None : self.hb_energies[index][i][1]+=hb_energy[i][1]
                            break
                        else :
                            try :
                                index=self.comp_seq.index((new_comp_str,direction),index+1) # look in the next part
                            except ValueError :
                                add_it=True
                                break
                else : # in this case just update the score
                    self.comp_scores[ self.comp_seq.index((new_comp_str,direction)) ]+=1
            if add_it : # add a new complementary to the list
                self.comp_seq += [ (new_comp_str,direction) ]
                self.comp_scores += [ 1 ]
                if hb_pattern is not None :
                    self.hb_patterns+= [ hb_pattern ]
                    if hb_energy is not None :
                        self.hb_energies += [ hb_energy ]
            if pdb_id is not None :
                if (new_comp_str,direction) not in self.pdb_ids : self.pdb_ids[(new_comp_str,direction)] = [pdb_id]
                elif pdb_id not in self.pdb_ids[(new_comp_str,direction)] : 
                    self.pdb_ids[(new_comp_str,direction)] += [pdb_id]
        else : # very basic without even strand kind.. (aka useless)
            if new_comp_str in self.comp_seq : # update score
                self.comp_scores[ self.comp_seq.index(new_comp_str) ]+=1
            else : # add it to list
                self.comp_seq += [ new_comp_str ]
                self.comp_scores += [ 1 ]
            if pdb_id is not None :
                if new_comp_str not in self.pdb_ids : self.pdb_ids[new_comp_str] = [pdb_id]
                elif pdb_id not in self.pdb_ids[new_comp_str] : self.pdb_ids[new_comp_str] += [pdb_id]
        return
    
    def finalize(self,free_others=False) :
        '''
        use after filling all to sort the comp_seq list (and consequently the comp_scores one) so that highest scores are at the beginning
        ''' 
        if self.hb_energies!=[] : # we do the average energies by dividing for the scores
            for i,en_strand in enumerate(self.hb_energies) :
                for j in range(0, len(en_strand)) :
                    if self.hb_energies[i][j][0]  is not None and self.hb_energies[i][j][0] !=0. : self.hb_energies[i][j][0] /=(1.*self.comp_scores[i])
                    if self.hb_energies[i][j][1]  is not None and self.hb_energies[i][j][1] !=0. : self.hb_energies[i][j][1] /=(1.*self.comp_scores[i])
        if len(self.comp_scores)==len(self.comp_seq)==len(self.hb_energies)==len(self.hb_patterns) :
            self.comp_scores,self.comp_seq,self.hb_energies,self.hb_patterns=list(zip(*sorted( zip(self.comp_scores,self.comp_seq,self.hb_energies,self.hb_patterns), reverse=True )))
        elif len(self.comp_scores)==len(self.comp_seq)==len(self.hb_patterns) :
            self.comp_scores,self.comp_seq,self.hb_patterns=list(zip(*sorted( zip(self.comp_scores,self.comp_seq,self.hb_patterns), reverse=True )))
        elif len(self.comp_scores)==len(self.comp_seq) and self.hb_patterns==[] :
            self.comp_scores,self.comp_seq=list(zip(*sorted( zip(self.comp_scores,self.comp_seq), reverse=True )))
        else :
            sys.stderr.write("**ERROR** in Beta_strand.finalize() no sorting performed as length are uncompatible for ref_seq %s\n\n" % (self.ref_seq))
        self.complementary=list(zip(self.comp_seq,self.comp_scores))
        if free_others :
            del self.comp_scores
            del self.comp_seq
        return
            
            
            
'''
**ERROR** in check_consistency_of_extracted_strand() option create_missing -, at wrong id 27, new compl id 567 key (539, 'OtoN') already present. DELETING complementary strand at index 27 1b70_['B', 'A']
**ERROR** in check_consistency_of_extracted_strand() option create_missing -, at wrong id 25, new compl id 303 key (324, 'OtoN') already present. DELETING complementary strand at index 25 1cjx_['A']

ERROR while matching hb patterns, found complementary start resnumbers (84) - of 4,GYYIS [76, 77, 78, 79, 80] complementary: FFT-I but whole complment doesnt map on sequence in file DSSP/1blj.dssp ch=['A']
'''
   




def analyze_dssp_folder(input_dir, chain_list_file=None,  only_longer_than=2,get_aa_abundance=False,debug=False, out_to_pickle=True,discard_strands_with_problematic_hb=True,save_fasta_of_processed_sequence=True,delete_strange_strands_shorter_than=4, warnfile_name='stderr'):
    '''
    this analyzes a whole folder of DSSP files, it is the key part of the BS database creation.
    it pickles output using protocol 0 in strand_dict.pkl and strange_strand_dict.pkl, which are the BS databses.
    '''
    if chain_list_file!='' and chain_list_file is not None : use_chain_list=True
    else : use_chain_list= False
    if input_dir[-1]!='/' : input_dir +='/'
    file_list=os.listdir(input_dir)
    N=len(file_list)
    sys.stdout.write('%d files found in folder %s\n\n' % (N,input_dir) )
    if warnfile_name=='stderr' : warnings=sys.stderr
    elif warnfile_name=='stdout' : warnings=sys.stdout
    else :
        warnings=open(warnfile_name,'w')
    warnings.write('warning file from analyze_dssp_folder(%s)\n\n' %(input_dir))
    if use_chain_list : # use list with pdbname and the chains we are interested in
        chain_list={}
        if type(chain_list_file) is str and os.path.isfile(chain_list_file): pass
        elif os.path.isfile(input_dir+chain_list_file) : chain_list_file=input_dir+chain_list_file
        chain_file=file(chain_list_file).read().splitlines()
        for line in chain_file :
            if line[0]!='#' and line[0]!=' ' and line[0]!='\n' :
                try :
                    chain_list[ line.split()[0].split('.')[0] ] = line.split()[1].split(';')[:-1] #save so that pdb_id is the key and the entry is a list of chain ids
                except IndexError :
                    sys.stderr.write('**Index error in chain list %s while processing line :\n\t|%s|\n' % (chain_list_file,line))
                    pass
        sys.stderr.write('%d files in chain_list\n\n' % len(chain_list) )
        N=len(chain_list)
    compute_abundance=False
    if get_aa_abundance==True or get_aa_abundance=={} : 
        get_aa_abundance={}
        for res in amino_list1 :
            get_aa_abundance[res]=0
        compute_abundance=True
    elif type(get_aa_abundance) is dict : compute_abundance=True
    if save_fasta_of_processed_sequence :
        outfasta=open('whole_processed_sequence.fasta','w')
    strand_dict={} # the final output, a dictionary of classes Beta_strands()
    strange_strand_dict={} # the final output, a dictionary of classes Beta_strands() (this contains the strange pairs, the one with extra non-pairing residues)
    run_on_chain=[ None ] # in this case it will run on everything unless a chain_list_file is given
    processed_pdb_id=[]
    n=0
    i=0
    skip_exception=0
    no_strand=0
    just_another_file_format=0
    not_in_pdblist=0
    for dssp_file in file_list :
        if '.swp' in dssp_file or '.dssp' not in dssp_file or dssp_file[0]=='.': 
            just_another_file_format+=1
            continue
        _,pdb_id,_=misc.get_file_path_and_extension(dssp_file)
        if use_chain_list :
            if pdb_id.upper() not in chain_list :
                not_in_pdblist+=1 
                continue
            run_on_chain=chain_list[pdb_id.upper()]
        if i%40==0 :
            sys.stdout.write('%4d ->approx %5.2lf %% |' % (i,100.*i/(1.*N)) )
            sys.stdout.flush()
        i+=1
        try :
            if run_on_chain==[None]: ch_id=None
            else : ch_id=run_on_chain
            #for ch_id in run_on_chain : # note that if ch_id=None get_HX_from_pdb will run on all the chains
            if debug :
                warnings.flush()
                sys.stdout.write( '\n  %d file %s chain %s ...' % (i,dssp_file,str(ch_id)) )
                sys.stdout.flush()
            if compute_abundance : 
                strand_seq,complementary_seq,direction,_,sequence,strange_strands,strand_ids,complementary_ids,hb_patterns,get_aa_abundance,processed_sequence = beta_from_dssp(input_dir+dssp_file, chain_to_grab=ch_id, only_longer_than=only_longer_than , warn_file=warnings, get_aa_abundance=get_aa_abundance)
            else :
                strand_seq,complementary_seq,direction,_,sequence,strange_strands,strand_ids,complementary_ids,hb_patterns,processed_sequence=beta_from_dssp(input_dir+dssp_file, chain_to_grab=ch_id, only_longer_than=only_longer_than , warn_file=warnings)
            if strand_seq==[] :
                no_strand+=1
                continue
            # processed_sequence is the sequence of the chain requested, which with good approximation is the whole sequence we processed (in truth we may have processed sequence from extra chain beta strands but that will be a striking minority)
            
            strand_seq, complementary_seq, direction, strange_strands, strand_ids, complementary_ids,hb_patterns \
              = check_consistency_of_extracted_strand(strand_seq,complementary_seq,direction,strange_strands,strand_ids,complementary_ids,hb_patterns,delete_strange_strands_shorter_than=delete_strange_strands_shorter_than, perform_homomer_check=True,create_missing=True, warn_id=pdb_id+'_'+str(ch_id))
            if strand_seq==[] : # we may have removed all in check_consistency_of_extracted_strand
                no_strand+=1
                continue
            hb_trace,hb_energy, pat_problems,problematic_patterns=process_hb_pattern(strand_seq,sequence,strand_ids,complementary_ids,hb_patterns,direction, strange_strands=strange_strands, complement_seq=complementary_seq)
            if pat_problems!='' :
                warnings.write("**WARNING** when analyzing hb patterns file %s chain %s\n%s" %  (dssp_file,str(ch_id),pat_problems))
                if discard_strands_with_problematic_hb :
# NB This may create asymmetries in the patterns (you may solve this by calling again the check_consistency_of_extracted_strand(), but at that point one needs to be careful or re-fixing hb_trace,hb_energy 
                    warnings.write("---> Discarding the corresponding strands! (%s)\n" % (str(problematic_patterns)))
                    for ij in sorted(problematic_patterns,reverse=True) :
                        del strand_ids[ij],strand_seq[ij],direction[ij],complementary_seq[ij],complementary_ids[ij]
                        if ij in strange_strands : strange_strands.remove(ij)
                        for ii,jj in enumerate(strange_strands) : # rescale other indeces (necessary even if this strand was not strange).
                            if jj>ij : strange_strands[ii]=jj-1  
            if debug :
                warnings.flush()
                sys.stdout.write('  found %d strands, now adding to strand_dict\n' % (len(strand_seq)))
                sys.stdout.flush()
            if save_fasta_of_processed_sequence :
                outfasta.write('> %s\n%s\n' % (pdb_id.lower()+'_'+''.join(ch_id),processed_sequence))
            # now fill the strand dictionaries
            for j,ref_seq in enumerate(strand_seq) :
                if j not in strange_strands : #avoid stange strands here
                    if ref_seq not in strand_dict : #add a new entry
                        n+=1
                        strand_dict[ref_seq]=Beta_strands() # init
                        strand_dict[ref_seq].ref_seq=ref_seq # save reference sequence
                    strand_dict[ref_seq].add_comp( complementary_seq[j],direction=direction[j],pdb_id=pdb_id.upper(),hb_pattern=hb_trace[j],hb_energy=hb_energy[j] ) # this updates the count score
                else :
                    if ref_seq not in strange_strand_dict :
                        strange_strand_dict[ref_seq]=Beta_strands() # init
                        strange_strand_dict[ref_seq].ref_seq=ref_seq # save reference sequence
                    strange_strand_dict[ref_seq].add_comp( complementary_seq[j],direction=direction[j],pdb_id=pdb_id.upper(),hb_pattern=hb_trace[j],hb_energy=hb_energy[j] )
            processed_pdb_id+=[ pdb_id.upper() ]
        
        except Exception :
            skip_exception+=1
            warnings.write( '\n**ERROR** while processing %s ch=%s\n File skipped\n' % (dssp_file,str(ch_id)))
            print_exc(file=warnings)
            warnings.flush()
            continue
    
    if save_fasta_of_processed_sequence : outfasta.close()
    # finalize and sort by score
    for entry in strand_dict : strand_dict[entry].finalize()
    for entry in strange_strand_dict : strange_strand_dict[entry].finalize()
    #pickle output using protocol 0.
    if out_to_pickle :
        output=open('strand_dict.pkl','wb')
        pickle.dump(strand_dict, output)
        output.close()
        output=open('strange_strand_dict.pkl','wb')
        pickle.dump(strange_strand_dict, output)
        output.close()
        if compute_abundance :
            output=open('aa_abundance.pkl','wb')
            pickle.dump(get_aa_abundance, output)
            output.close()
    sys.stdout.write('\n\n%d entries in pickle, corresponding to %d files processed\n' %(n,i))
    if len(processed_pdb_id)!=N :
        sys.stdout.write('%d file correctly analyzed among %d expected (%d in folder/chain_list but %d of different format)...\n  %d skipped because of exception, %d skipped because not containing beta strands. (not_in_pdblist = %d)\n' % (len(processed_pdb_id),N-just_another_file_format,N,just_another_file_format,skip_exception,no_strand,not_in_pdblist))
        if use_chain_list :
            sys.stdout.write('  skipped pdb ids are:\n')
            for pdb_id in chain_list :
                if pdb_id not in processed_pdb_id :
                    sys.stdout.write('%s ' % pdb_id)
            sys.stdout.write('\n\n')
    if warnfile_name!='stderr' and warnfile_name!='stdout' : warnings.close()
    if compute_abundance : return strand_dict,strange_strand_dict,get_aa_abundance
    return strand_dict,strange_strand_dict



def hb_pattern_to_tuple(hb_pattern):
    '''
    Auxiliary function used to make the pattern hashable.
    '''
    if type(hb_pattern) is tuple : return hb_pattern
    for i,p in enumerate(hb_pattern) :
        if type(p) is list :
            hb_pattern[i]=tuple(p)
    return tuple(hb_pattern)





class Comp():
    '''
    this stores information of each comlemenatary fragment to a given reference fragment. It is the values of the BSn databases.
    a list of these classes will be the values of the dictionary of complementary fragments contained in the Complement class (e.g. C.antiparallel_fragments[2] is a dict whose keys are complementary fragments of length 3)
    '''
    def __init__(self, complementary_sequence='',hb_pattern=None,count_score=-1,promiscuity_score=None,unique_reference_id=None,hb_energy=None):
        self.complementary_sequence=complementary_sequence
        self.hb_pattern=hb_pattern
        self.count_score=count_score
        self.promiscuity_score=promiscuity_score
        #self.strand_dict_reference_id=unique_reference_id # not implemented
        #self.hb_energy=hb_energy # not implemented 
    def __repr__(self):
        '''
        return '|Comp_seq %s count=%lf prom=%lf hb=%s|' % (self.complementary_sequence,float(self.count_score),str(self.hb_pattern))
        '''
        return '|Comp_seq %s count=%s prom=%s hb=%s|' % (self.complementary_sequence,repr(self.count_score),str(self.promiscuity_score),str(self.hb_pattern))
    def __str__(self):
        '''
        return '|Comp_seq %s count=%lf prom=%lf hb=%s|' % (self.complementary_sequence,float(self.count_score),str(self.hb_pattern))
        '''
        return '|Comp_seq %s count=%s prom=%s hb=%s|' % (self.complementary_sequence,repr(self.count_score),str(self.promiscuity_score),str(self.hb_pattern))
    def __len__(self):
        return len(self.complementary_sequence)
    def __getitem__(self, y) :
        '''
        we redefine this to mantain back compatibility as at the beginning instead of the Comp classes we had tuples in the 
        form (('IDA', ((None, None), (0, 0), (None, None))), 2)
        '''    
        if y==0 : return ( self.complementary_sequence, self.hb_pattern )
        if y==1 : return self.count_score
        if y==-1: return self.count_score
        raise IndexError

def parse(BSn_key_value_list, complementary_sequence, hb_pattern=None):
    rt=[]
    for comp in BSn_key_value_list :
        if comp.complementary_sequence==complementary_sequence :
            if hb_pattern is None : rt+=[comp]
            elif comp.hb_pattern==hb_pattern : rt+=[comp]
    if len(rt)==1 : return rt[0]
    else : return rt

def fragments_from_database(strand_dict_or_pickle, fragment_length=3, only_actual_n_plets=False, look_for_hb_pattern=True, strange_strand_dic_or_pickle=None, out_pickle_name=None, try_to_retrieve_lower_opposite_ch=False,add_promiscuity_score=True,calculate_enrichment_tables=True,debug=False):
    '''
    from the strand dict it generates the BSn databases, aka the fragment dictionaries of fixed length (one parallel and one antiparallel) that we use to create complementary pepetides.
     if only_actual_n_plets is enabled only the strands that are actually long fragment_length are considered. Otherwise strands longer than that are fragmented and accounted for (default)
      try_to_retrieve_lower_opposite_ch is an obscure option used to parase the strange_strand_dic_or_pickle when given. Since this contains the non-canonical beta-strands or the exceptions the
        try_to_retrieve_lower_opposite_ch==True will try to rescue some of them.
    '''
    if type(strand_dict_or_pickle) is not dict :
        strand_dict_or_pickle=pickle.load( open( strand_dict_or_pickle, "rb" ) ) # now it's a dict
    # generate all possible fragments
    frag_dic={} # this is a dict whose keys are the fragments, 
                # while being filled entries are dict as well, each has the complementary fragment (or (fragment, hb_pattern) as a key and its score (or [score, hb_energies]) as an entry
                # after entries will be list of sorted tuples that later on will become Comp classes. Each tuple will have the complementary fragments and its score (e.g. ('VTLV', 2) )
                # the sorting is such that higher score fragments will be first
    frag_dic_anti={} # same as above but for the anti_parallel beta strands (if info about direction is provided). 
    for strand_seq in strand_dict_or_pickle : # loop on the reference sequence in the dict database
        if len(strand_seq)>=fragment_length : # check that the reference sequence is at least as long as fragment_length
            if only_actual_n_plets and len(strand_seq)!=fragment_length :
                continue
            for j,comp_tuple in enumerate(strand_dict_or_pickle[strand_seq].complementary) : # loop on all the possible complementary sequences for the given sequence
                for i in range(0,len(strand_seq)-fragment_length+1) : # create window of size fragment_length that slides on the sequence
                    ref_seq=strand_seq[i:i+fragment_length] # this is the reference sequence.
                    if type(comp_tuple[0]) is tuple : #this means that we know whether it is parallel or antiparallel
                        if comp_tuple[0][1]=='+' : # it is parallel
                            if ref_seq not in frag_dic : 
                                frag_dic[ref_seq]={} # initialize as empty dictionary
                            if look_for_hb_pattern :
                                # define key as (comp_sequence, hb_pattern)
                                key = ( comp_tuple[0][0][i:i+fragment_length] , hb_pattern_to_tuple(strand_dict_or_pickle[strand_seq].hb_patterns[j][i:i+fragment_length]) )
                                if key[1]==() :
                                    print('*** skipping pair %s %s +  ==> missing hb_pattern' % (ref_seq,key[0])) 
                                    continue
                                if key in frag_dic[ref_seq] : 
                                    frag_dic[ref_seq][key] += comp_tuple[1] # update score
                                else :
                                    frag_dic[ref_seq][key] = comp_tuple[1] # save score
                            # Disused, we always have hb_pattern
                            else :   
                                if comp_tuple[0][0][i:i+fragment_length] in frag_dic[ref_seq] : 
                                    frag_dic[ref_seq][comp_tuple[0][0][i:i+fragment_length]] += comp_tuple[1] # update score
                                else :
                                    frag_dic[ref_seq][comp_tuple[0][0][i:i+fragment_length]] = comp_tuple[1] # save score 
                        elif comp_tuple[0][1]=='-' : # it is antiparallel
                            if ref_seq not in frag_dic_anti : 
                                frag_dic_anti[ref_seq]={} # initialize as empty dictionary 
                            if look_for_hb_pattern :
                                key = ( comp_tuple[0][0][i:i+fragment_length] , hb_pattern_to_tuple(strand_dict_or_pickle[strand_seq].hb_patterns[j][i:i+fragment_length]) )
                                if key[1]==() :
                                    print('*** skipping pair %s %s -  ==> missing hb_pattern' % (ref_seq,key[0])) 
                                    continue
                                if key in frag_dic_anti[ref_seq] :
                                    frag_dic_anti[ref_seq][key] += comp_tuple[1] # update score
                                else :
                                    frag_dic_anti[ref_seq][key] = comp_tuple[1] # save score
                            else : 
                                if comp_tuple[0][0][i:i+fragment_length] in frag_dic_anti[ref_seq] :
                                    frag_dic_anti[ref_seq][comp_tuple[0][0][i:i+fragment_length]] += comp_tuple[1] # update score
                                else :
                                    frag_dic_anti[ref_seq][comp_tuple[0][0][i:i+fragment_length]] = comp_tuple[1] # save score 
                        else :
                            sys.stderr.write("WARNING in fragments_from_database() strand %s doesn't have direction (%s) pdb_ids=%s\n" % (strand_seq,comp_tuple[0][1],str(strand_dict_or_pickle[strand_seq].pdb_ids)))  
                    # Disused bit:
                    else : # here we don't know if it is parallel or anti-parallel, hence we also ignore possible hb_patterns
                        sys.stderr.write("WARNING in fragments_from_database() strand %s doesn't have direction (%s) pdb_ids=%s\n" % (strand_seq,comp_tuple[0][1],str(strand_dict_or_pickle[strand_seq].pdb_ids)))
                        if ref_seq not in frag_dic : 
                            frag_dic[ref_seq]={} # initialize as empty dictionary
                        if comp_tuple[0][i:i+fragment_length] in frag_dic[ref_seq] : 
                            frag_dic[ref_seq][comp_tuple[0][i:i+fragment_length]] += comp_tuple[1] # update score
                        else :
                            frag_dic[ref_seq][comp_tuple[0][i:i+fragment_length]] = comp_tuple[1] # save score
    if strange_strand_dic_or_pickle is not None:
        del strand_seq
        if type(strange_strand_dic_or_pickle) is not dict :
            strange_strand_dic_or_pickle=pickle.load( open( strange_strand_dic_or_pickle, "rb" ) ) # now it's a dict
        
        for strand_seq in strange_strand_dic_or_pickle : # loop on the reference sequence in the dict database
            if len(strand_seq)>=fragment_length : # check tha the reference sequere is at least as long as fragment_length
                key_strand_seq=strand_seq
                for j,comp_tuple in enumerate(strange_strand_dic_or_pickle[key_strand_seq].complementary) : # loop on all the possible complementary sequences for the given sequence
                    if type(comp_tuple[0]) is tuple : # if we have directional information
                        comp_seq=comp_tuple[0][0]
                        dir_info=comp_tuple[0][1]
                        hb_p=''
                        if look_for_hb_pattern : 
                            hb_p = strange_strand_dic_or_pickle[key_strand_seq].hb_patterns[j]
                            if hb_p==() or hb_p is None or hb_p==[] : continue
                    else :
                        comp_seq=comp_tuple[0]
                        dir_info=None
                    if strand_seq.upper()!=strand_seq and comp_seq.upper()!=comp_seq: #check for lower case letter, if in both strand and complement they should be amino acids in strand that do not pair, keep them if they are of opposite charge
                        if try_to_retrieve_lower_opposite_ch :

                            ###print '  bef',strand_seq,comp_seq,hb_p
                            for k,res in enumerate(strand_seq) :
                                if res.islower() and res.upper() in amino_list1:
                                    if comp_seq[k].islower() and comp_seq[k].upper() in amino_list1:
                                        if mybio.aa_charge[comp_seq[k].upper()]*mybio.aa_charge[res.upper()] < 0 : # if they have opposite charge
                                            if k+1<len(strand_seq) : strand_seq=strand_seq[:k]+strand_seq[k].upper()+strand_seq[k+1:]
                                            else : strand_seq=strand_seq[:k]+strand_seq[k].upper()
                                            if k+1<len(comp_seq) : comp_seq=comp_seq[:k]+comp_seq[k].upper()+comp_seq[k+1:]
                                            else : comp_seq=comp_seq[:k]+comp_seq[k].upper() 
                                            ###print '  aft',k,strand_seq,comp_seq
                                        
                    strand_bits, comp_bits, bits = mybio.split_string_at(strand_seq, amino_list1, second_string=comp_seq, only_longer_than=fragment_length-1)
                    if look_for_hb_pattern : # in this case check that they are correctly annotate, otherwise delete/skip them.
                        try :
                            hb_p_bits,_= mybio.split_at_bits(bits, hb_p, only_longer_than=fragment_length-1)
                            if len(hb_p_bits)!=len(comp_bits) :
                                sys.stderr.write("*** Error in strange_strand_dic while splitting hb pattern for strand %s and its complementary strand %s,\n   pattern %s split to %s SKIPPING [Error len(hb_p_bits)!=len(comp_bits) %d %d ]\n" % (key_strand_seq,comp_seq,str(hb_p),str(hb_p_bits),len(hb_p_bits),len(comp_bits)))
                                continue 
                            else :
                                remove=[]
                                for jj,comp in enumerate(comp_bits) :
                                    if len(hb_p_bits[jj])!=len(comp) :
                                        sys.stderr.write("*** Warn in strange_strand_dic while splitting hb pattern for strand %s and its complementary strand %s, will be removed!\n   pattern %s split to %s [Error len(hb_p_bits[jj])!=len(comp) %d %d ]\n" % (key_strand_seq,comp_seq,str(hb_p),str(hb_p_bits),len(hb_p_bits[jj]),len(comp)))
                                        remove+=[jj]
                                if remove!=[] :
                                    remove.sort(reverse=True)
                                    for jj in remove :
                                        del hb_p_bits[jj],comp_bits[jj],strand_bits[jj]
                                    sys.stderr.write("*** Warn removing %d possible candidates that don't match hb_pattern qualitiy in strange_strand_dic while splitting hb pattern for strand %s and its complementary strand %s,\n   pattern %s split to %s \n" % (len(remove),key_strand_seq,comp_seq,str(hb_p),str(hb_p_bits)))
                                    if len(hb_p_bits)==0 : continue # if you have removed them all
                        except Exception :
                            sys.stderr.write("*** Error in strange_strand_dic while splitting hb pattern for strand %s and its complementary strand %s,\n   pattern %s split to %s \n" % (key_strand_seq,comp_seq,str(hb_p),str(hb_p_bits)))
                            print_exc(file=sys.stderr)
                            continue
                    
                    for k,strand_bit in enumerate(strand_bits) :
                        if only_actual_n_plets and len(strand_bit)!=fragment_length :
                            continue
                        for i in range(0,len(strand_bit)-fragment_length+1) : # create window of size fragment_length that slides on the sequence
                            ref_seq=strand_bit[i:i+fragment_length]
                            if dir_info=='+' or dir_info is None : # it's parallel
                                if ref_seq not in frag_dic : 
                                    frag_dic[ref_seq]={} # initialize as empty dictionary
                                if look_for_hb_pattern :
                                    key= ( comp_bits[k][i:i+fragment_length] , hb_pattern_to_tuple(hb_p_bits[k][i:i+fragment_length]) )
                                    if key in frag_dic[ref_seq] : # recall that frag_dic has already been filled with all possible fragments
                                        frag_dic[ref_seq][key] += comp_tuple[1] # update score
                                    else :
                                        frag_dic[ref_seq][key] = comp_tuple[1] # save score
                                else :
                                    if comp_bits[k][i:i+fragment_length] in frag_dic[ref_seq] : # recall that frag_dic has already been filled with all possible fragments
                                        frag_dic[ref_seq][comp_bits[k][i:i+fragment_length]] += comp_tuple[1] # update score
                                    else :
                                        frag_dic[ref_seq][comp_bits[k][i:i+fragment_length]] = comp_tuple[1] # save score
                            elif dir_info=='-' : # it's anti parallel
                                if ref_seq not in frag_dic_anti : 
                                    frag_dic_anti[ref_seq]={} # initialize as empty dictionary
                                if look_for_hb_pattern :
                                    key= ( comp_bits[k][i:i+fragment_length] , hb_pattern_to_tuple(hb_p_bits[k][i:i+fragment_length]) )
                                    if key in frag_dic_anti[ref_seq] : 
                                        frag_dic_anti[ref_seq][key] += comp_tuple[1] # update score
                                    else :
                                        frag_dic_anti[ref_seq][key] = comp_tuple[1] # save score
                                else :
                                    if comp_bits[k][i:i+fragment_length] in frag_dic_anti[ref_seq] : 
                                        frag_dic_anti[ref_seq][comp_bits[k][i:i+fragment_length]] += comp_tuple[1] # update score
                                    else :
                                        frag_dic_anti[ref_seq][comp_bits[k][i:i+fragment_length]] = comp_tuple[1] # save score
                            else :
                                sys.stderr.write("WARNING in fragments_from_database() strand %s frag %s doesn't have direction (%s) pdb_ids=%s\n" % (strand_seq,strand_bit[i:i+fragment_length],dir_info,str(strange_strand_dic_or_pickle[strand_seq].pdb_ids)))  

    # convert to list of tuple and sort it. Use Comp class if asked
    for frag in frag_dic :
        frag_dic[frag]= list(frag_dic[frag].items())     
        frag_dic[frag]= [ Comp(complementary_sequence=tup[0][0],hb_pattern=tup[0][1],count_score=tup[1]) for tup in frag_dic[frag] ]
    if add_promiscuity_score :
        if debug: print('adding promiscuity score to %d parallel reference fragments of length %d...' % (len(frag_dic),fragment_length))
        frag_dic=add_promiscuity_score_to_fragment_dict(frag_dic,is_antiparallel=False)
    if debug : print('now sorting parallel fragments...')
    for frag in frag_dic :
        if add_promiscuity_score :
            frag_dic[frag].sort( key = lambda x : x.promiscuity_score, reverse=False) #smaller promiscuity scores first
        frag_dic[frag].sort( key = lambda x : x.count_score, reverse=True) #higher scores first
        
    # picke_it
    if out_pickle_name is None : 
        out_pickle_name='frag_dic_'+str(fragment_length)+'.pkl'
    if frag_dic_anti!={} :
        # convert to list of tuple and sort it
        for frag in frag_dic_anti :
            frag_dic_anti[frag]= list(frag_dic_anti[frag].items())     
            frag_dic_anti[frag]= [ Comp(complementary_sequence=tup[0][0],hb_pattern=tup[0][1],count_score=tup[1]) for tup in frag_dic_anti[frag] ]
        if add_promiscuity_score :
            if debug: print('adding promiscuity score to %d antiparallel reference fragments of length %d' % (len(frag_dic_anti),fragment_length))
            frag_dic_anti=add_promiscuity_score_to_fragment_dict(frag_dic_anti,is_antiparallel=True)
        if debug:print('now sorting antiparallel fragments')
        for frag in frag_dic_anti :
            if add_promiscuity_score :
                frag_dic_anti[frag].sort( key = lambda x : x.promiscuity_score, reverse=False) #smaller promiscuity scores first
            frag_dic_anti[frag].sort( key = lambda x : x.count_score, reverse=True) #higher scores first
        
        if calculate_enrichment_tables and fragment_length==1 :
            compute_enrichment_scores(frag_dic_anti ,'-',get_score_distributions=False)
            print_enrichment_table(frag_dic_anti,'-')
            compute_enrichment_scores(frag_dic ,'+',get_score_distributions=False)
            print_enrichment_table(frag_dic,'+')
        
        output=open('para_'+out_pickle_name,'wb')
        pickle.dump(frag_dic, output)
        output.close()
        sys.stdout.write('\n%d fragments in output pickled to para_%s (max possible number=%d)\n\n' % (len(frag_dic),out_pickle_name,20**fragment_length))
        output=open('anti_'+out_pickle_name,'wb')
        pickle.dump(frag_dic_anti, output)
        output.close()
        sys.stdout.write('\n%d fragments in output pickled to anti_%s (max possible number=%d)\n\n' % (len(frag_dic_anti),out_pickle_name,20**fragment_length))
    else :
        output=open(out_pickle_name,'wb')
        pickle.dump(frag_dic, output)
        output.close()
        sys.stdout.write('\n%d fragments in output pickled to %s (max possible number=%d)\n\n' % (len(frag_dic),out_pickle_name,20**fragment_length))
    return frag_dic,frag_dic_anti




def compute_enrichment_scores(BSn, direction, hb_tolerance_level=0, save_contingency_table=False,get_score_distributions=True, debug=False ):
    '''
I added a patch for proline, which works only for fragments of size 1.
BUG: in parallel strand we are discarding pairs that are followed or preceeded by a proline as they won't have (-1,1)
NB thought for antiparalle, try to think if it need adaptation).
    BSn is a single BSn library, such as C.antiparallel_fragments[0] for fragments of lenght 1 and so on..)
    A word document named Beta pairings statistics.docx explains this procedure

enrichments,pvals=beta_strands.compute_enrichment_scores(C.antiparallel_fragments[2],'-',get_score_distributions=True)
    '''
    possible_patterns=create_possible_hb_patterns( len(list(BSn.keys())[0]) ) # possible hb pattern for this sequence size, return tuple(anti1), tuple(anti2), tuple(par1), tuple(par2)
    if direction=='-' : possible_patterns=possible_patterns[:2]
    elif direction=='+' :possible_patterns=possible_patterns[2:]
    else : 
        sys.stderr.write("**WARNING** direction %s not process (it should be either + or -)\n" % (str(direction)))
    if get_score_distributions :
        pvals=[]
        enrichments=[]
    print(direction,'possible_patterns',possible_patterns)
    # first calculate global counts necessary to get expected values and p values.
    Bi_counts={} # number of time a seq (with one HB pattern) is found in a beta strand, keys are HB patterns 
    for p in possible_patterns :
        Bi_counts[p]={} # new keys will be seq
    for seq in BSn :
        for comp in BSn[seq] :
            jp = map_fragment_hb_on_sequence_pattern( comp.hb_pattern, 0, possible_patterns, tolerance_level=hb_tolerance_level)
            if jp is None : # not matching HB pattern, not considering 
                if len(comp.complementary_sequence)==1 :
                    if direction =='-' :
                        if comp.complementary_sequence=='P' and comp.hb_pattern[0][0]==0 and comp.hb_pattern[0][1] is None :
                            #print 'Proline allowed in Bi',seq,comp 
                            jp=0 # as (0,0) prolines can only make C-O hb not N-H
                        elif seq=='P' and comp.hb_pattern[0][1]==0 and comp.hb_pattern[0][0] is None : jp=0
                        else : continue
                    elif seq=='P' and comp.hb_pattern[0][1]==1 and comp.hb_pattern[0][0] is None : jp=0 # forming HB, but P can only make one.
                    else : continue
                else : continue
            if seq not in Bi_counts[ possible_patterns[jp] ] : Bi_counts[ possible_patterns[jp] ][seq]=0
            Bi_counts[ possible_patterns[jp] ][seq]+=comp.count_score
    T_ij={} # total number of beta-pairs of this kind (n) for each HB pattern
    for p in Bi_counts :
        T_ij[p]=float(sum( Bi_counts[p].values() ))
        print(' Total',p,T_ij[p])
    # Actually calculate scores
    
    if direction=='+': off=1
    else : off=-1
    N=len(BSn)
    for j,seq in enumerate(BSn) :
        for i,comp in enumerate(BSn[seq]) :
            jp = map_fragment_hb_on_sequence_pattern( comp.hb_pattern, 0, possible_patterns, tolerance_level=hb_tolerance_level)
            # note that the below does not represent a problem in Tij because given the way it is constructed these were all in the count score of Bi (with the seq instead of the comp_seq HB pattern)
            if jp is None : # not matching HB pattern, not considering
                if len(comp.complementary_sequence)==1 : 
                    if direction =='-' :
                        if comp.complementary_sequence=='P' and comp.hb_pattern[0][0]==0 and comp.hb_pattern[0][1] is None :
                            #print 'Proline allowed in Eij 1',seq,comp,i ,comp.hb_pattern
                            jp=0 # as (0,0) prolines can only make C-O hb not N-H
                        elif seq=='P' and comp.hb_pattern[0][1]==0 and comp.hb_pattern[0][0] is None : jp=0
                        else : continue
                    elif seq=='P' and comp.hb_pattern[0][1]==1 and comp.hb_pattern[0][0] is None : jp=0 # forming HB, but P can only make one.
                    else : continue
                else : continue
# DO SOMETHING ABOUT THESE 
            try :
                flipped_pattern=flip_hb_pattern(comp.hb_pattern,direction)
                jp1 = map_fragment_hb_on_sequence_pattern( flipped_pattern, 0, possible_patterns, tolerance_level=hb_tolerance_level)
                if jp1 is None : # not matching HB pattern, not considering (unless proline)
                    if len(comp.complementary_sequence)==1 : # Hopefully the below for flipped Proline is correct..
                        if direction =='-' :
                            if comp.complementary_sequence=='P' and flipped_pattern[0][1]==0 and flipped_pattern[0][0] is None :
                                #print 'Proline allowed in Eij 2',seq,comp,i ,flipped_pattern
                                jp1=0 # as (0,0) prolines can only make C-O hb not N-H
                            elif seq=='P' and flipped_pattern[0][0]==0 and flipped_pattern[0][1] is None : jp1=0
                            else : continue
                        elif seq=='P' and flipped_pattern[0][1]==1 and flipped_pattern[0][0] is None : jp1=0 # forming HB, but P can only make one.
                        else : continue
                    else : continue
                li=i
                ljp=jp
                Oij= comp.count_score # count_score is the observed quantity
                Bi=Bi_counts[ possible_patterns[jp] ][seq]
                Bj=Bi_counts[ possible_patterns[jp1] ][ comp.complementary_sequence[::off] ]
                
                Eij = ( Bi * Bj )/ T_ij[ possible_patterns[jp] ] # expected value
                comp.enrichment_ratio = float(Oij)/Eij
    
                table=[ [ Oij , Bi - Oij ] , [ Bj - Oij , T_ij[possible_patterns[jp]] + Oij - Bi - Bj ] ]
                if save_contingency_table :
                    comp.enrichment_table=table
                if debug :
                    print(seq, comp.complementary_sequence,possible_patterns[jp],off)
                    print(table, end=' ')
                if Oij>5 : # chi-2 test of independence
                    chi2, p,dof, expected_freqs = scipy.stats.chi2_contingency( table)
                    if debug : print(chi2, p,dof,comp.enrichment_ratio,Oij,Eij)
                else : # fisher exact test 
                    _,p = scipy.stats.fisher_exact( table )
                    if debug : print(' Fisher ', p)
                comp.enrichment_p = p
                if get_score_distributions :
                    pvals+=[p]
                    enrichments+=[ comp.enrichment_ratio]
            except Exception :
                sys.stderr.write("\n**ERROR** raised by %s with %s\n" % (seq,str(comp)))
                sys.stderr.write("  Oij=%d Eij=%lf enrichment=%lf  Bi=%d  Bj=%d Tij=%d pat_k=%s\n" % (Oij,Eij,comp.enrichment_ratio, Bi,Bj,T_ij[possible_patterns[jp]],str(possible_patterns[jp])))
                sys.stderr.write("  table= %s\n" % (str(table)))
                raise
        if j%100==0 :
            sys.stdout.write('%d of %d -> %5.2lf %%\n' % (j,N,100.*j/N))
            sys.stdout.write("   last analysed is %s with %s\n" % (seq,str(BSn[seq][li])))
            sys.stdout.write("   Oij=%d Eij=%lf enrichment=%lf  Bi=%d  Bj=%d Tij=%d\n" % (Oij,Eij,BSn[seq][li].enrichment_ratio, Bi,Bj,T_ij[possible_patterns[ljp]]))
            sys.stdout.write("   table= %s  p=%s\n" % (str(table),repr(p)))
            sys.stdout.flush()
    print('\n\ndone')
    if get_score_distributions :
        return enrichments,pvals

def print_enrichment_table(BSn,direction, outfile=True,only_patter=None,hb_tolerance_level=0, save_pvalues=True, outfolder=''):
    '''
    both rows and column will be seqs in N to C term direction
      if outfile is True pvalues won't be saved and the table can be used to rank complementary peptides by looking at the rows 
      for the amino acid in the target sequence and the columns for that in the complementary peptide (this makes a difference only for parallel strands)
    '''
    #close_file=False
    #if type(outfile) is str :
    #    outfile=open(outfile,'w')
    #    close_file=True
    if outfolder is None or outfolder=='.' : outfolder=''
    elif outfolder!='' and outfolder[-1]!='/' : outfolder+='/'
    possible_patterns=create_possible_hb_patterns( len(list(BSn.keys())[0]) ,only_patterns_for_direction=direction)
    nl=sorted(BSn.keys())
    if type(outfile)==bool and outfile==True : save_pvalues=False
    if direction=='+': 
        off=1
    else : 
        off=-1
    pat_ids=[0,1]
    pairs_with_score=0
    pairs_without_score=0
    data={}
    for pi in pat_ids :
        print('Key',pi,'-->',possible_patterns[pi])
        data[pi]=csv_dict.Data()
        for j,s in enumerate(nl) :
            data[pi].hd[s]=j
            data[pi][s]=['' for i in range(len(nl))] # fill with empty table
        
    for seq in nl :
        for i,comp in enumerate(BSn[seq]) :
            if only_patter is not None :
                if not compare_patterns(comp.hb_pattern, only_patter, tolerance_level=hb_tolerance_level) : continue
            jp = map_fragment_hb_on_sequence_pattern( comp.hb_pattern, 0, possible_patterns, tolerance_level=hb_tolerance_level)
            if jp is None : 
                if len(comp.complementary_sequence)==1 : 
                    if direction =='-' :
                        if comp.complementary_sequence=='P' and comp.hb_pattern[0][0]==0 and comp.hb_pattern[0][1] is None :
                            #print 'Proline allowed in Eij 1',seq,comp,i ,comp.hb_pattern
                            jp=0 # as (0,0) prolines can only make C-O hb not N-H
                        elif seq=='P' and comp.hb_pattern[0][1]==0 and comp.hb_pattern[0][0] is None : jp=0
                        else :
                            pairs_without_score+=1 
                            continue
                    elif seq=='P' and comp.hb_pattern[0][1]==1 and comp.hb_pattern[0][0] is None : jp=0 # forming HB, but P can only make one.
                    else :
                        pairs_without_score+=1 
                        continue
                else : 
                    pairs_without_score+=1
                    continue
            elif jp not in data :
                print("Warning jp",jp,"not in data. From",comp, ' seq=',seq)
                continue
            if hasattr(comp, 'enrichment_ratio') :
                es=comp.enrichment_ratio
                if save_pvalues :
                    p=comp.enrichment_p
                    if p< 0.001 : p='{:.1e}'.format(float(p))
                    else : p=str(numpy.round(p,3))
                    data[jp][seq][data[jp].hd[  comp.complementary_sequence[::off]]]= '%.2lf p=%s' % (es,p)
                else :data[jp][seq][data[jp].hd[  comp.complementary_sequence[::off]]]=es
                if type(es) is float : pairs_with_score+=1
                else : pairs_without_score+=1
            else : pairs_without_score+=1
    key_column_hd_name=['k','k']
    if type(outfile) is bool and outfile==True :
        key_column_hd_name=[]
        if direction=='-' :outfile=['enrichment_table_antiparallel_', 'enrichment_table_antiparallel_']
        else : outfile=['enrichment_table_parallel_', 'enrichment_table_parallel_']
        for j,jp in enumerate(pat_ids) :
            if possible_patterns[jp][0]==(None,None) :
                key_column_hd_name+=[ str(possible_patterns[jp][0]) ]
                outfile[jp]+='0.tsv'
            else :
                key_column_hd_name+=[ str(possible_patterns[jp][0]) ]
                outfile[jp]+='1.tsv'
                
    for j,jp in enumerate(pat_ids) :
        if type(outfile) is str :
            pat,fname,ext=misc.get_file_path_and_extension(outfile)
            o=pat+fname+'_'+str(jp)+ext
        else : o=outfile[j]
        data[jp].Print(outfolder+o,key_column_hd_name=key_column_hd_name[jp])
        print(' table for ',possible_patterns[jp],' saved in ',o)
    #if close_file : outfile.close()
    print('printed scores for ',pairs_with_score, ' and found ',pairs_without_score,'without score')
    return data

def flip_hb_pattern(hb_pattern, direction, convert_to_tuple=True):
    '''
    works only for relatively decent HB pattern (e.g. it would fail if 1,-1 are found in pattern of '-' strand)
    '''
    if direction=='+' :
        flipped_hb=[[None,None] for i in range(len(hb_pattern))]
        for j,rp in enumerate(hb_pattern) :
            if j>0 and rp[0]==-1 : flipped_hb[j-1][1]=1
            if j+1<len(hb_pattern) and rp[1]==1 : flipped_hb[j+1][0]=-1
        if convert_to_tuple : flipped_hb=tuple([ tuple(rp) for rp in flipped_hb ])
    elif direction=='-' : 
        if convert_to_tuple : flipped_hb=tuple([ tuple(rp[::-1]) for rp in hb_pattern[::-1]])
        else : flipped_hb= [ rp[::-1] for rp in hb_pattern[::-1] ]
    else : flipped_hb=None
    return flipped_hb
                
            
def check_property_symmetry(BSn,direction, symmetric_property='enrichment_ratio', property_to_coprint=None,sensibility=0.03,print_warnings=True):
    '''
    debug function for compute_enrichment_scores()
beta_strands.check_property_symmetry(BS1,'-','enrichment_ratio', property_to_coprint=['enrichment_table','enrichment_p'])

5-08-15 Possible PROBLEM.
beta_strands.check_property_symmetry(C.antiparallel_fragments[3],'-','count_score',print_warnings=True) # gives some violation hopefully mostly due to dimers (see sentences inside function check_hb_trace_symmetry

fragments_from_database() without strange fragments in a debug folder and see if you get symmetry.. Check DONE they are all wrong so the problem is not specific to strange dict.
    '''
    not_found_count=0 
    symmetry_violations=0 #{'-':0,'+':0}
    processed=0
    
    for s in BSn :
        for comp in BSn[s] :
            if symmetric_property in dir(comp) : # should also guarantee for the quality of HB patterns
                hpt=flip_hb_pattern(comp.hb_pattern,direction)
                if direction=='-' : 
                    s2=comp.complementary_sequence[::-1]
                    st=s[::-1]
                else : 
                    s2=comp.complementary_sequence
                    st=s
                for comp2 in BSn[s2] :
                    processed+=1
                    if comp2.complementary_sequence==st and compare_patterns(comp2.hb_pattern , hpt, tolerance_level=0)  :
                        if symmetric_property not in dir(comp2) :
                            not_found_count+=1
                            if print_warnings :print('**WARNING** symmetric_property ',symmetric_property,' not in dir(comp2) key=',s2,comp2,'for ',s,comp)
                        elif not misc.CompareFloats( getattr(comp,symmetric_property), getattr(comp2,symmetric_property), sensibility=sensibility) :
                            symmetry_violations+=1
                            if print_warnings :print('* %s symmetry violation at %s facing %s with hb %s %s= %s vs %s [Deb: %s %s %s] ' % (symmetric_property,s2,st,str(hpt),symmetric_property,repr(getattr(comp2,symmetric_property)),repr(getattr(comp,symmetric_property)),comp.complementary_sequence,s,str(comp.hb_pattern) ))
                            if property_to_coprint is not None and property_to_coprint!=[] :
                                if type(property_to_coprint) is not list : property_to_coprint=[]
                                for pr in property_to_coprint :
                                    if pr in dir(comp) and print_warnings: print('    ',pr,':',getattr(comp2,pr),getattr(comp,pr), end=' ')
                                if print_warnings : print('') # endline
    print('processed',processed)
    print('without property:',not_found_count,100.*not_found_count/processed,'%')
    print('symmetry_violations',symmetry_violations,100.*symmetry_violations/processed,'%')
    return
                
def check_strand_dict_symmetry(strand_dict):
    not_found_count={'-':0,'+':0}
    symmetry_violations={'-':0,'+':0}
    processed_d={'-':0,'+':0}
    processed=0
    for k in strand_dict :
        for j,comp in enumerate(strand_dict[k].comp_seq) :
            try :
                processed+=1
                direction=comp[1]
                facing_seq=comp[0]
                count=strand_dict[k].comp_scores[j]
                pdb=strand_dict[k].pdb_ids[comp]
                hpt=flip_hb_pattern(strand_dict[k].hb_patterns[j],direction, convert_to_tuple=(type(strand_dict[k].hb_patterns[j]) is tuple ))
                processed_d[direction]+=1
                if direction=='-' : 
                    k2=facing_seq[::-1]
                    st=k[::-1]
                else : 
                    k2=facing_seq
                    st=k
                found=False
                for jj,comp2 in enumerate(strand_dict[k2].comp_seq) :
                    direction2=comp2[1]
                    facing_seq2=comp2[0]
                    if direction==direction2 and facing_seq2==st :
                        if not compare_patterns(strand_dict[k2].hb_patterns[jj], hpt, tolerance_level=0) : continue # hb pattern must be the same
                        found=True
                        if count!=strand_dict[k2].comp_scores[jj] :
                            symmetry_violations[direction]+=1
                            print(' * comp_score symmetry violation for %s facing %s %s; scores %d %d  [matched with %s %s] %s  hb: %s jj=%d j=%d' % (k,facing_seq,direction,count,strand_dict[k2].comp_scores[jj] ,k2,facing_seq2,str(pdb),str(strand_dict[k2].hb_patterns[jj]),jj,j))
                if not found :
                    not_found_count[direction]+=1
                    print('***WARNING*** not found Symmetric for %s facing %s %s from file %s j=%d' % (k,facing_seq,direction,str(pdb),j))
            except Exception :
                print('\n**ERROR** for key',k)
                print(' comp:',j,facing_seq,direction,'file',pdb)
                raise
    print('processed',processed,processed_d)
    print('not_found_count',not_found_count)
    print('symmetry_violations',symmetry_violations)
    for k in processed_d :
        print(k,'-> ',100.*not_found_count[k]/processed_d[k],'%% missing')
    for k in processed_d :
        print(k,'-> ',100.*symmetry_violations[k]/processed_d[k],'%% symmetry_violations')
    return
            
                
                
def generate_frag_dic_folder(strand_dict_or_pickle, length_range=list(range(1,9)), only_actual_n_plets=False, strange_strand_dic_or_pickle=None, calculate_enrichment_tables=True, output_folder='ALL_FRAG_DICT_hb' ):
    '''
    Save into a folder all the various pickles corresponding to frag_dictionaries of lengths contained in the list length_range
     if only_actual_n_plets is enabled only the strands that are actually long fragment_length are considered. Otherwise strands longer than that are fragmented and considered
    '''
    if type(strand_dict_or_pickle) is not dict :
        strand_dict_or_pickle=pickle.load( open( strand_dict_or_pickle, "rb" ) ) # now it's a dict
    if strange_strand_dic_or_pickle is not None:
        if type(strange_strand_dic_or_pickle) is not dict :
            strange_strand_dic_or_pickle=pickle.load( open( strange_strand_dic_or_pickle, "rb" ) ) # now it's a dict
    if output_folder[-1]!='/' : output_folder+='/'
    os.system('mkdir '+output_folder)
    length_range=[ l for l in length_range if l>0] # remove zero or smaller eventual entries
    working_dir=os.getcwd()
    try :
        os.chdir(output_folder)
        only=False
        for l in length_range : # save one length at the time
            if only_actual_n_plets and l<3 : only=False # the default for only_actual_n_plets is False anyway, but since there are no fragments smaller than 3 this is necessary for smaller lengths
            elif only_actual_n_plets : only=True
            _=fragments_from_database(strand_dict_or_pickle,l,strange_strand_dic_or_pickle=strange_strand_dic_or_pickle,only_actual_n_plets=only, calculate_enrichment_tables=calculate_enrichment_tables)
    except Exception :
        os.chdir(working_dir)
        raise
    os.chdir(working_dir)




def fragments_from_databaseOLD(strand_dict_or_pickle, fragment_length=3, only_actual_n_plets=False, look_for_hb_pattern=True, strange_strand_dic_or_pickle=None, out_pickle_name=None, try_to_retrieve_lower_opposite_ch=False,use_Comp_class=True):
    '''
    from the strand dict it generates the BSn databases, aka the fragment dictionaries of fixed length (one parallel and one antiparallel) that we use to create complementary pepetides.
     if only_actual_n_plets is enabled only the strands that are actually long fragment_length are considered. Otherwise strands longer than that are fragmented and accounted for (default)
      try_to_retrieve_lower_opposite_ch is an obscure option used to parase the strange_strand_dic_or_pickle when given. Since this contains the non-canonical beta-strands or the exceptions the
        try_to_retrieve_lower_opposite_ch==True will try to rescue some of them.
    '''
    if type(strand_dict_or_pickle) is not dict :
        strand_dict_or_pickle=pickle.load( open( strand_dict_or_pickle, "rb" ) ) # now it's a dict
    # generate all possible fragments
    frag_dic={} # this is a dict whose keys are the fragments, 
                # while being filled entries are dict as well, each has the complementary fragment (or (fragment, hb_pattern) as a key and its score (or [score, hb_energies]) as an entry
                # after entries will be list of sorted tuples. Each tuple will have the complementary fragments and its score (e.g. ('VTLV', 2) )
                # the sorting is such that higher score ones will be first
    frag_dic_anti={} # same as above but for the anti_parallel beta strands (if info about direction is provided). 
    for strand_seq in strand_dict_or_pickle : # loop on the reference sequence in the dict database
        if len(strand_seq)>=fragment_length : # check tha the reference sequence is at least as long as fragment_length
            if only_actual_n_plets and len(strand_seq)!=fragment_length :
                continue
            for j,comp_tuple in enumerate(strand_dict_or_pickle[strand_seq].complementary) : # loop on all the possible complementary sequences for the given sequence
                for i in range(0,len(strand_seq)-fragment_length+1) : # create window of size fragment_length that slides on the sequence
                    if type(comp_tuple[0]) is tuple : #this means that we know whether it is parallel or antiparalle
                        if comp_tuple[0][1]=='+' :
                            if strand_seq[i:i+fragment_length] not in frag_dic : 
                                frag_dic[strand_seq[i:i+fragment_length]]={} # initialize as empty dictionary
                            if look_for_hb_pattern :
                                key = ( comp_tuple[0][0][i:i+fragment_length] , hb_pattern_to_tuple(strand_dict_or_pickle[strand_seq].hb_patterns[j][i:i+fragment_length]) )
                                if key in frag_dic[strand_seq[i:i+fragment_length]] : 
                                    frag_dic[strand_seq[i:i+fragment_length]][key] += comp_tuple[1] # update score
                                else :
                                    frag_dic[strand_seq[i:i+fragment_length]][key] = comp_tuple[1] # save score
                            else :  
                                if comp_tuple[0][0][i:i+fragment_length] in frag_dic[strand_seq[i:i+fragment_length]] : 
                                    frag_dic[strand_seq[i:i+fragment_length]][comp_tuple[0][0][i:i+fragment_length]] += comp_tuple[1] # update score
                                else :
                                    frag_dic[strand_seq[i:i+fragment_length]][comp_tuple[0][0][i:i+fragment_length]] = comp_tuple[1] # save score 
                        elif comp_tuple[0][1]=='-' :
                            if strand_seq[i:i+fragment_length] not in frag_dic_anti : 
                                frag_dic_anti[strand_seq[i:i+fragment_length]]={} # initialize as empty dictionary 
                            if look_for_hb_pattern :
                                key = ( comp_tuple[0][0][i:i+fragment_length] , hb_pattern_to_tuple(strand_dict_or_pickle[strand_seq].hb_patterns[j][i:i+fragment_length]) )
                                if key in frag_dic_anti[strand_seq[i:i+fragment_length]] :
                                    frag_dic_anti[strand_seq[i:i+fragment_length]][key] += comp_tuple[1] # update score
                                else :
                                    frag_dic_anti[strand_seq[i:i+fragment_length]][key] = comp_tuple[1] # save score
                            else : 
                                if comp_tuple[0][0][i:i+fragment_length] in frag_dic_anti[strand_seq[i:i+fragment_length]] :
                                    frag_dic_anti[strand_seq[i:i+fragment_length]][comp_tuple[0][0][i:i+fragment_length]] += comp_tuple[1] # update score
                                else :
                                    frag_dic_anti[strand_seq[i:i+fragment_length]][comp_tuple[0][0][i:i+fragment_length]] = comp_tuple[1] # save score 
                        else :
                            sys.stderr.write("WARNING in fragments_from_database() strand %s doesn't have direction (%s) pdb_ids=%s\n" % (strand_seq,comp_tuple[0][1],str(strand_dict_or_pickle[strand_seq].pdb_ids)))  
                    else : # here we don't know if it is parallel or anti-parallel, hence we also ignore possible hb_patterns
                        if strand_seq[i:i+fragment_length] not in frag_dic : 
                            frag_dic[strand_seq[i:i+fragment_length]]={} # initialize as empty dictionary
                        if comp_tuple[0][i:i+fragment_length] in frag_dic[strand_seq[i:i+fragment_length]] : 
                            frag_dic[strand_seq[i:i+fragment_length]][comp_tuple[0][i:i+fragment_length]] += comp_tuple[1] # update score
                        else :
                            frag_dic[strand_seq[i:i+fragment_length]][comp_tuple[0][i:i+fragment_length]] = comp_tuple[1] # save score
    if strange_strand_dic_or_pickle is not None:
        del strand_seq
        if type(strange_strand_dic_or_pickle) is not dict :
            strange_strand_dic_or_pickle=pickle.load( open( strange_strand_dic_or_pickle, "rb" ) ) # now it's a dict
        
        for strand_seq in strange_strand_dic_or_pickle : # loop on the reference sequence in the dict database
            if len(strand_seq)>=fragment_length : # check tha the reference sequere is at least as long as fragment_length
                key_strand_seq=strand_seq
                for j,comp_tuple in enumerate(strange_strand_dic_or_pickle[strand_seq].complementary) : # loop on all the possible complementary sequences for the given sequence
                    if type(comp_tuple[0]) is tuple : # if we have directional information
                        comp_seq=comp_tuple[0][0]
                        dir_info=comp_tuple[0][1]
                        hb_p=''
                        if look_for_hb_pattern : 
                            hb_p = strange_strand_dic_or_pickle[key_strand_seq].hb_patterns[j]
                    else :
                        comp_seq=comp_tuple[0]
                        dir_info=None
                    if strand_seq.upper()!=strand_seq and comp_seq.upper()!=comp_seq: #check for lower case letter, if in both strand and complement they should be amino acids in strand that do not pair, keep them if they are of opposite charge
                        if try_to_retrieve_lower_opposite_ch :

                            ###print '  bef',strand_seq,comp_seq,hb_p
                            for k,res in enumerate(strand_seq) :
                                if res.islower() and res.upper() in amino_list1:
                                    if comp_seq[k].islower() and comp_seq[k].upper() in amino_list1:
                                        if mybio.aa_charge[comp_seq[k].upper()]*mybio.aa_charge[res.upper()] < 0 : # if they have opposite charge
                                            if k+1<len(strand_seq) : strand_seq=strand_seq[:k]+strand_seq[k].upper()+strand_seq[k+1:]
                                            else : strand_seq=strand_seq[:k]+strand_seq[k].upper()
                                            if k+1<len(comp_seq) : comp_seq=comp_seq[:k]+comp_seq[k].upper()+comp_seq[k+1:]
                                            else : comp_seq=comp_seq[:k]+comp_seq[k].upper()    
                                            ###print '  aft',k,strand_seq,comp_seq
                                        
                    strand_bits, comp_bits,bits = mybio.split_string_at(strand_seq, amino_list1, second_string=comp_seq, only_longer_than=fragment_length-1)
                    if look_for_hb_pattern : hb_p_bits,_= mybio.split_at_bits(bits, hb_p, only_longer_than=fragment_length-1)
                    for k,strand_bit in enumerate(strand_bits) :
                        if only_actual_n_plets and len(strand_bit)!=fragment_length :
                            continue
                        for i in range(0,len(strand_bit)-fragment_length+1) : # create window of size fragment_length that slides on the sequence
                            if dir_info=='+' or dir_info is None :
                                if strand_bit[i:i+fragment_length] not in frag_dic : 
                                    frag_dic[strand_bit[i:i+fragment_length]]={} # initialize as empty dictionary
                                if look_for_hb_pattern :
                                    key= ( comp_bits[k][i:i+fragment_length] , hb_pattern_to_tuple(hb_p_bits[k][i:i+fragment_length]) )
                                    if key in frag_dic[strand_bit[i:i+fragment_length]] : # recall that frag_dic has already been filled with all possible fragments
                                        frag_dic[strand_bit[i:i+fragment_length]][key] += comp_tuple[1] # update score
                                    else :
                                        frag_dic[strand_bit[i:i+fragment_length]][key] = comp_tuple[1] # save score
                                else :
                                    if comp_bits[k][i:i+fragment_length] in frag_dic[strand_bit[i:i+fragment_length]] : # recall that frag_dic has already been filled with all possible fragments
                                        frag_dic[strand_bit[i:i+fragment_length]][comp_bits[k][i:i+fragment_length]] += comp_tuple[1] # update score
                                    else :
                                        frag_dic[strand_bit[i:i+fragment_length]][comp_bits[k][i:i+fragment_length]] = comp_tuple[1] # save score
                            elif dir_info=='-' :
                                if strand_bit[i:i+fragment_length] not in frag_dic_anti : 
                                    frag_dic_anti[strand_bit[i:i+fragment_length]]={} # initialize as empty dictionary
                                if look_for_hb_pattern :
                                    key= ( comp_bits[k][i:i+fragment_length] , hb_pattern_to_tuple(hb_p_bits[k][i:i+fragment_length]) )
                                    if key in frag_dic_anti[strand_bit[i:i+fragment_length]] : # recall that frag_dic_anti has already been filled with all possible fragments
                                        frag_dic_anti[strand_bit[i:i+fragment_length]][key] += comp_tuple[1] # update score
                                    else :
                                        frag_dic_anti[strand_bit[i:i+fragment_length]][key] = comp_tuple[1] # save score
                                else :
                                    if comp_bits[k][i:i+fragment_length] in frag_dic_anti[strand_bit[i:i+fragment_length]] : # recall that frag_dic_anti has already been filled with all possible fragments
                                        frag_dic_anti[strand_bit[i:i+fragment_length]][comp_bits[k][i:i+fragment_length]] += comp_tuple[1] # update score
                                    else :
                                        frag_dic_anti[strand_bit[i:i+fragment_length]][comp_bits[k][i:i+fragment_length]] = comp_tuple[1] # save score
                            else :
                                sys.stderr.write("WARNING in fragments_from_database() strand %s frag %s doesn't have direction (%s) pdb_ids=%s\n" % (strand_seq,strand_bit[i:i+fragment_length],dir_info,str(strange_strand_dic_or_pickle[strand_seq].pdb_ids)))  

    # convert to list of tuple and sort it. Use Comp class if asked
    for frag in frag_dic :
        frag_dic[frag]= list(frag_dic[frag].items())
        if use_Comp_class : 
            frag_dic[frag]= [ Comp(complementary_sequence=tup[0][0],hb_pattern=tup[0][1],count_score=tup[1]) for tup in frag_dic[frag] ]
            frag_dic[frag].sort( key = lambda x : x.count_score, reverse=True) #higher scores first
        else :
            frag_dic[frag].sort(key=lambda pair : pair[1], reverse=True) #higher scores first
    # picke_it
    if out_pickle_name is None : 
        out_pickle_name='frag_dic_'+str(fragment_length)+'.pkl'
    if frag_dic_anti!={} :
        # convert to list of tuple and sort it
        for frag in frag_dic_anti :
            frag_dic_anti[frag]= list(frag_dic_anti[frag].items())
            if use_Comp_class : 
                frag_dic_anti[frag]= [ Comp(complementary_sequence=tup[0][0],hb_pattern=tup[0][1],count_score=tup[1]) for tup in frag_dic_anti[frag] ]
                frag_dic_anti[frag].sort( key = lambda x : x.count_score, reverse=True) #higher scores first
            else : 
                frag_dic_anti[frag].sort(key=lambda pair : pair[1], reverse=True) #higher scores first
        
        output=open('para_'+out_pickle_name,'wb')
        pickle.dump(frag_dic, output)
        output.close()
        sys.stdout.write('\n%d fragments in output pickled to para_%s (max possible number=%d)\n\n' % (len(frag_dic),out_pickle_name,20**fragment_length))
        output=open('anti_'+out_pickle_name,'wb')
        pickle.dump(frag_dic_anti, output)
        output.close()
        sys.stdout.write('\n%d fragments in output pickled to anti_%s (max possible number=%d)\n\n' % (len(frag_dic_anti),out_pickle_name,20**fragment_length))
    else :
        output=open(out_pickle_name,'wb')
        pickle.dump(frag_dic, output)
        output.close()
        sys.stdout.write('\n%d fragments in output pickled to %s (max possible number=%d)\n\n' % (len(frag_dic),out_pickle_name,20**fragment_length))
    return frag_dic,frag_dic_anti

def aa_ranking_from_database(strand_dict_or_pickle, out_pickle_name=None):
    '''
    used to compute beta strand propensity of single amino acids.. 
    (not really used, we use population_from_dssp_folder() instead.)
     it returns a dictionary with all amino acids and the number of time
     they are found in a beta strand (not used for complementary pepetide generation, 
     for example we used this to correct zyggregator)
    '''
    if type(strand_dict_or_pickle) is not dict :
        strand_dict_or_pickle= pickle.load( open( strand_dict_or_pickle, "rb" ) ) # now it's a dict
    aa_score={}
    for res in amino_list1 :
        aa_score[res]=0
    for strand_seq in strand_dict_or_pickle : #each key is a strand sequence
        for res in strand_seq :
            aa_score[res]+=1
    if type(out_pickle_name) is str :
        output=open(out_pickle_name,'wb')
        pickle.dump(aa_score, output)
        output.close()
    return aa_score
    
 
 



def create_possible_hb_patterns( ref_seq_length, only_patterns_for_direction=None):
    '''
    this auxiliary function generates the four possible hb pattern maps
    '''
    anti1= [(None,None) for i in range(0,ref_seq_length)]
    anti2=anti1[:]
    par1=anti1[:]
    par2=anti2[:]
    for i in range(0,ref_seq_length,2) :
        anti1[i]=(0,0)
        par1[i]=(-1,1)
    for i in range(1,ref_seq_length,2) :    
        anti2[i]=(0,0)
        par2[i]=(-1,1)
    if only_patterns_for_direction=='+' : return tuple(par1),tuple(par2)
    elif only_patterns_for_direction=='-' : return tuple(anti1),tuple(anti2) 
    return tuple(anti1),tuple(anti2),tuple(par1),tuple(par2)



def compare_patterns(pat1, pat2, tolerance_level=0, tolerated_matches={0:[] , 1:[((-1, None),(-1,1)), ((None,1),(-1,1))] , 2:[((-1, None),(-1,1)), ((None,1),(-1,1)), ((None,0),(0,0)), ((0,None),(0,0))]}):
    '''
    compare two hb_patterns of the same length with some possibilities to tune the strickness of the comparison (as from tolerance_level, if 0 True is returned only for identical patterns)
     keys of tolerated_matches should be the tolerance_level
    '''
    if tolerance_level==0 :
        if pat1==pat2 : 
            return True
        # for parallel strand we must tolerate small defect at the extrema, which essentially makes it equivalent ot the default for tolerance_level==1 besides the length check
        if len(pat1)>2 and pat1[1:-1]==pat2[1:-1] and pat1[0][1]==pat2[0][1] and pat1[-1][0]==pat2[-1][0] and (pat1[-1][0]==-1 or pat1[0][1]==1) :
            return True
        return False
    elif tolerance_level==1 : # we can tolerate sligths uncompatibilities at the extrema, but only for parallel strands.
        if pat1==pat2 :
            return True
        if pat1[1:-1]!=pat2[1:-1] : #exclude extrema 
            return False
        if ( (pat1[0],pat2[0]) in tolerated_matches[1] or (pat2[0],pat1[0]) in tolerated_matches[1] or pat1[0]==pat2[0]) and ( (pat1[-1],pat2[-1]) in tolerated_matches[1] or (pat2[-1],pat1[-1]) in tolerated_matches[1]  or pat1[-1]==pat2[-1]) :
            return True
        else : 
            return False
    elif tolerance_level==2 : # we can tolerate sligths uncompatibilities at the extrema, also for antiparallel strands
        if pat1==pat2 :
            return True
        if pat1[1:-1]!=pat2[1:-1] : #exclude extrema 
            return False
        if ( (pat1[0],pat2[0]) in tolerated_matches[tolerance_level] or (pat2[0],pat1[0]) in tolerated_matches[tolerance_level] or pat1[0]==pat2[0]) and ( (pat1[-1],pat2[-1]) in tolerated_matches[tolerance_level] or (pat2[-1],pat1[-1]) in tolerated_matches[1]  or pat1[-1]==pat2[-1]) :
            return True
        else : 
            return False
    elif tolerance_level in tolerated_matches : # generic application that requires custom tolerated_matches, does not distinguish between extrema and not
        for i in range(0,len(pat1)) :
            if not ( (pat1[i],pat2[i]) in tolerated_matches[tolerance_level] or (pat2[i],pat1[i]) in tolerated_matches[tolerance_level] ) :
                return False
        return True 
    else :
        sys.stderr.write("\n**ERROR** invalid tolerance passed to compare_patterns()\n\n")
        return None



def classify_pattern(pattern, possible_patterns ,direction=None,direction_dic={'+':[2,3],'-':[0,1]},hb_energy=None):
    '''
    similar to the function compare_patterns() but just yields a score that reflects the pattern it is more similar to, among the possible ones
     it compare with the beginning of the possible patterns!
     direction and direction_dic can be used in case there are equivalent matches to extrapolate the correct one
    ''' 
    warn_str=''
    if len(pattern)> len(possible_patterns[0]) :
        sys.stderr.write('WARNING in classify_pattern, pattern longer than templates! %d %d respectively --> using only beginning to classify\n' % (len(pattern),len(possible_patterns[0])))
        pattern=pattern[:len(possible_patterns[0])]
    classification_scores=[0 for p in possible_patterns]
    for i,template in enumerate(possible_patterns):
        if direction is not None and direction in direction_dic:
            if i not in direction_dic[direction] : classification_scores[i]-=2 # enought to avoid choosing patterns not compatible with the direction
        for j,pat in enumerate(pattern) :
            if pat[0]==template[j][0] : classification_scores[i]+= 0.5 # useful also because sometimes one is tuple and one list
            if pat[1]==template[j][1] : classification_scores[i]+= 0.5
    maxs=max(classification_scores)
    classi_inde=classification_scores.index(maxs)
    if classification_scores.count(maxs) > 1 :
        first_hb=None
        if hb_energy is not None and direction is not None:
            first_hb=classify_pattern_from_energy(hb_energy )
            if first_hb is not None :
                for ind in direction_dic[direction] :
                    if possible_patterns[ind][0][0] is None :template_first_hb=False
                    else : template_first_hb=True
                    if template_first_hb==first_hb : classi_inde=ind
        if first_hb is None : 
            warn_str+=" ERROR in classify_pattern() ambigous classification to index %d , more than one max %s input %s %s %s " % (classi_inde,str(classification_scores),str(pattern),str(direction),str(hb_energy))
    return classi_inde, classification_scores,warn_str
        
    
def classify_pattern_from_energy(hb_energy ):
    '''
    also for positive strand the energy is associated to those that make bonds:
    hb_trace[direction.index('+')]
    [[None, None], [-1, 1], [None, None]]
    >>> hb_energy[direction.index('+')]
    [[None, None], [-2.5, -2.7], [None, None]]
    '''
    tot_en=[]
    for a,b in hb_energy :
        if a is not None and b is not None : tot_en+=[a+b]
        elif a is None and b is None : tot_en+=[None]
        elif a is not None :tot_en+=[a]
        else : tot_en+=[b]
    
    predicted_hbs=[ a is not None and (tot_en[j+1] is None or a<tot_en[j+1]) for j,a in enumerate(tot_en[:-1])]
    hb1=sum(predicted_hbs[::2])
    hb2=sum(predicted_hbs[1:][::2])
    if hb1==hb2 : return None
    return hb1>hb2
    

def map_fragment_hb_on_sequence_pattern( fragment_hb_pattern, fragment_start_pos, possible_patterns, tolerance_level=0):
    '''
     map a pattern from a sscore_pos_comp on one of the possible ones
     used by the cascade method
    returns the index of the first matching pattern among those given in possible_patterns if you are mapping short fragments (never happens in cascade) give only the pattern relevant
      to that direction (parallel/antiparallel) otherwise None,None will be misassigned to the first occurrence.
    '''
    for j,p in enumerate(possible_patterns) :
        if compare_patterns( p[fragment_start_pos : fragment_start_pos+len(fragment_hb_pattern)] , fragment_hb_pattern , tolerance_level=tolerance_level) :
            return j
    return None
































"""

FUNCTIONS TO GENERATE COMPLEMENTARY PEPTIDES (look expecially at the class Complement()

"""



def get_enrichment_score_of_peptide(target_seq,complementary_peptide, hb_pattern, direction_data_enrichment_scores,return_ernichment_str=True,print_warinings=False):
    '''
    at the moment calculates the mean enrichment of facing amino acids.
    direction_data_enrichment_scores must be a dictionary for the relevant direction (either parallel or antiparallel)
     containing two keys, 0 to be used when the residue in the target sequence does not form HB and 1 in the opposite case.
     values are two data classes as 20*20 tables, where amino acids are both the keys and the hd keys.
    
    on whole symmetric matrix I get
    numpy.percentile(enrichments_ANTIPARALLEL,[5,33.333,50,66.666,95])
    [0.53826868160594743, 0.77595825818294628, 0.92647328255995365, 1.1069343273715384, 1.8428596021803128]
    numpy.percentile(enrichmentsPARALLEL,[5,33.333,50,66.666,95])
    [0.5458825647137896, 0.88626343984923839, 1.0155921472045002, 1.1876697804787284, 2.154752396450994]
    '''
    enrichment_str=''
    mean_enrichment=0.
    none=0
    for j,aa in enumerate(target_seq) :
        if hb_pattern[j]==(None,None) : pat_k=0
        else : pat_k=1
        e = direction_data_enrichment_scores[pat_k][aa][ direction_data_enrichment_scores[pat_k].hd[complementary_peptide[j] ] ]
        if type(e) is not float :
            if print_warinings : print('ERROR in get_enrichment_score_of_peptide() pair ',j,pat_k,aa,complementary_peptide[j],'gives','|',complementary_peptide[j],'|')
            enrichment_str+='n'
            none+=1
            continue
        elif return_ernichment_str :
            if e>1.1 : enrichment_str+='+'
            elif e<0.9 : enrichment_str+='-'
            else : enrichment_str+='='
        mean_enrichment+=e
    if return_ernichment_str :
        return mean_enrichment/(len(target_seq)+none),enrichment_str
    return mean_enrichment/(len(target_seq)+none)




def split_score_pos_comp_with_hb_pattern(score_pos_comp, ref_seq_length, tolerance_level=0, number_of_possible_patterns=4,logfile=sys.stderr):
    '''
    this analyzes a score_pos_comp list as returned by get_score_pos_comp (or a merged of more score_pos_comp)
    and filters its content by hb_pattern compatibility
      clusters,possible_patterns,cluster_coverage,longer_fragment,coverage_map=split_score_pos_comp_with_hb_pattern(score_pos_comp, ref_seq_length, remove_hb_pattern_from_clusters=True, tolerance_level=0, number_of_possible_patterns=4)
    note that a single score_pos_comp generally comes either from parallel strands or from antiparallel strand hence the actual number of full clusters is generally 2 (not 4!)
    '''
    possible_patterns=create_possible_hb_patterns( ref_seq_length) # initializers of different clusters
    longer_fragment=[ 0 for i in range(0,number_of_possible_patterns)] # the length of the longest fragment in each cluster
    cluster_coverage=[0. for i in range(0,number_of_possible_patterns)] # a float with the coverage (in terms of length of the original sequence) of each cluster
    coverage_map=[ [0 for j in range(0,ref_seq_length)] for i in range(0,number_of_possible_patterns) ]
    clusters=[ [] for i in range(0,number_of_possible_patterns) ]
    None_count=0
    for C in score_pos_comp :
        k=map_fragment_hb_on_sequence_pattern(C[2].hb_pattern, C[1] , possible_patterns, tolerance_level=tolerance_level )
        if k is not None :
            seg_len=len(C[2].complementary_sequence)    
            clusters[k] += [ C ]
            coverage_map[k][C[1]:C[1]+seg_len]=[1]*seg_len
            if seg_len> longer_fragment[k] : longer_fragment[k]=seg_len
        else :
            None_count+=1
    for i in range(0,number_of_possible_patterns) :
        cluster_coverage[i]= 1.*sum(coverage_map[i])/len(coverage_map[i])
    if None_count>0 and logfile is not None:
        logfile.write("   Warning in split_score_pos_comp_with_hb_pattern found %d None mapping fragments among %d.\n" % (None_count,len(score_pos_comp)))
    return clusters,possible_patterns,cluster_coverage,longer_fragment,coverage_map
        
        
        
        
def hb_pattern_to_char(hb_pattern_of_residue):
    if hb_pattern_of_residue==(None,None) or hb_pattern_of_residue==[None,None] : return ' ' # not participating in Hbond
    if hb_pattern_of_residue==(0,0) or hb_pattern_of_residue==[0,0] : return ':' # middle of anti-parallel strand
    if hb_pattern_of_residue==(-1,1) or hb_pattern_of_residue==[-1,1] : return ';' #middle of parallel strand
    if hb_pattern_of_residue==(-1,None) or hb_pattern_of_residue==[-1,None] : return ',' # begin of parallel strand
    if hb_pattern_of_residue==(None,1) or hb_pattern_of_residue==[None,1] : return '.' #end of parallel strand
    if len(hb_pattern_of_residue)!=2 or hb_pattern_of_residue[0] not in [0,1,-1,None] or hb_pattern_of_residue[1] not in [0,1,-1,None] : return '?' #invalid combo, wrong type!
    return 'z' # valid combo but we don't know what it represents

 

def get_score_pos_comp(ref_sequence, frag_dict, sort=True, add_end_pos=False) :
    '''
    it scans ref_sequence for all possible complementary fragments
    if add_end_pos the postion of the last residue will be added at the end of each tuple
    it returns score_pos_comp a list of tuple in the format 
    (score(energy), start_position of fragment relative to ref_sequence, sequence of complementary fragment, if add_normalize_score normalized_score , if add_end_pos end_position of fragment)
    Note that if hb_patterns are in frag_dict than the comp element of each tuple is (com_sequence,corresponding_hb_pattern)
    '''
    target_len=len(list(frag_dict.keys())[0])
    ref_len=len(ref_sequence)
    score_pos_comp=[] # a list of tuples, each with score of complementary fragment, starting position (in ref_sequence) of complementary fragment, complementary fragment
    if target_len > ref_len : return []
    for i in range(0, (ref_len-target_len+1) ) :
        if ref_sequence[i:i+target_len] in frag_dict :
            for j,C in enumerate(frag_dict[ ref_sequence[i:i+target_len] ]) :
                temp_tup=(C.count_score, i, C )
                if add_end_pos :
                    temp_tup+=(i+target_len,)
                score_pos_comp += [ temp_tup  ]
    if sort : 
        score_pos_comp.sort(key=lambda x : x[2].promiscuity_score ) # first sort by promiscuity score. (smaller to larger)
        score_pos_comp.sort(key=lambda x : x[2].count_score, reverse=True) # then sort by count score (larger to smaller)
    return score_pos_comp



        
def get_constituents_from_score_pos_comp(ref_seq, score_pos_comp, complementary_seq):
    '''
    assumes that score_pos_comp has been generated using only the correct hb pattern. 
    '''
    score_pos_comp_contributing=[]
    ref_len=len(ref_seq)
    comp_len=len(complementary_seq) # can be used also for complementary sequences shorter thant the reference_sequence used to generate the score_pos_comp
    for j,C in enumerate(score_pos_comp) :
        comp_st= C[2].complementary_sequence # 2 is a Comp class now
        if comp_len==ref_len :
            if -(ref_len-C[1]-len(comp_st))==0 : off=None
            else : off=-(ref_len-C[1]-len(comp_st))
            if comp_st == complementary_seq[C[1]:off] :
                score_pos_comp_contributing +=[ C ]    
        elif comp_st in complementary_seq :
            score_pos_comp_contributing +=[ C ] 
    return score_pos_comp_contributing





def read_composition_file(comp_file,full_index_key='index_of_first_seq_res:',full_len_key='len=',remove_range_from_key=False):
    '''
    # keys will look like ('YVKIIETSM', '-', 0) (complementary peptide, direction, first_res_has_HB)
    # if the file corresponds to a coverage file, where lines containing full_index_key are present
    #  a composition dict for each different index_of_first_seq_res will be returned
    '''
    index_dict={}
    key=None
    tmp=None
    ind='comp'
    for line in open(comp_file) :
        if full_index_key in line :
            if key is not None :
                if remove_range_from_key :
                    remove_first_double1,firstbar1,remove_last_double1,lastbar1 = find_remove_range_from_composition_line(tmp[0])
                    if remove_first_double1 or remove_last_double1 : # we eploit mean enrichment line that has been added knowing the intial range (otherwise we shoul read all lines)
                        mel= tmp[1].split()[-2]
                        if '+' in mel or '=' in mel or '-' in mel : # check there is a mean enrichment line
                            st=tmp[1].find(mel)
                            p1= tmp[0][ st:st+len(mel)]
                            if p1 in key[0] : key=(p1,key[1],key[2])
                            else : sys.stderr.write(' **Warning** in read_composition_file() in remove_range_from_key region p1=%s not in key[0]=%s ind=%s mel=%s\n' % (p1,key[0],str(ind),mel))
                        else : sys.stderr.write(' **Warning** in read_composition_file() in remove_range_from_key mean enrichment line not found key=%s ind=%s mel=%s\n' % (str(key),str(ind),mel))
                if key in index_dict[ind] : # should not happen
                    if not remove_range_from_key : sys.stderr.write('**Warning** in read_composition_file() Leaving old one for ind %s key %s  - this can happen if remove_range_from_key is true (now %s)\n' % (str(ind),str(key),str(remove_range_from_key)))
                    else : 
                        for l in tmp[2:] : # try to copy different line so that final composition give a sense of all possibilities. (skip target and hb)
                            if l not in index_dict[ind][key] : index_dict[ind][key]+=[l] 
                else : index_dict[ind][key]=tmp[:]
            ind= int(line.split(full_index_key)[1].split()[0])
            if full_len_key in line :
                ind=(ind, int(line.split(full_len_key)[1].split()[0]) )
            if ind not in index_dict :
                index_dict[ind]={}
        elif '<---> ' in line :
            tmp=[line]
            comp_seq,direction= line.split('<---> ')[1].split()[:2]
            if direction[:2]=='pa' : direction='+'
            elif direction[:2]=='an' : direction='-'
        elif ':' in line or ';' in line :
            t_hb=line.split('|')
            if t_hb[1]=='' : t_hb=t_hb[2][0] # double || in line
            else : t_hb=t_hb[1][0] # single | in line 
            if t_hb==' ' : key=(comp_seq, direction, 0 )
            else : key=(comp_seq, direction, 1 )
            tmp+=[line]
        elif line!='\n' and line!='' :
            tmp+=[line]
    
    if key is not None : # save the last one
        if remove_range_from_key :
            remove_first_double1,firstbar1,remove_last_double1,lastbar1 = find_remove_range_from_composition_line(tmp[0])
            if remove_first_double1 or remove_last_double1 : # we eploit mean enrichment line that has been added knowing the intial range (otherwise we shoul read all lines)
                mel= tmp[1].split()[-2]
                if '+' in mel or '=' in mel or '-' in mel : # check there is a mean enrichment line
                    st=tmp[1].find(mel)
                    p1= tmp[0][ st:st+len(mel)]
                    if p1 in key[0] : key=(p1,key[1],key[2])
                    else : sys.stderr.write(' **Warning** in read_composition_file() in remove_range_from_key region p1=%s not in key[0]=%s ind=%s mel=%s\n' % (p1,key[0],str(ind),mel))
                else : sys.stderr.write(' **Warning** in read_composition_file() in remove_range_from_key mean enrichment line not found key=%s ind=%s mel=%s\n' % (str(key),str(ind),mel))
        if key in index_dict[ind] : # should not happen
            if not remove_range_from_key : sys.stderr.write('**Warning** in read_composition_file() Leaving old one for ind %s key %s  - this can happen if remove_range_from_key is true (now %s)\n' % (str(ind),str(key),str(remove_range_from_key)))
            else : 
                for l in tmp[2:] : # try to copy different line so that final composition give a sense of all possibilities. (skip target and hb)
                    if l not in index_dict[ind][key] : index_dict[ind][key]+=[l] 
        else : index_dict[ind][key]=tmp[:]
    if len(index_dict)==1 and list(index_dict.keys())[0]=='comp' : index_dict=index_dict['comp']       
    return index_dict





def print_compositions_form_len_dic_score_pos_comp(ref_seq, len_dic_score_pos_comp, relevant_ids, template_hb_pattern=None, interested_in_range=None,extra_str=None, extra_str_pat=None, print_hb_pattern=True):
    '''
    relevant_ids must be a list of tuples like [(4, 17), (3, 194), (3, 140)] where first tuple is the length (key of len_dic_score_pos_comp) and second the index in len_dic_score_pos_comp[length]
     extra_str can be e.g. '<-->  DSVAEQKK (Z=-0.906094)' with the complementary sequence and its aggregation propensity score.
    '''
    whole_string=''
    ref_len=len(ref_seq)
    tmp= '|'.join( ref_seq[i:i+1] for i in range(0,ref_len))
    if extra_str is not None : tmp+='  '+extra_str
    if interested_in_range is not None : 
        stap,endp=interested_in_range
        stap=2*stap
        endp=2*endp-1
        tmp=tmp[:stap]+'|'+tmp[stap:endp]+'|'+tmp[endp:]
    whole_string+=tmp+'\n'
    if template_hb_pattern is not None :
        temp_pat=' '.join( hb_pattern_to_char(template_hb_pattern[i:i+1][0]) for i in range(0,ref_len) )
        if interested_in_range is not None :
            temp_pat=temp_pat[:stap]+'|'+temp_pat[stap:endp]+'|'+temp_pat[endp:]
        if extra_str_pat is not None : temp_pat+='  '+extra_str_pat
        whole_string+=temp_pat+'\n'
    # loop on complementary_sequence composition
    for length,ind in relevant_ids :
        C=len_dic_score_pos_comp[length][ind]
        comp_st=C[2].complementary_sequence
        temp_st=(C[1]*' ') +comp_st+ (' '*(ref_len-C[1]-len(comp_st)))
        temp_st2 = '|'.join( temp_st[i:i+1] for i in range(0,ref_len))+'  '
        if interested_in_range is not None :
            temp_st2=temp_st2[:stap]+'|'+temp_st2[stap:endp]+'|'+temp_st2[endp:]
        if print_hb_pattern : 
            temp_pat=''
            for pa in C[2].hb_pattern : temp_pat+=hb_pattern_to_char(pa)
            temp_pat= (C[1]*' ')+temp_pat+((ref_len-C[1]-len(comp_st))*' ')
            temp_patt=' '.join( temp_pat[i:i+1] for i in range(0,ref_len))+'  '
            if interested_in_range is not None :
                temp_patt=temp_patt[:stap]+'|'+temp_patt[stap:endp]+'|'+temp_patt[endp:]
            whole_string+= temp_patt+'\n'
        if type(C[0]) is int : tmp='count=%2d' % (C[0]) # C[0] is the score (whether normalized or not)
        else : tmp='count=%7.5lf' % (float(C[0])) # C[0] is the score (whether normalized or not)
        tmp+='  promiscuity=%3d' % ( C[2].promiscuity_score)
        temp_st2+=tmp
        whole_string+= temp_st2+'\n'
    return whole_string
        
def print_score_pos_comp(ref_seq, score_pos_comp, out_file=sys.stdout, template_hb_pattern=None, interested_in_range=None, mark_entries_in='',extra_str='', print_only_marked=False, marked_flag='***', print_hb_pattern=True):
    if type(score_pos_comp) is dict :
        len_l=list(score_pos_comp.keys())
        len_l.sort(reverse=True)
        tmp_list=[]
        for l in len_l :
            if l>0 :
                tmp_list+=score_pos_comp[l]
            else :
                sys.stderr.write("ERROR in print_score_pos_comp type is dict but key %s not recognized\n" % (str(l)))
                return -1
        score_pos_comp=tmp_list
    if print_only_marked : marked_ids=[]
    else : marked_ids=None
    
    if score_pos_comp==[] : 
        sys.stderr.write('score_pos_comp for %s is empty\n' % (ref_seq) )
        return
    ref_len=len(ref_seq)
    tmp= '|'.join( ref_seq[i:i+1] for i in range(0,ref_len))
    if mark_entries_in!='' and mark_entries_in is not None :
        tmp+='  <-->  '+mark_entries_in+extra_str
    if interested_in_range is not None : 
        stap,endp=interested_in_range
        stap=2*stap
        endp=2*endp-1
        tmp=tmp[:stap]+'|'+tmp[stap:endp]+'|'+tmp[endp:]
    out_file.write(tmp+'\n')
    if template_hb_pattern is not None :
        temp_pat=' '.join( hb_pattern_to_char(template_hb_pattern[i:i+1][0]) for i in range(0,ref_len) )
        if interested_in_range is not None :
            temp_pat=temp_pat[:stap]+'|'+temp_pat[stap:endp]+'|'+temp_pat[endp:]
        out_file.write(temp_pat+'\n')
    # loop on complementary_sequences
    for j,C in enumerate(score_pos_comp) :
         
        comp_st=C[2].complementary_sequence
    
        temp_st=(C[1]*' ') +comp_st+ (' '*(ref_len-C[1]-len(comp_st)))
        temp_st2 = '|'.join( temp_st[i:i+1] for i in range(0,ref_len))+'  '
        if print_hb_pattern : 
            temp_pat=''
            for pa in C[2].hb_pattern : temp_pat+=hb_pattern_to_char(pa)
            temp_pat= (C[1]*' ')+temp_pat+((ref_len-C[1]-len(comp_st))*' ')
            temp_patt=' '.join( temp_pat[i:i+1] for i in range(0,ref_len))+'  '
        if type(C[0]) is int : tmp='count=%2d' % (C[0]) # C[0] is the score (whether normalized or not)
        else : tmp='count=%7.5lf' % (float(C[0])) # C[0] is the score (whether normalized or not)
        
        tmp+='  promiscuity=%3d' % ( C[2].promiscuity_score)
            
        temp_st2+=tmp
        
        tmp='id=%-4d' % (j)
        temp_st2+='  length='+str(len(comp_st))+'   '+tmp
        if mark_entries_in!='' and mark_entries_in is not None :
            if len(mark_entries_in)==ref_len :
                if -(ref_len-C[1]-len(comp_st))==0 : off=None
                else : off=-(ref_len-C[1]-len(comp_st))
                if comp_st == mark_entries_in[C[1]:off] :
                    temp_st2+='\t'+marked_flag
                    if print_only_marked : 
                        if interested_in_range is not None :
                            temp_st2=temp_st2[:stap]+'|'+temp_st2[stap:endp]+'|'+temp_st2[endp:] 
                        out_file.write( temp_st2+'\n')
                        marked_ids+=[j]
                        continue
                elif print_only_marked : continue
            elif comp_st in mark_entries_in :
                temp_st2+='\t'+marked_flag+' match=?'
                if print_only_marked : 
                    if interested_in_range is not None :
                        temp_st2=temp_st2[:stap]+'|'+temp_st2[stap:endp]+'|'+temp_st2[endp:]   
                    out_file.write( temp_st2+'\n')
                    marked_ids+=[j]
                    continue
            elif print_only_marked : continue
        #if print_only_marked :  marked_ids+=[j]
        if interested_in_range is not None :
            temp_st2=temp_st2[:stap]+'|'+temp_st2[stap:endp]+'|'+temp_st2[endp:]
        out_file.write( temp_st2+'\n')
        if print_hb_pattern :
            if interested_in_range is not None :
                temp_patt=temp_patt[:stap]+'|'+temp_patt[stap:endp]+'|'+temp_patt[endp:]
            out_file.write( temp_patt+'\n')
    return marked_ids
 

def mirror_parallel_hb_pattern(pattern):
    # now I use flip_hb_pattern instead
    mirror=[[None,None] for i in range(len(pattern))]
    for j,tup in enumerate(pattern) :
        for i,x in enumerate(tup) :
            if x==1 and j<len(pattern)-1: mirror[j+1][i-1] = -1 # here i should be 1
            elif x==-1 and j>0  :mirror[j-1][i-1]=1 # here i should be 0, so 0-1=-1 is the last position.
    return hb_pattern_to_tuple(mirror)
        
def add_promiscuity_score_to_fragment_dict(frag_dict, is_antiparallel=False,sum_count_scores=True,hb_tolerance_level=1,debug=True, overwrite_existing_scores=True):
    '''
    frag_dict must be a dictionary from the BSn library whose keys are reference fragments and values list with all possible complementaries (as Comp classes)
    '''
    # default patterns:
    
    for ref_seq in frag_dict : # loop on reference fragments
        if len(ref_seq)<3 : hb_tolerance_level=0
        for j,comp in enumerate(frag_dict[ref_seq]) : # loop on all its possible complementary fragments
            if overwrite_existing_scores or comp.promiscuity_score is None :
                if is_antiparallel : 
                    tpat=flip_hb_pattern(comp.hb_pattern, '-', convert_to_tuple=type(comp.hb_pattern) is tuple) # it will be stored from the N to C terminus in the frag_dict as a reference sequence
                    com=comp.complementary_sequence[::-1]
                else :
                    tpat=flip_hb_pattern(comp.hb_pattern, '+', convert_to_tuple=type(comp.hb_pattern) is tuple) # we offset the parallel hb pattern since if the ref_seq binds the complementary like ;_; the complementary will bind the ref_seq like _;_ 
                    com=comp.complementary_sequence
                if sum_count_scores : affinity_score_to_subtract=comp.count_score
                else : affinity_score_to_subtract=1
                try :
                    tot=0
                    for other_comp in frag_dict[com] :
                        if not compare_patterns(other_comp.hb_pattern, tpat, tolerance_level=hb_tolerance_level ) : continue
                        if sum_count_scores : tot+=other_comp.count_score
                        else : tot+=1
                except Exception :
                    sys.stderr.write( '***ERROR*** in add_promiscuity_score_to_fragment_dict() for ref seq %s with comp %s   antiparallel=%s\n' % (ref_seq,str(comp),str(is_antiparallel)) )
                    sys.stderr.flush()
                    raise
                frag_dict[ref_seq][j].promiscuity_score=tot-affinity_score_to_subtract
    return frag_dict


          
def get_promiscuity_score(comp_frag,ref_frag_dict, first_has_hb,affinity_score_to_subtract, is_antiparallel=False,sum_affinities=True,debug=True):
    '''
    # We now use add_promiscuity_score_to_fragment_dict() instead
    computes the promiscuity score, ref_frag_dict is either parallel_fragments or anti_parallel_fragments (and give is_antiparallel accordingly)
     first_has_hb indicates whether the first residue is involved in HB bonds (residue of the target sequence, in other words indicates if there are : or ; on top of the first residue).
     if None is given hb are not accounted for.
    if you set sum_affinities to False than you probably want to give affinity_score_to_subtract=1 every time
    
   BUG BUG (solved in add_promiscuity_score)
when promiscuity score for parallel fragment check that the pattern of hb should be reverted.
  if your fragment binds ;_; then the complmentary will do _;_

    '''
    if first_has_hb is not None :
        antipat= ((0, 0),(None, None), (0, 0),(None, None), (0, 0),(None, None), (0, 0),(None, None), (0, 0),(None, None), (0, 0),(None, None), (0, 0))
        parapat= ((-1, 1), (None, None),(-1, 1), (None, None),(-1, 1), (None, None),(-1, 1), (None, None),(-1, 1), (None, None),(-1, 1), (None, None),(-1, 1) )
        off=int(not first_has_hb)
        if is_antiparallel : tpat=  antipat[off:off+len(comp_frag)][::-1]
        else : tpat=  parapat[off:off+len(comp_frag)]
    if is_antiparallel : com=comp_frag[::-1]
    else : com=comp_frag
    tot=0
    
    for comp in ref_frag_dict[com] :
        #print comp
        if first_has_hb is not None :
            try :
                if debug and type(comp[0][1][0]) is not tuple or len(comp[0][1][0])!=2 :
                    print('ERROR in get_promiscuity_score when looking at HBs.. expecting hb pattern but found %s in value %s for comp %s anti=%s' % (str(comp[0][1][0]),str(comp),com,str(is_antiparallel)))
                
                if not compare_patterns(comp[0][1], tpat, tolerance_level=0 ) : continue
            
            except Exception :
                sys.stderr.write( '***ERROR*** in get_promiscuity_score for seq %s   antiparallel=%s\n' % (comp_frag,str(is_antiparallel)) )
                sys.stderr.flush()
                print('comp=',comp)
                print('comp[0]=',comp[0])
                raise
        #print '*'
        if sum_affinities : tot+=comp[-1]
        else : tot+=1
    return tot-affinity_score_to_subtract



def check_quality_of_frag_dic(ref_frag_dict,check_hb=True, remove_wrong_ones=True, flag=''):
    '''
    work in progress at the moment it only checks if the length of the HB pattern is equal to the length of the fragment...
     it now also check that the count_score is positive
    '''
    for frag in ref_frag_dict :
        remove=[]
        for j,comp in enumerate(ref_frag_dict[frag]) :
            try :
                if len(frag)!=len(comp[0][1]) : # compare sequence length to hb pattern
                    sys.stderr.write("** WARNING for frag %s in ref_dict complementary %s %s has wrongly formatted HB patterns\n" % (frag,flag,str(comp)))
                    sys.stderr.flush()
                    remove+=[j]
                if comp.count_score<1 :
                    sys.stderr.write("** WARNING for frag %s in ref_dict complementary %s %s negative count_score: %s [not removing]\n" % (frag,flag,str(comp),str(comp.count_score) ))
                    sys.stderr.flush()
            except Exception :
                remove+=[j]
                sys.stderr.write("** ERROR for frag %s in ref_dict complementary %s %s generated error:\n" % (frag,flag,str(comp)))
                print_exc(file=sys.stderr)
                sys.stderr.flush()
                pass
        if remove_wrong_ones and remove!=[] :
            remove.sort(reverse=True)
            for j in remove :
                del ref_frag_dict[frag][j]
            sys.stderr.write(" REMOVING %d complementaries from %s %s\n" % (len(remove),frag,flag))
            if ref_frag_dict[frag]==[] or ref_frag_dict[frag]==() : # if it is now empty
                del ref_frag_dict[frag]
        


def init_database(folder_or_list_of_frag_pickles,logfile=sys.stdout):
    '''
    this loads the whole fragment database, pickle names are important!!
     if it is a folder it looks for keyword frag and .pkl in file. 
     pickles with keyword anti will be saved in anti_frag_dict
    '''
    #from beta_strands import Comp,Beta_strands # does not work
    frag_dic_list=[]
    anti_frag_dic_list=[]
    anti_enrichment_scores={}
    para_enrichment_scores={}
    if type(folder_or_list_of_frag_pickles) is list : # if it is a list of pickles
        for dic in folder_or_list_of_frag_pickles : 
            if type(dic) is not dict :
                if 'anti' in dic :
                    anti_frag_dic_list += [ pickle.load( open( dic, "rb" ) ) ]# now it's a dict
                    if type(list(anti_frag_dic_list[-1].values())[0][0][0]) is tuple :
                        logfile.write("loaded anti_frag_dict from %s, %d fragments of lenght %d.\n" % (dic,len(anti_frag_dic_list[-1]),len(list(anti_frag_dic_list[-1].values())[0][0][0][0])))
                    else :
                        logfile.write("loaded anti_frag_dict from %s, %d fragments of lenght %d.\n" % (dic,len(anti_frag_dic_list[-1]),len(list(anti_frag_dic_list[-1].values())[0][0][0])))
                else :
                    frag_dic_list += [ pickle.load( open( dic, "rb" ) ) ]# now it's a dict
                    if type(list(frag_dic_list[-1].values())[0][0][0]) is tuple :
                        logfile.write("loaded frag_dict from %s, %d fragments of lenght %d.\n" % (dic,len(frag_dic_list[-1]),len(list(frag_dic_list[-1].values())[0][0][0][0])))
                    else :
                        logfile.write("loaded frag_dict from %s, %d fragments of lenght %d.\n" % (dic,len(frag_dic_list[-1]),len(list(frag_dic_list[-1].values())[0][0][0])))
    elif os.path.isdir(folder_or_list_of_frag_pickles) : # if it is a directory
        if folder_or_list_of_frag_pickles[-1]!='/' : folder_or_list_of_frag_pickles+='/'
        file_list=os.listdir(folder_or_list_of_frag_pickles)
        for fil in file_list :
            if 'frag' in fil.lower() and '.pkl' in fil :
                if 'anti' in fil.lower() :
                    anti_frag_dic_list += [ pickle.load( open( folder_or_list_of_frag_pickles+fil, "rb" ) ) ]# now it's a dict
                    if type(list(anti_frag_dic_list[-1].values())[0][0][0]) is tuple :
                        logfile.write("loaded anti_frag_dict from %s %5d fragments of lenght %d.\n" % (fil,len(anti_frag_dic_list[-1]),len(list(anti_frag_dic_list[-1].values())[0][0][0][0])))
                    else :
                        logfile.write("loaded anti_frag_dict from %s %5d fragments of lenght %d.\n" % (fil,len(anti_frag_dic_list[-1]),len(list(anti_frag_dic_list[-1].values())[0][0][0])))

                else :
                    frag_dic_list += [ pickle.load( open( folder_or_list_of_frag_pickles+fil, "rb" ) ) ]# now it's a dict
                    if type(list(frag_dic_list[-1].values())[0][0][0]) is tuple :
                        logfile.write("loaded frag_dict from %s %5d fragments of lenght %d.\n" % (fil,len(frag_dic_list[-1]),len(list(frag_dic_list[-1].values())[0][0][0][0])))
                    else :
                        logfile.write("loaded frag_dict from %s %5d fragments of lenght %d.\n" % (fil,len(frag_dic_list[-1]),len(list(frag_dic_list[-1].values())[0][0][0])))
            elif 'enrichment_table' in fil :
                if 'anti' in fil :
                    if '0' in fil :
                        anti_enrichment_scores[0]=csv_dict.Data()
                        anti_enrichment_scores[0].load(folder_or_list_of_frag_pickles+fil, verbose=False)
                        logfile.write('  loaded enrichment_table from %s to antiparallel scores with target aa not-forming HB\n' % (fil))
                    elif '1' in fil :
                        anti_enrichment_scores[1]=csv_dict.Data()
                        anti_enrichment_scores[1].load(folder_or_list_of_frag_pickles+fil, verbose=False)
                        logfile.write('  loaded enrichment_table from %s to antiparallel scores with target aa forming HB\n' % (fil))
                    else :
                        sys.stderr.write("**WARNING** when trying to load enrichment_tables file %s not recognised HB int identifier\n" % (fil))
                elif 'par' in fil :
                    if '0' in fil :
                        para_enrichment_scores[0]=csv_dict.Data()
                        para_enrichment_scores[0].load(folder_or_list_of_frag_pickles+fil, verbose=False)
                        logfile.write('  loaded enrichment_table from %s to parallel scores with target aa not-forming HB\n' % (fil))
                    elif '1' in fil :
                        para_enrichment_scores[1]=csv_dict.Data()
                        para_enrichment_scores[1].load(folder_or_list_of_frag_pickles+fil, verbose=False)
                        logfile.write('  loaded enrichment_table from %s to parallel scores with target aa forming HB\n' % (fil))
                    else :
                        sys.stderr.write("**WARNING** when trying to load enrichment_tables file %s not recognised HB int identifier\n" % (fil))
    else :
        sys.stderr.write("**ERROR** in init_database folder_or_list_of_frag_pickles=%s not recognized\n" % (str(folder_or_list_of_frag_pickles)))
        return -1
    logfile.flush()
    return frag_dic_list,anti_frag_dic_list,para_enrichment_scores,anti_enrichment_scores
    











# len_list should be sorted from longer to shorter
def one_step_right(init_score_pos_comp, len_dic_score_pos_comp, len_list, limit_new_candidates=None, track_continuations=None):
    '''
    do one cascade step to the right
    '''
    init_start=init_score_pos_comp[1]
    init_end=init_score_pos_comp[-1] # for the cascade method the score_pos_comp must be created with the end position as last element
    init_seq= init_score_pos_comp[2].complementary_sequence
    candidates=[]    
    for k in len_list : # len list is sorted, we look at longer length first
        for i in range(0, len(len_dic_score_pos_comp[k])) : # since the list is sorted we consider first high score continuations
            if len_dic_score_pos_comp[k][i][-1]> init_end and len_dic_score_pos_comp[k][i][1]< init_end and len_dic_score_pos_comp[k][i][1]> init_start : # if it ends after the end and begins before it. I.e. if tehy overlap
                t_start=len_dic_score_pos_comp[k][i][1]
                t_end=len_dic_score_pos_comp[k][i][-1]
                t_seq=len_dic_score_pos_comp[k][i][2].complementary_sequence
                if init_seq[t_start-init_start:]==t_seq[0:-(t_end-init_end)] :
                    candidates+= [ len_dic_score_pos_comp[k][i] ]
                    if type(track_continuations) is list : track_continuations+=[ (k,i) ]
                    if limit_new_candidates is not None and len(candidates)>=limit_new_candidates : # the idea is that since the score_pos_comp are soreted by score the other candidates are unimportant
                        return candidates
    return candidates

# len_list should be sorted from longer to shorter
def one_step_left(init_score_pos_comp, len_dic_score_pos_comp, len_list, limit_new_candidates=None, track_continuations=None):
    '''
    do one cascade step to the left
     track_continuations can be a list, which is filled with (k,i) where k represent the length of the fragment chosen as a candidate for this step and i its index in len_dic_score_pos_comp[k]
    '''
    init_start=init_score_pos_comp[1]
    init_end=init_score_pos_comp[-1]
    init_seq= init_score_pos_comp[2].complementary_sequence
    candidates=[]    
    for k in len_list : # len list is sorted, we look at longer length first     
        for i in range(0, len(len_dic_score_pos_comp[k])) : # since the list is sorte we consider first high score continuations
            if len_dic_score_pos_comp[k][i][-1]> init_start and len_dic_score_pos_comp[k][i][1]< init_start and len_dic_score_pos_comp[k][i][-1]< init_end: # if it ends after the start and begins before it. I.e. if tehy overlap
                t_start=len_dic_score_pos_comp[k][i][1]
                t_end=len_dic_score_pos_comp[k][i][-1]
                t_seq=len_dic_score_pos_comp[k][i][2].complementary_sequence
                if init_seq[0:t_end-init_start]==t_seq[init_start-t_start:] :
                    candidates+= [ len_dic_score_pos_comp[k][i] ]
                    if type(track_continuations) is list : track_continuations+=[ (k,i) ]
                    if limit_new_candidates is not None and len(candidates)>=limit_new_candidates : # the idea is that since the score_pos_comp are soreted by score the other candidates are unimportant
                        return candidates
    return candidates  

def grow_to_right(init_score_pos_comp, len_dic_score_pos_comp, len_list, len_ref, limit_new_candidates=None,logfile=sys.stdout, initializer_tracker=None):
    '''
    grow the peptide all the way to the right
    '''
    carry_on=True
    track_composition=None
    right_designs= [ (init_score_pos_comp[2].complementary_sequence, False, init_score_pos_comp) ] #initialize with complementary initializer sequence and flag it as incomplete 
    if type(initializer_tracker) is tuple :
        composition=[ [ initializer_tracker ] ]
        track_composition=[]
    while carry_on :
        new_right_designs=[]
        if track_composition is not None : new_composition=[]
        for i in range(0,len(right_designs)) :
            if right_designs[i][1]==False : # if it has not reach a point at which we can't find new candidates.
                if track_composition is not None : track_composition=[]
                candidates= one_step_right(right_designs[i][2], len_dic_score_pos_comp, len_list, limit_new_candidates=limit_new_candidates, track_continuations=track_composition)
                if candidates==[] : #label as finished, we can't find more candidates to carry on
                    new_right_designs+= [ ( right_designs[i][0],True, right_designs[i][2] ) ]
                    if track_composition is not None : new_composition+=[ composition[i] ]
                else : # merge new candidates for next step
                    for j,C in enumerate(candidates) :
                        finish=False
                        if C[-1]==len_ref : 
                            finish=True
                        new_right_designs+= [ ( right_designs[i][0] + C[2].complementary_sequence[-(C[-1]-right_designs[i][2][-1]):], finish, C )  ]
                        if track_composition is not None : new_composition+=[ composition[i]+[track_composition[j]] ]
            else :
                new_right_designs+= [ right_designs[i] ]
                if track_composition is not None : new_composition+=[ composition[i] ]
        right_designs=new_right_designs # update move
        if track_composition is not None : composition=new_composition
       
        
        carry_on=False  # check exit condition
        actually_ended=False # finish can be given also when there are no candidates to grow the cascade
        tcount=0
        for c in right_designs :
            if c[1]==False : carry_on=True # c[1] is finish
            else : tcount+=1
            if c[2][-1]==len_ref :actually_ended =True
        if logfile is not None :    
            logfile.write( '-----> %s\t%2d  finished=%2d at least 1 has reached end:%s\n' % (right_designs[0][0],len(right_designs),tcount,str(actually_ended)) )
    if track_composition is not None :
        return right_designs, init_score_pos_comp[1],composition # return the canditates and the original starting position and the candidate composition
    return right_designs, init_score_pos_comp[1] # return the canditates and the original starting position

def grow_to_left(init_score_pos_comp, len_dic_score_pos_comp, len_list, len_ref,limit_new_candidates=None,logfile=sys.stdout, initializer_tracker=None):
    '''
    grow the peptide all the way to the left
    '''
    carry_on=True
    track_composition=None
    left_designs= [ (init_score_pos_comp[2].complementary_sequence,False,init_score_pos_comp) ] #initialize with complementary initializer sequence and flag it as incomplete
    if type(initializer_tracker) is tuple :
        composition=[ [ initializer_tracker ] ]
        track_composition=[] 
    while carry_on :
        new_left_designs=[]
        if track_composition is not None : 
            new_composition=[]
        for i in range(0,len(left_designs)) :
            if left_designs[i][1]==False :
                if track_composition is not None : track_composition=[]
                candidates= one_step_left(left_designs[i][2], len_dic_score_pos_comp, len_list,limit_new_candidates=limit_new_candidates,track_continuations=track_composition)
                if candidates==[] : # label as finished, we can't find more candidates to carry on
                    new_left_designs+= [ ( left_designs[i][0],True, left_designs[i][2] ) ] 
                    if track_composition is not None : new_composition+=[ composition[i] ]
                else : # merge new candidates for next step
                    for j,C in enumerate(candidates) :
                        finish=False
                        if C[1]==0 : finish=True
                        new_left_designs+= [ ( C[2].complementary_sequence[:-(C[-1]-left_designs[i][2][1])] + left_designs[i][0], finish, C )  ]
                        if track_composition is not None : new_composition+=[ [track_composition[j]] +composition[i] ]
            else :
                new_left_designs+=[ left_designs[i] ]# save it again because outside the for loop left_designs get overwritten by new_left_design
                if track_composition is not None : new_composition+=[ composition[i] ]
        left_designs=new_left_designs # update move
        if track_composition is not None : composition=new_composition
        
        carry_on=False  # check exit condition
        actually_ended=False # finish can be given also when there are no candidates to grow the cascade
        tcount=0
        for c in left_designs :
            if c[1]==False : carry_on=True # c[1] is finish
            else : tcount+=1
            if c[2][1]==0 : actually_ended=True
        if logfile is not None :
            logfile.write( '<----- %s\t%2d finished=%2d at least 1 has reached end:%s\n' % (left_designs[0][0],len(left_designs),tcount,str(actually_ended)) )
    if track_composition is not None :
        return left_designs, init_score_pos_comp[-1],composition # return the canditates end the original ending position and the candidate composition
    return left_designs, init_score_pos_comp[-1] # return the canditates end the original ending position


def merge_right_left(right_designs, left_designs, len_init_seq , len_ref, right_composition=None,left_composition=None):
    '''
    it returns two lists of lists. Inner in the form [ start_pos, sequence, end_pos ]
    the first contain canditates as long as the reference sequence and the second one shorter candidates.
    '''
    complete_design=[] # list of lists. Inner in the form [ start_pos, sequence, end_pos ]
    partial_design=[]
    if left_composition is not None :
        if len(right_designs)!=len(right_composition) or len(left_designs)!=len(left_composition) :
            sys.stderr.write("**WARNING** in cascade merge_right_left len(right_designs)!=len(right_composition) %d %d OR len(left_designs)!=len(left_composition) %d %d\n" % (len(right_designs),len(right_composition),len(left_designs),len(left_composition)))
        complete_composition=[]
        incomplete_composition=[]
    # double loop to get all possible combinations
    for i,r in enumerate(right_designs) :
        for j,l in enumerate(left_designs) :
            # recall that both r and l are in the form (complementary_sequence , Bool:[has not reached the end ?], score_pos_comp like entry)
            tmp=[ l[2][1] , l[0]+r[0][len_init_seq:], r[2][-1] ]
            if len(tmp[1])==len_ref :
                complete_design+= [ tmp[:] ]
                if right_composition is not None :
                    complete_composition+=[  right_composition[i][:1]+left_composition[j][:-1]+right_composition[i][1:] ]  # put initializer at the beginning
            else :
                partial_design+= [ tmp[:] ]
                if right_composition is not None :
                    incomplete_composition+=[  right_composition[i][:1]+left_composition[j][:-1]+right_composition[i][1:] ]  # put initializer at the beginning
            del tmp
    if right_composition is not None :
        return complete_design,partial_design,complete_composition,incomplete_composition
    return complete_design,partial_design
    
  
def cascade(init_score_pos_comp, len_dic_score_pos_comp, len_list, len_ref,limit_new_candidates=None,logfile=sys.stdout, initializer_tracker=None):
    '''
    starting from a given initializer it perform the cascade and returns possible candidates
    len_list should be sorted from longer to shorter
    it returns two lists of lists. Inner in the form [ start_pos, sequence, end_pos ]
    the first contain canditates of the same length of the reference sequence and the second one shorter candidates.
    '''
    if initializer_tracker  is not None :
        right_designs,_ ,right_composition= grow_to_right(init_score_pos_comp, len_dic_score_pos_comp, len_list, len_ref,limit_new_candidates=limit_new_candidates,logfile=logfile,initializer_tracker=initializer_tracker)
        left_designs,_ ,left_composition= grow_to_left(init_score_pos_comp, len_dic_score_pos_comp, len_list, len_ref,limit_new_candidates=limit_new_candidates,logfile=logfile,initializer_tracker=initializer_tracker)
        complete_design,partial_design,complete_composition,incomplete_composition = merge_right_left(right_designs, left_designs, len(init_score_pos_comp[2].complementary_sequence), len_ref,right_composition=right_composition,left_composition=left_composition)
        return complete_design, partial_design,complete_composition,incomplete_composition
    else : 
        right_designs,_ = grow_to_right(init_score_pos_comp, len_dic_score_pos_comp, len_list, len_ref,limit_new_candidates=limit_new_candidates,logfile=logfile)
        left_designs,_ = grow_to_left(init_score_pos_comp, len_dic_score_pos_comp, len_list, len_ref,limit_new_candidates=limit_new_candidates,logfile=logfile)
        complete_design,partial_design = merge_right_left(right_designs, left_designs, len(init_score_pos_comp[2].complementary_sequence), len_ref)
        return complete_design, partial_design


def determine_cascade_intializers(score_pos_comp, begin_at_score=None, consider_extra=0,count_all_initializers_as_extra=False):
    '''
    # returns a list of score_pos_comp entries corresponding to initializers (the requirement is that they have equal energy score)
    # and a float corresponding to the next lower maximum score.
    consider_extra adds additional intializers (in addition to the ones of the highest count score already considered)
    '''
    if type(score_pos_comp[0][-1]) is not int : 
        sys.stderr.write("ERROR uncompatible input in determine_cascade_intializers() score_pos_comp should be genearted with add_end_pos=True\n\n")
        return -1
    if score_pos_comp==[] : 
        sys.stderr.write("ERROR empty score_pos_comp in determine_cascade_intializers()\n\n")
        return -1
    if begin_at_score is None :
        start_count_score=score_pos_comp[0][0] # get the score of the first, as they are supposed to be ranked
    else : 
        start_count_score=begin_at_score
    initializers=[]
    initializer_indices=[]
    next_start_count_score=None # it will be the next lower maximum score
    extras=0
    #print '--------- count_all_initializers_as_extra',count_all_initializers_as_extra,consider_extra
    for i,C in enumerate(score_pos_comp) :
        if CompareFloats( C[0],start_count_score ) :
            initializers+= [ C ]
            initializer_indices+=[i]
            if count_all_initializers_as_extra :
                extras+=1
                if type(consider_extra) is int and extras>=consider_extra :
                    if next_start_count_score is None : next_start_count_score=C[0]
                    break
        elif C[0] < start_count_score :
            if type(consider_extra) is int and extras<consider_extra :    
                initializers+= [ C ]
                initializer_indices+=[i]
                extras+=1
            elif next_start_count_score is None : 
                next_start_count_score=C[0]
                break # they are supposed to be ranked there is no point in doing more cycles
            elif next_start_count_score < C[0] : next_start_count_score=C[0] # if the fragment count score is greater than the saved one (next_start_count_score)
    return initializers,initializer_indices, next_start_count_score, extras



def frag_list_to_len_dict(list_of_frag_dict_or_pickle, use_hb_patterns=True, fast_indexing=True,logfile=sys.stdout):
    '''
    convert a list of frag_dictionaries into a len_frag_dictionaries whose outer keys are length of fragments
    in this case use_hb_patterns is meaningful only if you read from pickle, to correctly guess the length of the fragments
    fast_indexing check for some standard formats and then does the thing much faster! (more error prone)
    '''
    len_list=[0]*len(list_of_frag_dict_or_pickle) # initialize list with the length
    len_dic_frag={} # dictionary whose keys are lengths and values the corresponding frag_dic
    if fast_indexing and type(list_of_frag_dict_or_pickle[0]) is dict and len(list(list_of_frag_dict_or_pickle[0].keys())[0])==1 :
        len_list=list(range(1,len(list_of_frag_dict_or_pickle)+1))
        for l in len_list :
            len_dic_frag[l]=list_of_frag_dict_or_pickle[l-1]
    else :
        if fast_indexing and logfile is not None :
            logfile.write('in frag_list_to_len_dict() doing slow indexing as format not recognized!\n')
        for i,dic in enumerate(list_of_frag_dict_or_pickle) : 
            if type(dic) is not dict :
                list_of_frag_dict_or_pickle[i]=pickle.load( open( dic, "rb" ) ) # now it's a dict
            if use_hb_patterns : len_list[i]=len( list(list_of_frag_dict_or_pickle[i].values())[0][0][0][0] ) # len of the fragments saved in dictionary
            else : len_list[i]=len( list(list_of_frag_dict_or_pickle[i].values())[0][0][0] ) # len of the fragments saved in dictionary
            len_dic_frag[ len_list[i] ] = list_of_frag_dict_or_pickle[i]
    len_list.sort(reverse=True) # sort longer to smaller
    if logfile is not None : logfile.write("   found "+str(len(len_list))+" frag dict with length: "+str(len_list)+"\n")
    return  list_of_frag_dict_or_pickle,len_dic_frag,len_list

# complete_candidates,incomplete_candidates,len_score_pos_comp= beta_strands.complementaries_from_cascade(frag_dic_list, "EQVTNVGGAVVTGVTA", fragment_length_lower_threshold=3,limit_new_candidates=10, already_reverse_expanded=True, already_normalized=True)





def get_score_pos_comp_with_clustered_hb(len_dic_frag,len_list,target_peptide,fragment_length_lower_threshold=3,use_hb_patterns=True, hb_tolerance_level=1,logfile=sys.stdout):
    '''
    used by complementaries_from_cascade this function separates the entries in len_dic_frag according to the hb pattern they align to
    it returns len_score_pos_comp,hb_clusters_id,possible_patterns where len_score_pos_comp is a dictionary 
      that has as outher keys the elements in hb_clusters_id and as values the values (and keys) of len_dic_frag that align to that hb pattern
    '''
    len_score_pos_comp=[{} for i in range(0,4)] # 4 is the number of possible clusters according to hb_patterns
    len_ref=len(target_peptide)
    len_list=[ l for l in len_list if l>=fragment_length_lower_threshold ] # this way shorter fragments won't be considered
    if use_hb_patterns :
        hb_clusters_id=[]
        tmp_score_pos_comp = get_score_pos_comp(target_peptide, len_dic_frag[len_list[-1]], sort=True, add_end_pos=True) # use the shortest length to get the HB ids, shorter fragments are more and thus yield better coverage
        clusters,possible_patterns,cluster_coverage,_,_=split_score_pos_comp_with_hb_pattern(tmp_score_pos_comp, len_ref, tolerance_level=hb_tolerance_level, number_of_possible_patterns=4,logfile=logfile)
        # clusters is a list of list, outer list are one per possible hb pattern, each inner list contains entries from score_pos_comp that map to that hb_pattern
        cluster_size=[len(c) for c in clusters]
        
        # determine which cluster is worth considering
        for j,cov in enumerate(cluster_coverage) : 
            if cov>0.1 : # at least 10% of coverage 
                hb_clusters_id+=[j]
                
        if logfile is not None :
            logfile.write( ' HB clustering of lower length fragments '+str(len_list[-1])+'.. found ids '+str(hb_clusters_id)+' coverage '+str(cluster_coverage)+'\n')
            logfile.write( '   cluster sizes: %s\n' % ( str(cluster_size)))
        
            # save the shortes length clusters
        for hb_id in hb_clusters_id :
            len_score_pos_comp[hb_id][len_list[-1]]=clusters[hb_id][:] #actually copy it!
        del tmp_score_pos_comp,clusters
        # now save all other length, exclude the already done from the list
        for le in len_list[:-1] : 
            tmp_score_pos_comp =  get_score_pos_comp(target_peptide, len_dic_frag[le], sort=True, add_end_pos=True)
            clusters,_,_,_,_=split_score_pos_comp_with_hb_pattern(tmp_score_pos_comp, len_ref, tolerance_level=hb_tolerance_level, number_of_possible_patterns=4,logfile=logfile)
            for hb_id in hb_clusters_id :
                len_score_pos_comp[hb_id][le]=clusters[hb_id][:] #actually copy it!
            del tmp_score_pos_comp,clusters
    else :
        hb_clusters_id=[0] 
        possible_patterns=[None]       
        for le in len_list :
            len_score_pos_comp[0][le]=get_score_pos_comp(target_peptide, len_dic_frag[le], sort=True, add_end_pos=True)
    return len_score_pos_comp,hb_clusters_id,possible_patterns,len_list


def complementaries_from_cascade(len_dic_frag,len_list, target_peptide, fragment_length_lower_threshold=3,initializer_minimum_length=None,limit_new_candidates=None, use_hb_patterns=True, hb_tolerance_level=1,consider_extra_initializers=0,track_composition=False,logfile=sys.stdout):
    '''
    get the complementaries with the cascade method. (use as in the class Complement, by giving only parallel or antiparallel candidates).
    complete_candidates is a list of lists. Inner in the form [ start_pos, complementary_sequence, end_pos ]. start_pos should be 0 and end_pos len(target_peptide)
    incomplete_candidates is the same but these are shorter complementaries thus start_pos can be > 0 and end_pos < len(target_peptide)
    USE list_of_frag_dict_or_pickle,len_dic_frag,len_list= frag_list_to_len_dict(list_of_frag_dict_or_pickle, use_hb_patterns=use_hb_patterns, fast_indexing=True,logfile=logfile)  # len_list is sorted longer to smaller
    '''
    if initializer_minimum_length is None or type(initializer_minimum_length) is not int :
        initializer_minimum_length=fragment_length_lower_threshold
    len_ref=len(target_peptide)
    # cluster according to different possible hb patterns:
    len_score_pos_comp,hb_clusters_id,possible_patterns,len_list = get_score_pos_comp_with_clustered_hb(len_dic_frag,len_list,target_peptide, fragment_length_lower_threshold=fragment_length_lower_threshold,use_hb_patterns=use_hb_patterns, hb_tolerance_level=hb_tolerance_level,logfile=logfile)
    #output_list_of_tuples= [None for i in xrange(len(cluster_coverage))]
    #print output_list_of_tuples
    output_list_of_tuples=[]
    for j,hb_id in enumerate(hb_clusters_id) : # loop on allowed kinds of hydrogen bond patterns
        complete_candidates=[]
        incomplete_candidates=[]
        if track_composition :
            complete_candidate_compositions=[]
            incomplete_candidate_compositions=[]
        first_len=0
        this_len_list=len_list[:]
        #init_len_list=len_list[:]
        while complete_candidates==[] and len(this_len_list)>=1 : #and len(init_len_list)>=1 :
            remove=[]
            initializer_indices,initializers=[],[]
            cons_extra=consider_extra_initializers
            for j,le in enumerate(this_len_list) :
                if len_score_pos_comp[hb_id][le]!=[] :
                    inits,initializer_inds, _, extras = determine_cascade_intializers(len_score_pos_comp[hb_id][le], begin_at_score=None ,consider_extra=cons_extra, count_all_initializers_as_extra=(cons_extra<consider_extra_initializers)) # None<None --> False
                    
                    #print '\nDEB: ',hb_id,le,len(initializers),consider_extra_initializers,initializers
                    #out=open('DEB_score_pos_comp_hb_id_%d_le_%d.pkl' % (hb_id,le),'wb')
                    #pickle.dump(len_score_pos_comp[hb_id][le],out)
                    #out.close()
                    initializers+=inits
                    initializer_indices += [ (le,ii) for ii in initializer_inds]
                    if first_len<=0 : first_len=le # so that we do not overwrite if we enter the condtion below
                    if consider_extra_initializers>0 : # if given. This condition makes the cascade method slightly slower if no new candidates are found, as length are removed one by one at the end of this function (thus some initializer lengths may be tried multiple times).
                        if len(inits)>=cons_extra or le<=initializer_minimum_length : # if they are all found at this length or they would be too short at the next round.
                            break
                        else : 
                            cons_extra -= extras
                            continue
                    else :
                        break
                    
                else :
                    remove+=[j] 
                    del len_score_pos_comp[hb_id][le] # delete the empties
            # exit the while cycle if the initializers found are too short:
            if first_len< initializer_minimum_length : break
            
            if remove!=[] :
                for j in sorted(remove,reverse=True) : 
                    del this_len_list[j]#delete empty lengths (this_len_list is sorted, we are removing the high lengths corresponding to empty complementary fragments list)
            
            if logfile is not None :
                logfile.write("initializing cascade at length %d (hb_id=%d). %d initializer framgents selected (among %d)\n" % (first_len,hb_id, len(initializers), len(len_score_pos_comp[hb_id][le])))
                logfile.write( '  intializers: '+str(initializers)+' length used '+str(this_len_list)+'\n')
                logfile.write( '  limit_new_candidates: '+str(limit_new_candidates)+' consider_extra_initializers: '+str(consider_extra_initializers)+'\n')
            for i,initializer in enumerate(initializers) :
                init_info=(initializer[2],len(initializer[2].complementary_sequence))
                if track_composition :
                    complete_design, partial_design, complete_composition, incomplete_composition = cascade(initializer, len_score_pos_comp[hb_id], this_len_list, len_ref,limit_new_candidates=limit_new_candidates, logfile=logfile, initializer_tracker=initializer_indices[i])
                    complete_candidate_compositions += complete_composition
                    incomplete_candidate_compositions += incomplete_composition
                else :
                    complete_design, partial_design = cascade(initializer, len_score_pos_comp[hb_id], this_len_list, len_ref,limit_new_candidates=limit_new_candidates, logfile=logfile, initializer_tracker=None)
                complete_candidates+= [c+[init_info] for c in complete_design ] # complete_design[:]+[init_info]
                incomplete_candidates+= [c+[init_info] for c in partial_design ] #partial_design[:]+[init_info]
                del complete_design, partial_design
            this_len_list=this_len_list[1:] # in case no complete candidates have been found remove the longest length, so that new initializers will be shorter, note that we have alredy tried all initializers of this length for this hb pattern
            # you still need to score the candidates with energy
        #output_list_of_tuples[hb_id] = (complete_candidates,incomplete_candidates,len_score_pos_comp[hb_id],len_dic_frag,possible_patterns[hb_id])
        if track_composition : output_list_of_tuples += [(complete_candidates,incomplete_candidates,len_score_pos_comp[hb_id],len_dic_frag,possible_patterns[hb_id], complete_candidate_compositions, incomplete_candidate_compositions) ]
        else : output_list_of_tuples += [(complete_candidates,incomplete_candidates,len_score_pos_comp[hb_id],len_dic_frag,possible_patterns[hb_id]) ] 
    return output_list_of_tuples

   
    
    
#TO DO: try to init cascade also from lower length fragments




def pre_screen_intializers(whole_target_sequence, fragment_list, is_antiparallel=False, up_to_N_initializers=8,consider_extra_initializers=8, out_file=sys.stdout):
    '''
    auxiliary function of screen_initializers()
    This is particular useful when one is interested in a large epitope and wants to find out where there are 
     higher chances of designing a good complementary peptide.
    it scans a whole sequence (or part of it) to identify regions against which a peptide can more effectively be design  
    assumes that fragment list is sorted so that longer length fragments are contained at the end of the list..   
    '''
    # quickly convert to dictionary of lenght
    fragment_list,len_dic_frag,len_list= frag_list_to_len_dict(fragment_list, use_hb_patterns=True,fast_indexing=True)
    candidate_intializers={} # dic of list, outer keys are frag len
    exitl=False
    total_count=0
    for i in len_list :
        score_pos_comp=get_score_pos_comp(whole_target_sequence, len_dic_frag[i], sort=True, add_end_pos=True)
        if score_pos_comp!=[] :
            count=0
            first_score=None
            while count < len(score_pos_comp) :
                initializers,initializer_indices, first_score, _ = determine_cascade_intializers(score_pos_comp, begin_at_score=first_score, consider_extra=consider_extra_initializers ) # begin_at_score doesn't do anything as it is not implemented
                if initializers==[] : break
                count+=len(initializers) # add the number of equivalent initializers found at this length.
                total_count+=len(initializers)
                if initializers!=[] :
                    if i not in candidate_intializers : # save the length as a key
                        candidate_intializers[i]=initializers
                    else :
                        candidate_intializers[i]+=initializers
                    if total_count>up_to_N_initializers :
                        exitl=True
                        break
        if exitl:
            break
    if out_file is not None :
        out_file.write("   %d  initializers found\n" % ( total_count ))
        out_file.write( '  intializers: '+str(candidate_intializers)+'\n')
        for i in len_list :
            if i in candidate_intializers :
                out_file.write('\nlength %d N=%d\n' % (i,len(candidate_intializers[i])))
                if is_antiparallel :
                    print_score_pos_comp(whole_target_sequence, candidate_intializers[i], out_file=out_file \
                                                            , template_hb_pattern=None, interested_in_range=None, mark_entries_in='', print_only_marked=False, marked_flag='***', print_hb_pattern=True)
                else :
                    print_score_pos_comp(whole_target_sequence, candidate_intializers[i], out_file=out_file \
                                                        , template_hb_pattern=None, interested_in_range=None, mark_entries_in='', print_only_marked=False, marked_flag='***', print_hb_pattern=True)
    return candidate_intializers
   
    
         
def screen_initializers(whole_target_sequence, complement_class, up_to_N_initializers=15, out_file=sys.stdout):
    '''
    This is particular useful when one is interested in a large epitope and wants to find out where there are 
     higher chances of designing a good complementary peptide.
    it scans a whole sequence (or part of it) to identify regions against which a peptide can more effectively be design  
    assumes that fragment list is sorted so that longer length fragments are contained at the end of the list..   
parallel_hb_0, parallel_hb_1 , anti_hb_0 , anti_hb_1,outlines=beta_strands.screen_initializers(target,C)
    '''
    # First do antiparallel
    antiparallel_initializers=pre_screen_intializers(whole_target_sequence, complement_class.antiparallel_fragments, is_antiparallel=True, up_to_N_initializers=up_to_N_initializers, out_file=None)
    # the do parallel
    parallel_initializers=pre_screen_intializers(whole_target_sequence, complement_class.parallel_fragments, is_antiparallel=False, up_to_N_initializers=up_to_N_initializers, out_file=None)
    # now separate the found initializers in subclasses according to whether or not their first amino acid is involved in backbone hydrogen bonding.
    parallel_hb_0={} # these are dictionary whose keys are index on the whole_target_sequence
    parallel_hb_1={}
    anti_hb_0={}
    anti_hb_1={}
    for i in range(len(whole_target_sequence)) :
        parallel_hb_0[i]=[]
        parallel_hb_1[i]=[]
        anti_hb_0[i]=[]
        anti_hb_1[i]=[]
    for length in antiparallel_initializers :
        for comp in antiparallel_initializers[length] :
            start_ind= comp[1]
            if comp[2].hb_pattern[0]==(None,None) : anti_hb_0[start_ind]=[ [ comp[2].complementary_sequence+','+str(comp[0]), '-', length, comp[2].complementary_sequence ] ] # [ sequence of initializer, paral/antiparal, length, sequence of initializer ] # it is repeated to keep the format the same as protein_coverage 
            else :  anti_hb_1[start_ind]=[ [ comp[2].complementary_sequence+','+str(comp[0]), '-', length, comp[2].complementary_sequence ] ]
    for length in parallel_initializers :
        for comp in parallel_initializers[length] :
            start_ind= comp[1]
            if comp[2].hb_pattern[0]==(None,None) : parallel_hb_0[start_ind]=[ [ comp[2].complementary_sequence+','+str(comp[0]), '+', length, comp[2].complementary_sequence ] ] # [ sequence of initializer, paral/antiparal, length, sequence of initializer ] # it is repeated to keep the format the same as protein_coverage 
            else :  parallel_hb_1[start_ind] = [ [ comp[2].complementary_sequence+','+str(comp[0]), '-', length, comp[2].complementary_sequence ] ]
    if out_file is not None :
        closefile=False
        if type(out_file) is str :
            closefile=True
            out_file=open(out_file,'w')
        outlines = visualize_coverage( whole_target_sequence, (parallel_hb_0, parallel_hb_1 , anti_hb_0 , anti_hb_1), outfile=out_file,out_as_first_seq=True, out_as_best_seq=False,split_every=None)
        if closefile : out_file.close()
    return parallel_hb_0, parallel_hb_1 , anti_hb_0 , anti_hb_1,outlines

"""
New Aug 15 libraries (antiparallels

Length 3
C.antiparallel_fragments[2] 7908 fragments for
    processed 155570576  pairs
    without count_score: 0 0.0 %
    symmetry_violations 30 1.92838522369e-05 %
OLD: 7900
    processed 135413583
    without property: 0 0.0 %
    symmetry_violations 1980 0.00146218714263 %
Length 4
C.antiparallel_fragments[3] 93398
    processed 6737332
    without property: 0 0.0 %
    symmetry_violations 10 0.000148426706595 %
OLD 90656
    processed 5823721
    without property: 0 0.0 %
    symmetry_violations 408 0.00700583012133 %
"""



"""
#USAGE:

import beta_strands,sys,os
from beta_strands import Comp,Beta_strands
dbfolder= '/Users/francesco/Dropbox/CAM/Grafting_project/DB_DSSP/BSn_libraries'   #OLD:'/Users/francesco/Dropbox/CAM/Grafting_project/DB_DSSP/ALL_FRAG_DICT_hb'

# loads the db, since it fills the RAM and takes a while it is better to init a class and then define new ones using copy=c to obtain complentary
C=beta_strands.Complement(folder_name_with_frag_pickles=dbfolder,out_file_composition='complementary_composition.txt',out_filename='complementary_peptides.tsv') 

# this is probably a small bug (should not be there with the newer libraries) in the db creation thus after loading you should run: 
#beta_strands.check_quality_of_frag_dic(C.antiparallel_fragments[2],remove_wrong_ones=True) # Resolved in new libraries, useful to run this every time an update is done.

# BASIC USAGE
target='VVTGVTA' # part of alpha-syn nac region target='FQEAVSG'
d=beta_strands.Complement( target, copy=C, interested_in_range=None,consider_extra_initializers=4)  # if you give a target it authomatically execute the .make()

# ADVANCED USAGE with interested_in_range=[1,-1]
target='AVVTGVTAV'
interested_in_range=[1,len(target)-1]
out_tag=target[1:len(target)-1]
d = beta_strands.Complement(target, copy=C, interested_in_range=interested_in_range,initializer_minimum_length=4,consider_extra_initializers=12, out_file_composition='complementary_composition_'+out_tag+'.txt',out_filename='complementary_peptides_'+out_tag+'.tsv')

# then one can screen the results to look for peptides with favorable properties, for example:
key_list=beta_strands.screen_complementary_dic(d, lambda x,cd,hd : cd[x][hd['len_initializer']]>4 and cd[x][hd['Solubility_score']]>-0.6 and 'W' in cd[x][hd['complementary_sequence']])
# and visulaize the results, also opening in a TMP excel window if asked to
d.print_keys(key_list,outfile=sys.stdout,excel=False)



# BLAST the fragments on the human proteome
target='SFQEAVSGE' # number_of_seq_with_one_hit={0: 2, 1: 72, 2: 109} tot=109 (note that 163 is returned with range [1,8]. but 109 are unique in interesting range) on E.Coli {0: 1, 1: 30, 2: 105}
d =beta_strands.Complement(target, copy=C, interested_in_range=[1,8],consider_extra_initializers=4, out_file_composition='complementary_composition_'+target[1:8]+'.txt',out_filename='complementary_peptides_'+target[1:8]+'.tsv')
overall, _, number_of_seq_with_one_hit, tot= beta_strands.blast_coverage(d ,blastDB='/Users/francesco/Documents/proteomes/human_reference_proteome_filtered.fasta')


# Complementary Peptides coverage of a target protein or protein region
prot_seq=''
pep_len=[9,11]
outid='complementary_peptide'
start_rescount=1 # residue id of first residue in sequence
outfile=outid+'_coverage.tsv'
out_file_composition=outid+'_composition.txt'
cov_tup = beta_strands.screen_complementary_peptides(prot_seq,pep_len,start_rescount=start_rescount,initializer_minimum_length=4, condition= lambda x,cd,hd : cd[x][hd['len_initializer']]>=4 and cd[x][hd['Solubility_score']]>-0.5,save_coverage_results=outfile,copy=C,out_file_composition=out_file_composition)
_=beta_strands.visualize_coverage( prot_seq,cov_tup, outfile=sys.stdout,start_rescount=start_rescount,min_initializer_len=5)
matrix,pl=beta_strands.plot_coverage(prot_seq, cov_tup , plot=True,start_rescount=start_rescount,use_log_plus_one=False,negatives_according_to_top_bottom_of_seq=True,separate_parallel_antiparallel=True,stacked=False,just_count_start=False,min_initializer_len=None,min_complementarity_score=None,min_Solubility_score=None,figsize=(14,8),print_all_sequence=True,save=None)



# Protein Coverage (old, all coverage)
import sys,mybio
prot_seq=mybio.asyn
complementary_peptide_len=8
cov_tup= beta_strands.protein_coverage(prot_seq,complementary_peptide_len,copy=C)
_=beta_strands.visualize_coverage( prot_seq,cov_tup, outfile=sys.stdout)
parallel_hb_0, parallel_hb_1 , anti_hb_0 , anti_hb_1=cov_tup



# consider also FOR NAC2 of asyn :
target='AVVTGVTAV' # part of alpha-syn nac region
d=beta_strands.Complement( target, copy=C, interested_in_range=[1,8],consider_extra_initializers=0)  # if you give a target it authomatically execute the .make()

t2='EQVTNVGG'
d=beta_strands.Complement( t2, copy=C, interested_in_range=None)  # if you give a target it authomatically execute the .make()



[likely solved] In print_score_comp we evaluate the promiscuity score for every fragment (and even hb patterns and more) even when we print only the 'marked' fragments. Waste of CPU!
Compute organism specific specificity score: with your candidate sequence regenerate complementary candidates and see if they are found in the proteome of the organism we are interested in.
"""

class Complement() :
    '''
    the main class for generating complementary peptides.
    '''
    target_seq=''  # the target sequence we want to obtain complementaries for..
    complementary_dic={} # auxiliary dict, filled with complementaries and the printed to output csv file
    complementary_hd=[] # header of complementary_dic 
    folder_name_with_frag_database='' # folder that contains the database
    parallel_fragments=None
    antiparallel_fragments=None
    warning_file=sys.stderr
    logfile=sys.stdout
    out_file_composition=sys.stdout  # outfile where the composition of the complementary fragments generated is printed
    out_filename='complementary_peptides.tsv' # name of the output file
    fragment_length_lower_threshold=3 # smaleest size of complemntary fragments considered to generate the complementary peptides
    limit_new_candidates=50 # maximum number of new candidates per cascade step.
    score_candidates=True
    use_hb_patterns=True # if True hydrogen bonds patterns are carefully examine and only complementary peptides with the correct hydrogne bonds pattern are generated.
    hb_tolerance_level=1 # 0 does not tolerate any mismatch, 1 tolerates some mismatches at the extrema of parallel strands while 2 also of antiparalle strands. 
    debug=False
    def __init__(self, target_seq=None,out_filename='complementary_peptides.tsv' , out_file_composition='complementary_composition.txt', copy=None,interested_in_range=None, folder_name_with_frag_pickles='/Users/francesco/Dropbox/CAM/Grafting_project/DB_DSSP/BSn_libraries/',list_of_parallel_fragments=None,list_of_anti_parallel_fragments=None \
                               , fragment_length_lower_threshold=3,initializer_minimum_length=None,limit_new_candidates=50 ,pH=7. \
                              , score_candidates=True, use_hb_patterns=True, hb_tolerance_level=1,consider_extra_initializers=12,logfile=sys.stdout,consider_sub_fragments_in_scoring=False,debug=False):
        self.out_filename=out_filename
        self.consider_sub_fragments_in_scoring=consider_sub_fragments_in_scoring
        self.out_file_composition=out_file_composition
        self.interested_in_range=interested_in_range
        self.complementary_dic={}
        self.complementary_hd=[]
        self.save_memory=True
        self.len_dic_frag_antiparallel=None
        self.antiparallel_len_list=None
        self.len_dic_frag_parallel=None
        self.parallel_len_list=None
        self.get_compositon_dic=False # this is automatically switched to true if out_file_composition is not None
        self.initializer_minimum_length=initializer_minimum_length
        self.consider_extra_initializers=0 # saved after so that it overwrites copy if necessary
        self.parallel_enrichment_scores={} # two keys, 0 when row amino acid is not involved in HB, 1 when that is, each value is a 20x20 matrix with enrichment scores.
        self.antiparallel_enrichment_scores={} # two keys, 0 when row amino acid is not involved in HB, 1 when that is, each value is a 20x20 matrix with enrichment scores.
        self.debug=debug
        if copy is not None : # saves time as there's no need to reaload the whole database of complementary fragments
            # If copy is given it must be another Complementary class, which is used to initialize the new one.
            self.copy(copy, copy_also_outfilenames=False)
            self.debug=debug # overwrited debug
        else :
            self.warning_file=sys.stderr
            self.logfile=logfile
            if type(self.logfile) is str :
                self.logfile=open(self.logfile,'w')
            if folder_name_with_frag_pickles is None :
                folder_name_with_frag_pickles='ALL_FRAG_DICT_hb'
            self.folder_name_with_frag_database=folder_name_with_frag_pickles
            if list_of_parallel_fragments is None :
                list_of_parallel_fragments,list_of_anti_parallel_fragments,self.parallel_enrichment_scores,self.antiparallel_enrichment_scores=init_database(folder_name_with_frag_pickles,logfile=self.logfile)
                if list_of_parallel_fragments==[] :
                    self.warn_message('empty list list_of_parallel_fragments from database %s' % (folder_name_with_frag_pickles))
                if list_of_anti_parallel_fragments==[] :
                    self.warn_message('empty list list_of_anti_parallel_fragments from database %s' % (folder_name_with_frag_pickles))
                ### correct possible bug by checking fragments and removing 'wrong ones'. It prints warnings if there are wrong framgents, in latest libary versions it should print nothing and do nothing
                for ii in range(2,len(list_of_anti_parallel_fragments)) :
                    check_quality_of_frag_dic(list_of_anti_parallel_fragments[ii],remove_wrong_ones=True,flag='-')
                for ii in range(2,len(list_of_parallel_fragments)) :
                    check_quality_of_frag_dic(list_of_parallel_fragments[ii],remove_wrong_ones=True,flag='-')
            self.antiparallel_fragments=list_of_anti_parallel_fragments
            self.parallel_fragments=list_of_parallel_fragments
            self.fragment_length_lower_threshold=fragment_length_lower_threshold
            self.limit_new_candidates=limit_new_candidates
            self.score_candidates=score_candidates
            self.hb_tolerance_level=hb_tolerance_level
            self.use_hb_patterns=use_hb_patterns
        if consider_extra_initializers is not None and consider_extra_initializers >= self.consider_extra_initializers : # overwrite copy in case
            self.consider_extra_initializers = consider_extra_initializers # number of initializers that have a smaller score than the best(s) initializer candidate but that we wish to consider anyway (yields a better coverage).
        if initializer_minimum_length is not None and initializer_minimum_length>=self.initializer_minimum_length : # note that any int is > than None
            self.initializer_minimum_length=initializer_minimum_length
            self.log_message('consider_extra_initializers= %d' % (self.consider_extra_initializers))
        # this must be at the end
        if target_seq is not None :
            self.target_seq=target_seq
            self.make(target_seq)
        
    def warn_message(self,warn_string):
        if self.warning_file is None : return
        self.warning_file.write('**WARNING** %s\n'  % warn_string)
        self.warning_file.flush()
        return
    def log_message(self,log_string):
        if self.logfile is None : return
        self.logfile.write(log_string+'\n')
        self.logfile.flush()
        return
    def make(self,target_sequence=None):
        if target_sequence is not None : self.target_seq=target_sequence
        try :
            if self.consider_sub_fragments_in_scoring : self.complementary_dic, self.complementary_hd= self.get_complementariesSLOW(self.target_seq, interested_in_range=self.interested_in_range, debug=self.debug )
            else : self.complementary_dic, self.composition_dic= self.get_complementaries(self.target_seq, interested_in_range=self.interested_in_range,get_compositon_dic=self.get_compositon_dic, debug=self.debug )
        except Exception :
            sys.stderr.write('***ERROR*** in get_complementaries() for |%s|\n' % (self.target_seq))
            raise
        return
    def get_complementaries(self, target_peptide, interested_in_range=None, pH=7.,outfile_separator='\t',get_compositon_dic=False,print_hb_pattern_of_each=False,debug=False):
        '''
        look for all the complementaries with the cascade method and genearte a .csv file with all the scores (maybe).
        
        interested_in_range allows to input a longer sequence than the one we are interested in and to declare what we are interested
        through this keyword. if one gives interested_in_range=[2,8] (or (2,8)) complementaries to target_peptide[2:8] will be returned
        giving a longer target_peptides increase the number of results
        It also print the composition of the complementaries peptides generated in self.out_file_composition (if it is not None) 
        '''    
        solubility_score=-9999
        extra_str=''
        mean_enrichment=''
        # check whether we have a file or a filename
        close_out_file_composition=False # outfile could be stdout, hence no closing.
        if type(self.out_file_composition) is str :
            out_file_composition=open(self.out_file_composition,'w')
            close_out_file_composition=True
        else : out_file_composition=self.out_file_composition
        if out_file_composition is not None : get_compositon_dic=True
        start=0
        sta=0
        if interested_in_range is not None :
            start,end=interested_in_range
            self.log_message("\n  interested_in_range option given\n  --> actually targeting range %d,%d --> %s.\n" % (start,end,target_peptide[start:end]))
        # get complementary candidates with the cascade method
        self.log_message("\ncascade on parallel fragments...")
        if self.len_dic_frag_parallel is None : 
            _,self.len_dic_frag_parallel,self.parallel_len_list= frag_list_to_len_dict(self.parallel_fragments, use_hb_patterns=self.use_hb_patterns, fast_indexing=True,logfile=self.logfile)  # parallel_len_list is sorted longer to smaller
        
        self.parallel_list_of_tup = complementaries_from_cascade(self.len_dic_frag_parallel,self.parallel_len_list, target_peptide, initializer_minimum_length=self.initializer_minimum_length, fragment_length_lower_threshold=self.fragment_length_lower_threshold,limit_new_candidates=self.limit_new_candidates, use_hb_patterns=self.use_hb_patterns,hb_tolerance_level=self.hb_tolerance_level, consider_extra_initializers=self.consider_extra_initializers,logfile=self.logfile, track_composition=True)
        
        if self.antiparallel_fragments is not None and self.antiparallel_fragments!=[] :
            self.log_message("\ncascade on anti_parallel fragments...\n")
            if self.len_dic_frag_antiparallel is None : 
                _,self.len_dic_frag_antiparallel,self.antiparallel_len_list= frag_list_to_len_dict(self.antiparallel_fragments, use_hb_patterns=self.use_hb_patterns, fast_indexing=True,logfile=self.logfile)  # antiparallel_len_list is sorted longer to smaller
            self.antiparallel_list_of_tup =   complementaries_from_cascade(self.len_dic_frag_antiparallel,self.antiparallel_len_list, target_peptide, initializer_minimum_length=self.initializer_minimum_length, fragment_length_lower_threshold=self.fragment_length_lower_threshold,limit_new_candidates=self.limit_new_candidates, use_hb_patterns=self.use_hb_patterns,hb_tolerance_level=self.hb_tolerance_level, consider_extra_initializers=self.consider_extra_initializers, logfile=self.logfile, track_composition=True)
        
        if debug :
            print('DEB: parallel_list_of_tup',len(self.parallel_list_of_tup))
            print('DEB: antiparallel_list_of_tup',len(self.antiparallel_list_of_tup),'the two length should be the same!!')
        #if len(self.antiparallel_list_of_tup)!=len(self.parallel_list_of_tup) :
        #    self.warn_message('Something is wrong, len(self.antiparallel_list_of_tup)!=len(self.parallel_list_of_tup) %d %d, aborting...' % (len(self.antiparallel_list_of_tup),len(self.parallel_list_of_tup)))
        #    return None,None
        self.log_message("\ncascade terminated, now getting scores and writing output...\n")
        
        complementary_dic={} # store the results, keys will be (complementary_sequence, parallel/antiparallel, FirstResHB) so that it is granted that they are unique
        composition_dic={} # (only with get_compositon_dic=True) store the compositions of the complementary sequences, keys are the same of complementary_dic and values are string, containing the composition.
        hb_tag=''
        if  self.score_candidates :
            complementary_hd=['CandidateIndex', 'Complementary_sequence','parallel/antiparallel','len_initializer','initializer','FirstResHB','NumberOfPaths','Complementarity_score','Solubility_score','Solubility_bin','Count_score','Promiscuity_score','Mean_enrichment' ]
        else :
            complementary_hd=['CandidateIndex', 'Complementary_sequence','parallel/antiparallel','len_initializer','initializer','FirstResHB','NumberOfPaths','Solubility_score','Solubility_bin' ]
        if interested_in_range is not None : 
            complementary_hd[2:2]=['FullRangeCompSeq[%d:%d]' % (start,end)]
        self.complementary_hd={}
        for j,el in enumerate(complementary_hd) :
            self.complementary_hd[el]=j
        # first do parallel 
        for k,out_tuple in enumerate(self.parallel_list_of_tup) : # there should be two, one per hb_pattern
            parallel_complete, parallel_incomplete, parallel_len_score_pos_comp,_,parallel_hb_pattern,parallel_complete_composition,parallel_incomplete_composition = out_tuple
            if interested_in_range is not None :
                for jca,ca in enumerate(parallel_incomplete) :
                    if ca[0]<=start and ca[2]>=end :
                        parallel_complete+=[ca]
                        parallel_complete_composition+=[ parallel_incomplete_composition[jca] ]
            if parallel_complete==[]:
                self.warn_message("no complete parallel strand found for %s! in tup %d with hb_pattern %s\n    (consider targeting a longer or shorter bit of sequence or giving hb_tolerance_level=1)\n" % (target_peptide,k,str(parallel_hb_pattern))) 
            
            # now save the output
            multi_paths_candidates={} # this cotains the IDs of the candidates of a given complementary sequence that can be built using more than one cascade path
            for i,candidate in enumerate(parallel_complete) :
                seq=candidate[1]
                if interested_in_range is not None : 
                    off_initi_len=0
                    if candidate[-1][0].complementary_sequence not in seq[start-candidate[0]:end-candidate[0]] : # if the initializer is not fully contained in the complementary peptide after the extra aas have been removed.
                        off_initi_len=-1 # we save CPU by just subtracting 1 . This will be accurate in most cases
                    lis=[i,seq[start-candidate[0]:end-candidate[0]],seq,'+',candidate[-1][1]+off_initi_len,candidate[-1][0].complementary_sequence] # the last two are initializer length and initializer seq
                else : lis=[i,seq,'+',candidate[-1][1],candidate[-1][0].complementary_sequence] # the last two are initializer length and initializer seq
                
                # add a tag to describe HB pattern
                if self.use_hb_patterns : # note that the parallel_hb_pattern is the template pattern and is the same for all fragments in the out tuple.
                    if interested_in_range is not None :
                        if parallel_hb_pattern[start]==(None,None) or parallel_hb_pattern[start]==[None,None] : hb_tag=0 # meaning that the first residue in the complementary is not involved in HB
                        else: hb_tag=1 # meaning that the first residue in the complementary is involved in HB
                    else :
                        if parallel_hb_pattern[0]==(None,None) or parallel_hb_pattern[0]==[None,None] : hb_tag=0 # meaning that the first residue in the complementary is not involved in HB
                        else: hb_tag=1 # meaning that the first residue in the complementary is involved in HB
                lis+=[hb_tag, 1]
                if debug>3 : print('DEB: parallel candidate',candidate)
                extra_str=' <---> '+seq+' parallel'
                if HAVE_camsol :
                    if len(seq)>=7 :
                        if interested_in_range is not None : 
                            _,solubility_score = camsol_intrinsic(seq[start-candidate[0]:end-candidate[0]], pH=pH, use_old_table=RUN_CAMSOL_PUBLISHED_VERSION)
                            extra_str+=' (CamSol(range)=%lf)' % (solubility_score)
                        else : 
                            _,solubility_score = camsol_intrinsic( seq, pH=pH, use_old_table=RUN_CAMSOL_PUBLISHED_VERSION )
                            extra_str+=' (CamSol=%lf)' % (solubility_score)
                        solubility_bin=get_solubility_bin(solubility_score) # bin is 0 for <=-2 and grows by one every 0.5
                    else : 
                        solubility_score=-9999
                        solubility_bin=-9999
                key=(seq,'+',hb_tag)
                if key not in complementary_dic :
                    complementary_dic[key]= lis
                else :
                    # below function uses ind=complementary_dic[key][0] for index in parallel_complete_composition
                    multi_paths_candidates, parallel_complete_composition = update_multi_paths_candidates(key,i,complementary_dic[key][0],  parallel_complete_composition,parallel_len_score_pos_comp, multi_paths_candidates)
                    
                    #if key not in multi_paths_candidates : # handle identical peptides that have been built through different paths
                    #    multi_paths_candidates[key]=[  i , 2 ] # save the first occurrence as the reference one, 2 is the count
                    #else : multi_paths_candidates[key][1]+=1 # update the count
                    #for tup_id in parallel_complete_composition[i] :
                    #    if tup_id not in parallel_complete_composition[multi_paths_candidates[key][0]] :
                    #        parallel_complete_composition[multi_paths_candidates[key][0]]+=[tup_id] # we do this so that the final scores and the composition we print take into account all possible cascade paths that lead to this candidate.
                    continue
                extra_str_pat=None
                if  self.score_candidates :
                    score, count_score, promiscuity_score= score_candidate_from_len_dic(parallel_len_score_pos_comp, parallel_complete_composition[i])
                    
                    if self.parallel_enrichment_scores!={} and self.antiparallel_enrichment_scores!={} :
                        if interested_in_range is not None : 
                            complementary_peptide=seq[start-candidate[0]:end-candidate[0]]
                            targ=target_peptide[start:end]
                        else : 
                            complementary_peptide=seq
                            targ=target_peptide
                        mean_enrichment,enrichment_str = get_enrichment_score_of_peptide(targ, complementary_peptide, parallel_hb_pattern[start:], self.parallel_enrichment_scores,  return_ernichment_str=True)
                        extra_str_pat='       '+(' '*(start-candidate[0]))+enrichment_str+' (ME=%lf)' % (mean_enrichment)
                    
                    complementary_dic[key]+=[score,solubility_score,solubility_bin, count_score, promiscuity_score, mean_enrichment]
                    
                    
                else : complementary_dic[key]+=[ solubility_score,solubility_bin ]
                if get_compositon_dic:
                    composition_string= print_compositions_form_len_dic_score_pos_comp(target_peptide,parallel_len_score_pos_comp,parallel_complete_composition[i], extra_str=extra_str, extra_str_pat=extra_str_pat,interested_in_range=interested_in_range,template_hb_pattern=parallel_hb_pattern,print_hb_pattern=print_hb_pattern_of_each)    
                    if key in composition_dic :
                        self.warn_message("%s already present in composition_dic; OVERWRITING!" % (str(key)))
                    composition_dic[key]=composition_string # in the case of multi_paths_candidates this will be overwritten each time one is found, but at the end (below) it will be overwritten again using an updated anti_complete_composition (udated above)
            # now add back the candidates that appeared more than once
            if debug :
                if k==0 : 
                    self.parallel_multi_paths_candidates_hb0=multi_paths_candidates.copy()
                    self.parallel_complete_composition_hb0= parallel_complete_composition[:]
                    self.parallel_len_score_pos_comp_hb0=parallel_len_score_pos_comp.copy()
                else : 
                    self.parallel_multi_paths_candidates_hb1=multi_paths_candidates.copy()
                    self.parallel_complete_composition_hb1= parallel_complete_composition[:]
                    self.parallel_len_score_pos_comp_hb1=parallel_len_score_pos_comp.copy()
# BUG SOMEWHERE - scores get errouneously corrected and multiple paths are sometimes accounted for from shorter components of initialiser rather than actual new paths.
            for key in multi_paths_candidates :
                # key is already in composition_dic and in complementary_dic we overwrite it
                complementary_dic[key][self.complementary_hd['NumberOfPaths']]= multi_paths_candidates[key][1] # update the number of pathways
                i= multi_paths_candidates[key][0]
                if  self.score_candidates : 
                    score, count_score, promiscuity_score = score_candidate_from_len_dic(parallel_len_score_pos_comp, parallel_complete_composition[i])
                    if self.parallel_enrichment_scores!={} and self.antiparallel_enrichment_scores!={} : # recalculate this as we haven't saved the enrichment_str
                        seq=key[0] #key=(seq,'+',hb_tag)
                        if interested_in_range is not None : 
                            complementary_peptide=complementary_peptide=complementary_dic[key][1] # previously saved
                            sta= complementary_dic[key][2].find(complementary_peptide)
                            targ=target_peptide[start:end]
                        else : 
                            complementary_peptide=seq
                            targ=target_peptide
                        mean_enrichment,enrichment_str = get_enrichment_score_of_peptide(targ, complementary_peptide, parallel_hb_pattern[start:], self.parallel_enrichment_scores,  return_ernichment_str=True)
                        extra_str_pat='       '+(' '*sta)+enrichment_str+' (ME=%lf)' % (mean_enrichment)
                    complementary_dic[key][self.complementary_hd['Complementarity_score']]=score # update complementarity score
                    complementary_dic[key][self.complementary_hd['Count_score']]=count_score # index 8 should be aggregation propensity
                    complementary_dic[key][self.complementary_hd['Promiscuity_score']]= promiscuity_score 
                    complementary_dic[key][self.complementary_hd['Mean_enrichment']]= mean_enrichment
                    extra_str=' <---> '+key[0]+' parallel'
                    if complementary_dic[key][self.complementary_hd['Solubility_score']]!='' :  
                        if interested_in_range is not None : extra_str += ' (CamSol(range)=%lf)' % (complementary_dic[key][self.complementary_hd['Solubility_score']])
                        else : extra_str += ' (CamSol=%lf)' % (complementary_dic[key][self.complementary_hd['Solubility_score']])
                        
                if get_compositon_dic :
                    composition_dic[key]= print_compositions_form_len_dic_score_pos_comp(target_peptide,parallel_len_score_pos_comp,parallel_complete_composition[i], extra_str=extra_str, extra_str_pat=extra_str_pat,interested_in_range=interested_in_range,template_hb_pattern=parallel_hb_pattern,print_hb_pattern=print_hb_pattern_of_each)
            
        
        # now do antiparallel, the code is identical as for parallel but the quantities change.
        if self.antiparallel_fragments is not None and self.antiparallel_fragments!=[] :
            for k,out_tuple in enumerate(self.antiparallel_list_of_tup) :   
                anti_complete, anti_incomplete, anti_len_score_pos_comp,_,antiparallel_hb_pattern,anti_complete_composition,anti_incomplete_composition = out_tuple
                if interested_in_range is not None :
                    for jca,ca in enumerate(anti_incomplete) :
                        if ca[0]<=start and ca[2]>=end :
                            anti_complete+=[ca]
                            anti_complete_composition+=[ anti_incomplete_composition[jca] ]
                if anti_complete==[]:
                    self.warn_message("no complete anti parallel strand found for %s! in tup %d with hb_pattern %s\n    (consider targeting a longer or shorter bit of sequence or giving hb_tolerance_level=1)\n" % (target_peptide,k,str(antiparallel_hb_pattern)))
                
                # now save the output      
                multi_paths_candidates={} # this cotains the IDs of the candidates of a given complementary sequence that can be built using more than one cascade path
                for i,candidate in enumerate(anti_complete) :
                    seq=candidate[1]
                    if interested_in_range is not None : 
                        off_initi_len=0
                        if candidate[-1][0].complementary_sequence not in seq[start-candidate[0]:end-candidate[0]] : # if the initializer is not fully contained in the complementary peptide after the extra aas have been removed.
                            off_initi_len=-1 # we save CPU by just subtracting 1 . This will be accurate in most cases
                        lis=[i,seq[start-candidate[0]:end-candidate[0]],seq,'-',candidate[-1][1]+off_initi_len,candidate[-1][0].complementary_sequence] # the last two are initializer length and initializer seq
                    else : lis=[i,seq,'-',candidate[-1][1],candidate[-1][0].complementary_sequence] # the last two are initializer length and initializer seq
                    # add a tag to describe HB pattern
                    if self.use_hb_patterns :
                        if interested_in_range is not None : # note that the antiparallel_hb_pattern is the template pattern and is the same for all fragments in the out tuple. It is as long as the target sequence disregarding interested_in_range
                            if antiparallel_hb_pattern[start]==(None,None) or antiparallel_hb_pattern[start]==[None,None] : hb_tag=0 # meaning that the first residue in the complementary is not involved in HB
                            else: hb_tag=1 # meaning that the first residue in the complementary is involved in HB
                        else :
                            if antiparallel_hb_pattern[0]==(None,None) or antiparallel_hb_pattern[0]==[None,None] : hb_tag=0 # meaning that the first residue in the complementary is not involved in HB
                            else: hb_tag=1 # meaning that the first residue in the complementary is involved in HB
                    lis+=[hb_tag, 1]
                    if debug >3 : print('DEB: antiparallel candidate',candidate)
                    extra_str=' <---> '+seq+' antiparallel'
                    if HAVE_camsol :
                        if len(seq)>=7 :
                            if interested_in_range is not None : 
                                _,solubility_score = camsol_intrinsic(seq[start-candidate[0]:end-candidate[0]], pH=pH, use_old_table=RUN_CAMSOL_PUBLISHED_VERSION)
                                extra_str+=' (CamSol(range)=%lf)' % (solubility_score)
                            else : 
                                _,solubility_score = camsol_intrinsic(seq, pH=pH, use_old_table=RUN_CAMSOL_PUBLISHED_VERSION)
                                extra_str+=' (CamSol=%lf)' % (solubility_score)
                            #if seq=='EFLLLSVITCSV' or seq=='EFLLLSVITCSV'[::-1] : print '****** %s bit %s solubility_score %lf ph=%lf' % (seq,seq[start-candidate[0]:end-candidate[0]],solubility_score,pH)
                            solubility_bin=get_solubility_bin(solubility_score) # bin is 0 for <=-2 and grows by one every 0.5
                        else : 
                            solubility_score=-9999
                            solubility_bin=-9999
                    key=(seq,'-',hb_tag)
                    if key not in complementary_dic:
                        complementary_dic[key] = lis
                    else :
                        multi_paths_candidates, anti_complete_composition = update_multi_paths_candidates(key,i,complementary_dic[key][0],  anti_complete_composition,anti_len_score_pos_comp, multi_paths_candidates)
                        #if key not in multi_paths_candidates :
                        #    multi_paths_candidates[key]=[  complementary_dic[key][0] , 2 ] # save the first occurrence as the reference one, 2 is the count
                        #else : multi_paths_candidates[key][1]+=1 # update the count
                        #for tup_id in anti_complete_composition[i] :
                        #    if tup_id not in anti_complete_composition[multi_paths_candidates[key][0]] :
                        #        anti_complete_composition[multi_paths_candidates[key][0]]+=[tup_id] # we do this so that the final scores and the composition we print take into account all possible cascade paths that lead to this candidate.
                        continue
                    extra_str_pat=None
                    if  self.score_candidates :
                        score, count_score, promiscuity_score= score_candidate_from_len_dic(anti_len_score_pos_comp, anti_complete_composition[i])
                        if self.parallel_enrichment_scores!={} and self.antiparallel_enrichment_scores!={} :
                            if interested_in_range is not None : 
                                complementary_peptide=seq[start-candidate[0]:end-candidate[0]]
                                targ=target_peptide[start:end]
                            else : 
                                complementary_peptide=seq
                                targ=target_peptide
                            mean_enrichment,enrichment_str = get_enrichment_score_of_peptide(targ, complementary_peptide, antiparallel_hb_pattern[start:], self.antiparallel_enrichment_scores,  return_ernichment_str=True)
                            extra_str_pat='       '+(' '*(start-candidate[0]))+enrichment_str+' (ME=%lf)' % (mean_enrichment)
                        complementary_dic[key]+=[score,solubility_score,solubility_bin, count_score, promiscuity_score, mean_enrichment]
                    else : complementary_dic[key]+=[ solubility_score,solubility_bin ]
                    if get_compositon_dic :
                        composition_string= print_compositions_form_len_dic_score_pos_comp(target_peptide,anti_len_score_pos_comp,anti_complete_composition[i], extra_str=extra_str, extra_str_pat=extra_str_pat,interested_in_range=interested_in_range,template_hb_pattern=antiparallel_hb_pattern,print_hb_pattern=print_hb_pattern_of_each)    
                        if key in composition_dic :
                            self.warn_message("%s already present in composition_dic; OVERWRITING!" % (str(key)))
                        composition_dic[key]=composition_string # in the case of multi_paths_candidates this will be overwritten each time one is found, but at the end (below) it will be overwritten again using an updated anti_complete_composition (udated above)
                if debug :
                    if k==0 : 
                        self.antiparallel_multi_paths_candidates_hb0=multi_paths_candidates.copy()
                        self.anti_complete_composition_hb0=anti_complete_composition[:]
                        self.parallel_len_score_pos_comp_hb0=anti_len_score_pos_comp.copy()
                    else : 
                        self.antiparallel_multi_paths_candidates_hb1=multi_paths_candidates.copy()
                        self.anti_complete_composition_hb1=anti_complete_composition[:]
                        self.parallel_len_score_pos_comp_hb0=anti_len_score_pos_comp.copy()
                # now add back in the composition candidates that appeared more than once
                for key in multi_paths_candidates :
                    # key is already in composition_dic and in complementary_dic we overwrite it
                    complementary_dic[key][self.complementary_hd['NumberOfPaths']]= multi_paths_candidates[key][1] # update the number of pathways
                    i= multi_paths_candidates[key][0]
                    if  self.score_candidates :
                        score, count_score, promiscuity_score= score_candidate_from_len_dic(anti_len_score_pos_comp, anti_complete_composition[i])
                        if self.parallel_enrichment_scores!={} and self.antiparallel_enrichment_scores!={} : # recalculate this as we want the enrichment_str
                            seq=key[0] #key=(seq,'+',hb_tag)
                            if interested_in_range is not None : 
                                complementary_peptide=complementary_dic[key][1] # previously saved
                                sta= complementary_dic[key][2].find(complementary_peptide)
                                targ=target_peptide[start:end]
                            else : 
                                complementary_peptide=seq
                                targ=target_peptide
                            if len(targ)!=len(complementary_peptide) :
                                print("\nERROR in get_enrichment_score_of_peptide area interested_in_range",len(targ),len(complementary_peptide),interested_in_range,start,end,'targ',targ,'complementary_peptide',complementary_peptide,'antiparallel_hb_pattern',antiparallel_hb_pattern,'seq',seq)
                            mean_enrichment,enrichment_str = get_enrichment_score_of_peptide(targ, complementary_peptide, antiparallel_hb_pattern[start:], self.antiparallel_enrichment_scores,  return_ernichment_str=True)
                            extra_str_pat='       '+(' '*sta)+enrichment_str+' (ME=%lf)' % (mean_enrichment)
                        complementary_dic[key][self.complementary_hd['Complementarity_score']]=score # update complementarity score
                        complementary_dic[key][self.complementary_hd['Count_score']]=count_score # index 8 should be aggregation propensity
                        complementary_dic[key][self.complementary_hd['Promiscuity_score']]= promiscuity_score
                        complementary_dic[key][self.complementary_hd['Mean_enrichment']]= mean_enrichment
                        extra_str=' <---> '+key[0]+' antiparallel'
                        if complementary_dic[key][self.complementary_hd['Solubility_score']]!='' :  
                            if interested_in_range is not None : extra_str += ' (CamSol(range)=%lf)' % (complementary_dic[key][self.complementary_hd['Solubility_score']])
                            else : extra_str += ' (CamSol=%lf)' % (complementary_dic[key][self.complementary_hd['Solubility_score']])
                    if get_compositon_dic :
                        composition_dic[key]= print_compositions_form_len_dic_score_pos_comp(target_peptide,anti_len_score_pos_comp,anti_complete_composition[i], extra_str=extra_str,extra_str_pat=extra_str_pat,interested_in_range=interested_in_range,template_hb_pattern=antiparallel_hb_pattern,print_hb_pattern=print_hb_pattern_of_each)
        
        
        # print the output
        if self.out_filename is not None :
            self.print_table(complementary_dic, complementary_hd=complementary_hd,outfile=self.out_filename,outfile_separator=outfile_separator)
            self.log_message("\nOutput table printed in %s\n" % (str(self.out_filename)))
        if out_file_composition is not None:
            for key in complementary_dic :
                out_file_composition.write('%s\n' % (composition_dic[key]))
            
            if close_out_file_composition : out_file_composition.close()
            self.log_message("Peptide composition printed in %s\n" % (str(self.out_file_composition)))
        if self.save_memory and not debug:
            del self.antiparallel_list_of_tup,self.parallel_list_of_tup
        return complementary_dic,composition_dic
    
    def print_table(self,complementary_dic,only_these_keys=None,complementary_hd=None,outfile=sys.stdout,outfile_separator='\t'):
        if complementary_hd is None : complementary_hd=self.complementary_hd
        if type(complementary_hd) is dict :
            hd=['?']*len(complementary_hd)
            for title in complementary_hd :
                if type(complementary_hd[title]) is int and complementary_hd[title]<len(complementary_hd):
                    hd[complementary_hd[title]]=title
        else : 
            hd=complementary_hd
        closeit=False
        if type(outfile) is str :
            closeit=True
            outfile=open(outfile,'w')
        outfile.write('%s\n' % ( outfile_separator.join( map(str, hd[1:] ) ))) # remove index 0 since that is just the index in len_dic_score_pos_comp (internal useless reference)
        for key in complementary_dic :
            if only_these_keys is not None and key not in only_these_keys : continue
            outfile.write('%s\n' % ( outfile_separator.join( map(str, complementary_dic[key][1:] ) ))) # remove index 0 since that is just the index in len_dic_score_pos_comp (internal useless reference)
        if closeit : outfile.close()
        return
            
    def get_complementariesSLOW(self, target_peptide, interested_in_range=None, pH=7.,debug=False):
        '''
        look for all the complementaries with the cascade method and genearte a .csv file with all the scores (maybe).
        
        interested_in_range allows to input a longer sequence than the one we are interested in and to declare what we are interested
        through this keyword. if one gives interested_in_range=[2,8] (or (2,8)) complementaries to target_peptide[2:8] will be returned
        giving a longer target_peptides increase the number of results
        It also print the composition of the complementaries peptides generated in self.out_file_composition (if it is not None) 
        '''    
        solubility_score=''
        extra_str=''
        # check whether we have a file or a filename
        close_out_file_composition=False # outfile could be stdout, hence no closing.
        if type(self.out_file_composition) is str :
            out_file_composition=open(self.out_file_composition,'w')
            close_out_file_composition=True
        
        if interested_in_range is not None :
            start,end=interested_in_range
            self.log_message("\n  interested_in_range option given\n  --> actually targeting range %d,%d --> %s.\n" % (start,end,target_peptide[start:end]))
        # get complementary candidates with the cascade method
        self.log_message("\ncascade on parallel fragments...")
        _,self.len_dic_frag_parallel,len_list= frag_list_to_len_dict(self.parallel_fragments, use_hb_patterns=self.use_hb_patterns, fast_indexing=True,logfile=self.logfile)  # len_list is sorted longer to smaller
        self.parallel_list_of_tup = complementaries_from_cascade(self.len_dic_frag_parallel,len_list, target_peptide, initializer_minimum_length=self.initializer_minimum_length, fragment_length_lower_threshold=self.fragment_length_lower_threshold,limit_new_candidates=self.limit_new_candidates, use_hb_patterns=self.use_hb_patterns,hb_tolerance_level=self.hb_tolerance_level, consider_extra_initializers=self.consider_extra_initializers,logfile=self.logfile, track_composition=False)
        
        if self.antiparallel_fragments is not None and self.antiparallel_fragments!=[] :
            self.log_message("\ncascade on anti_parallel fragments...\n")
            _,self.len_dic_frag_antiparallel,len_list= frag_list_to_len_dict(self.antiparallel_fragments, use_hb_patterns=self.use_hb_patterns, fast_indexing=True,logfile=self.logfile)  # len_list is sorted longer to smaller
            self.antiparallel_list_of_tup =   complementaries_from_cascade(self.len_dic_frag_antiparallel,len_list, target_peptide, initializer_minimum_length=self.initializer_minimum_length, fragment_length_lower_threshold=self.fragment_length_lower_threshold,limit_new_candidates=self.limit_new_candidates, use_hb_patterns=self.use_hb_patterns,hb_tolerance_level=self.hb_tolerance_level, consider_extra_initializers=self.consider_extra_initializers,logfile=self.logfile, track_composition=False)
        
        if debug :
            print('DEB: parallel_list_of_tup',len(self.parallel_list_of_tup))
            print('DEB: antiparallel_list_of_tup',len(self.antiparallel_list_of_tup),'the two length should be the same!!')
        #if len(self.antiparallel_list_of_tup)!=len(self.parallel_list_of_tup) :
        #    self.warn_message('Something is wrong, len(self.antiparallel_list_of_tup)!=len(self.parallel_list_of_tup) %d %d, aborting...' % (len(self.antiparallel_list_of_tup),len(self.parallel_list_of_tup)))
        #    return None,None
        self.log_message("\ncascade terminated, now getting scores and writing output...\n")
        
        
        complementary_dic={}
        complementary_hd=[ 'parallel/antiparallel','len_initializer','initializer','FirstResHB','NumberOfPaths','Complmentarity_score','Solubility_score','Solubility_bin' \
                          ,'Count_score','Promiscuity_score' ]
            
        if not self.score_candidates :
            score, count_score,promiscuity_score='','',''
        for k,out_tuple in enumerate(self.parallel_list_of_tup) : # there should be two, one per hb_pattern
            parallel_complete, parallel_incomplete, parallel_len_score_pos_comp,_,parallel_hb_pattern = out_tuple
            if interested_in_range is not None :
                for ca in parallel_incomplete :
                    if ca[0]<=start and ca[2]>=end :
                        parallel_complete+=[ca]
            if parallel_complete==[]:
                self.warn_message("no complete parallel strand found for %s! in tup %d with hb_pattern %s\n    (consider targeting a longer or shorter bit of sequence or giving hb_tolerance_level=1)\n" % (target_peptide,k,str(parallel_hb_pattern))) 
            
                
            # now we need to calculate energy and specificity for the candidates we have obtained.
            # and to save the data
            
            if type(parallel_len_score_pos_comp) is dict :
                len_l=list(parallel_len_score_pos_comp.keys())
                len_l.sort(reverse=True)
                score_pos_comp=[]
                for l in len_l :
                    if l>0 :
                        score_pos_comp+=parallel_len_score_pos_comp[l]
                    else :
                        sys.stderr.write("ERROR unrecognized key in parallel_len_score_pos_comp (type is dict but key %s not recognized)\n" % (str(l)))
                        return -1
            else :
                score_pos_comp=parallel_len_score_pos_comp
            
            # first do parallel 
            for candidate in parallel_complete :
                lis=['+',candidate[-1][1],candidate[-1][0]] # the last two are initializer length and initializer seq
                # add a tag to describe HB pattern
                if parallel_hb_pattern[0]==(None,None) or parallel_hb_pattern[0]==[None,None] : hb_tag=0 # meaning that the first residue in the complementary is not involved in HB
                else: hb_tag=1 # meaning that the first residue in the complementary is involved in HB
                lis+=[hb_tag, 1]
                if debug : print('DEB: parallel candidate',candidate)
                seq=candidate[1]
                if HAVE_camsol :
                    if len(seq)>=7 :
                        _,solubility_score = camsol_intrinsic(seq, pH=pH, use_old_table=RUN_CAMSOL_PUBLISHED_VERSION)
                        solubility_bin=get_solubility_bin(solubility_score) # bin is 0 for <=-2 and grows by one every 0.5
                        extra_str=' (CamSol=%lf)' % (solubility_score)
                    else : 
                        solubility_score='too_short'
                        solubility_bin='too_short'
                        extra_str=''
                # save candidate composition
                if out_file_composition is not None :
                    relevant_indices= print_score_pos_comp(target_peptide, score_pos_comp, out_file=out_file_composition \
                                         ,interested_in_range=interested_in_range,template_hb_pattern=parallel_hb_pattern, mark_entries_in=seq,extra_str=extra_str, print_only_marked=True, marked_flag='par')
                    out_file_composition.write('\n')
                    out_file_composition.flush()
                    if  self.score_candidates :
                        score, count_score, promiscuity_score= score_candidate(score_pos_comp, relevant_indices)
                
                # below here the code for parallel and antiparallel should be precisely the same
                
                # save scores, which so far will be zero aside from solubility_score
                lis+=[score,solubility_score,solubility_bin, count_score, promiscuity_score]
                if seq not in complementary_dic :
                    complementary_dic[seq]= lis
                else :
                    if complementary_dic[seq][3]== hb_tag and complementary_dic[seq][0]=='+' : # compare direction and HB pattern
                        complementary_dic[seq][4]+=1 # add one to number of possible cascade paths
                        #self.warn_message("complement + %s found more than once, keeping only first found!!" % seq)
                    else : self.warn_message("complement + %s found more than once as a sequence, but present with differnt hb/direction, keeping only first found!!" % seq)
                del lis
          
        
        if self.antiparallel_fragments is not None and self.antiparallel_fragments!=[] :
            for k,out_tuple in enumerate(self.antiparallel_list_of_tup) :   
                anti_complete, anti_incomplete, anti_len_score_pos_comp,_,antiparallel_hb_pattern = out_tuple
                if interested_in_range is not None :
                    for ca in anti_incomplete :
                        if ca[0]<=start and ca[2]>=end :
                            anti_complete+=[ca]
                if anti_complete==[]:
                    self.warn_message("no complete anti parallel strand found for %s! in tup %d with hb_pattern %s\n    (consider targeting a longer or shorter bit of sequence or giving hb_tolerance_level=1)\n" % (target_peptide,k,str(antiparallel_hb_pattern)))
                
                if type(anti_len_score_pos_comp) is dict :
                    len_l=list(anti_len_score_pos_comp.keys())
                    len_l.sort(reverse=True)
                    score_pos_comp=[]
                    for l in len_l :
                        if l>0 :
                            score_pos_comp+=anti_len_score_pos_comp[l]
                        else :
                            sys.stderr.write("ERROR unrecognized key in anti_len_score_pos_comp (type is dict but key %s not recognized)\n" % (str(l)))
                            return -1
                else :
                    score_pos_comp=anti_len_score_pos_comp
                    
                # now do antiparallel
                for candidate in anti_complete :
                    if debug : print('DEB: anti candidate',candidate)
                    lis=['-',candidate[-1][1],candidate[-1][0]] # the last two are initializer length and initializer seq
                    # add a tag to describe HB pattern
                    if antiparallel_hb_pattern[0]==(None,None) or antiparallel_hb_pattern[0]==[None,None] : hb_tag=0 # meaning that the first residue in the complementary is not involved in HB
                    else: hb_tag=1 # meaning that the first residue in the complementary is involved in HB
                    lis+=[hb_tag,1]
                    seq=candidate[1]
                    if HAVE_camsol :
                        if len(seq)>=7 :
                            _,solubility_score = camsol_intrinsic(seq, pH=pH, use_old_table=RUN_CAMSOL_PUBLISHED_VERSION)
                            solubility_bin=get_solubility_bin(solubility_score) # bin is 0 for <=-2 and grows by one every 0.5
                            extra_str=' (CamSol=%lf)' % (solubility_score)
                        else : 
                            solubility_score='too_short'
                            solubility_bin='too_short'
                            extra_str=''
                        
                    if out_file_composition is not None : 
                        relevant_indices= print_score_pos_comp(target_peptide, score_pos_comp, out_file=out_file_composition \
                                             ,interested_in_range=interested_in_range,template_hb_pattern=antiparallel_hb_pattern, mark_entries_in=seq,extra_str=extra_str, print_only_marked=True, marked_flag='ant')
                        out_file_composition.write('\n')
                        out_file_composition.flush()
                        if  self.score_candidates :
                            score, count_score, promiscuity_score= score_candidate(score_pos_comp, relevant_indices)
                    
                    # below here the code for parallel and antiparallel should be precisely the same
                    
                    # save the entry
                    lis+=[score,solubility_score,solubility_bin, count_score, promiscuity_score]
                    if seq not in complementary_dic :
                        complementary_dic[seq]= lis
                    else : 
                        if complementary_dic[seq][3]== hb_tag and complementary_dic[seq][0]=='-' : # compare direction and HB pattern
                            complementary_dic[seq][4]+=1 # add one to number of possible cascade paths
                            #self.warn_message("complement - %s found more than once, keeping only first found!!" % seq)
                        else : self.warn_message("complement - %s found more than once as a sequence, but present with differnt hb/direction, keeping only first found!!" % seq)
                    del lis
                
        if close_out_file_composition : out_file_composition.close()
            # save final output
        hd={}
        for j,c in enumerate(complementary_hd) : hd[c]=j
        csv_dict.csvPrintDictionary(self.out_filename, complementary_dic,key_column=1,DELIMITER=None,HEADER=hd,key_column_header_name='comp:'+target_peptide)
        
        self.log_message("\nOutput printed in %s and peptide composition printed in %s\n" % (str(self.out_filename),str(self.out_file_composition)))
        if self.save_memory and not debug:
            del self.antiparallel_list_of_tup,self.parallel_list_of_tup
        return complementary_dic,hd
    
    def print_keys(self,key_list,outfile=sys.stdout,excel=False):
        """
        after make has been called it can be used to print in tab separated some keys belonging to the complementary_dic
        """
        if excel :
            outfile='tmp%d.tsv' % (numpy.random.randint(1000,9999))
        self.print_table(self.complementary_dic,only_these_keys=key_list,complementary_hd=self.complementary_hd,outfile=outfile,outfile_separator='\t')
        if excel :
            os.system('excel %s' % (outfile))
            os.system('sleep 5') # allow five seconds to open the file before deleting it
            os.system('rm -f %s' % (outfile))
        return
        
    def copy(self,another_Complement_class,copy_also_outfilenames=True):
        """
        it copies fragments database and other parameters from another class of the same kind
        """
        self.antiparallel_fragments=another_Complement_class.antiparallel_fragments
        self.parallel_fragments=another_Complement_class.parallel_fragments
        self.antiparallel_enrichment_scores= another_Complement_class.antiparallel_enrichment_scores
        self.parallel_enrichment_scores =another_Complement_class.parallel_enrichment_scores
        self.use_hb_patterns=another_Complement_class.use_hb_patterns
        self.folder_name_with_frag_database=another_Complement_class.folder_name_with_frag_database
        self.fragment_length_lower_threshold=another_Complement_class.fragment_length_lower_threshold
        self.score_candidates=another_Complement_class.score_candidates
        self.hb_tolerance_level=another_Complement_class.hb_tolerance_level
        self.limit_new_candidates=another_Complement_class.limit_new_candidates
        self.warning_file=another_Complement_class.warning_file
        self.logfile=another_Complement_class.logfile
        self.debug=another_Complement_class.debug
        self.consider_extra_initializers=another_Complement_class.consider_extra_initializers
        self.consider_sub_fragments_in_scoring= another_Complement_class.consider_sub_fragments_in_scoring
        self.save_memory=another_Complement_class.save_memory
        self.initializer_minimum_length=another_Complement_class.initializer_minimum_length
        if copy_also_outfilenames :
            self.out_file_composition=another_Complement_class.out_file_composition
            self.out_filename=another_Complement_class.out_filename
        return
            




def screen_complementary_dic(complementary_dic_or_class, condition= lambda x,cd,hd : cd[x][hd['len_initializer']]>4 and cd[x][hd['Solubility_score']]>-0.9 ,consider_only_keys=None, debug=True,hd=None):
    """
    example, d is a complementary class where make has been called:
m=beta_strands.screen_complementary_dic(d,lambda x,cd,hd : cd[x][hd['len_initializer']]>4 and cd[x][hd['Solubility_score']]>-0.6 and cd[x][hd['FirstResHB']]==1)
    otherwise d can be the complementary dict (also loaded from file)
    Example:
ab=csv_dict.Data()
ab.load('abetaC_complementary_peptide_coverage.tsv')
keys=beta_strands.screen_complementary_dic(ab, condition= lambda x,cd,hd : cd[x][hd['binds_from_top']]==True and cd[x][hd['parallel/antiparallel']] =='-' and cd[x][hd['Target Sequence']][ cd[x][hd['Complementary_sequence']].find(cd[x][hd['initializer']]):cd[x][hd['Complementary_sequence']].find(cd[x][hd['initializer']])+4 ]=='IGLM' )
    
    """
    if type(complementary_dic_or_class) is not dict and not isinstance(complementary_dic_or_class, csv_dict.Data):#and not isinstance(complementary_dic_or_class, OrderedDict) :
        cd=complementary_dic_or_class.complementary_dic
        hd=complementary_dic_or_class.complementary_hd
    else :
        cd=complementary_dic_or_class
        if hd  is None :
            if isinstance(complementary_dic_or_class, csv_dict.Data) : hd=complementary_dic_or_class.hd
            else : raise Exception("Error in screen_complementary_dic(): If you don't give a class as complementary_dic_or_class then you must give the kwarg hd (the third arg) (now hd=None)\n")
    if consider_only_keys is None : consider_only_keys=list(cd.keys())
    if condition is None :
        return consider_only_keys
    keys_matching_condition=[]
    failed=0
    for x in consider_only_keys :
        try :
            if condition(x,cd,hd) :
                keys_matching_condition+=[x]
        except Exception :
            failed+=1
    if failed>0 : 
        sys.stderr.write("**WARNING** assessment of condition yielded exception in %d cases\n" % (failed))
    if debug : print('keys_matching_condition: ',len(keys_matching_condition),list(hd.keys()))
    return keys_matching_condition

def sort_complementary_dic(complementary_dic_or_class,consider_only_keys=None, sort_keys=['Complementarity_score','Solubility_bin'],reverse=[True,True],hd=None):
    '''
    entries in sort_keys are used to sort the  complementary_dic_or_class in the order specified by entries in reverse.
     sort_keys are applied from first to last, so that entries are ultimately ordered according to sort_keys[-1] in the order reverse[-1]
    it does not sort the dictionary (as this could in principle be a dict rather than an OrderedDict) but it 
    return sorted_keys
    '''
    if type(complementary_dic_or_class) is not dict and not isinstance(complementary_dic_or_class, csv_dict.Data) and not isinstance(complementary_dic_or_class, OrderedDict) :
        cd=complementary_dic_or_class.complementary_dic
        hd=complementary_dic_or_class.complementary_hd
    else :
        cd=complementary_dic_or_class
        if hd is None :
            if isinstance(complementary_dic_or_class, csv_dict.Data) : hd=complementary_dic_or_class.hd
            else : raise Exception("Error in screen_complementary_dic(): If you don't give a class as complementary_dic_or_class then you must give the kwarg hd (the third arg) (now hd=None)\n")
    if consider_only_keys is None : consider_only_keys=list(cd.keys())
    sorted_keys=[ [k]+[complementary_dic_or_class[k][hd[sk]] for sk in sort_keys ] for k in consider_only_keys ]
    if len(sorted_keys)<1 :
        if consider_only_keys is not None : sys.stderr.write("potential **WARNING** in sort_complementary_dic() len(sorted_keys)<1 len(complementary_dic_or_class)=%d len(consider_only_keys)=%d\n" % (len(cd),len(consider_only_keys)))
        else : sys.stderr.write("potential **WARNING** in sort_complementary_dic() len(sorted_keys)<1 len(complementary_dic_or_class)=%d type(consider_only_keys)=%s\n" % (len(cd),str(type(consider_only_keys))))
        return []    
    for j,sk in enumerate(sort_keys) :
        if reverse[j] : f=-1.
        else : f=1.
        sorted_keys.sort( key= lambda x : f*x[j+1] )
    return list(zip(*sorted_keys))[0]

def sort_complementary_dic_each_initialiser(complementary_dic_or_class,consider_only_keys=None, sort_keys=['Complementarity_score','Solubility_bin'],reverse=[True,True],hd=None):
    '''
    entries in sort_keys are used to sort the  complementary_dic_or_class in the order specified by entries in reverse.
     sort_keys are applied from first to last, so that entries are ultimately ordered according to sort_keys[-1] in the order reverse[-1]
    it does not sort the dictionary (as this could in principle be a dict rather than an OrderedDict) but it 
    return initialiser_sorted_keys # a dict with the different initialisers and its correspodning sorted keys
    '''
    if type(complementary_dic_or_class) is not dict and not isinstance(complementary_dic_or_class, csv_dict.Data) and not isinstance(complementary_dic_or_class, OrderedDict) :
        cd=complementary_dic_or_class.complementary_dic
        hd=complementary_dic_or_class.complementary_hd
    else :
        cd=complementary_dic_or_class
        if hd is None :
            if isinstance(complementary_dic_or_class, csv_dict.Data) : hd=complementary_dic_or_class.hd
            else : raise Exception("Error in screen_complementary_dic(): If you don't give a class as complementary_dic_or_class then you must give the kwarg hd (the third arg) (now hd=None)\n")
    if consider_only_keys is None : consider_only_keys=list(cd.keys())
    
    if len(consider_only_keys)<1 : 
        sys.stderr.write("potential **WARNING** in sort_complementary_dic_each_initialiser() len(sorted_keys)<1 len(complementary_dic_or_class)=%d len(consider_only_keys)=%d\n" % (len(cd),len(consider_only_keys)))
        return {}
    
    initialiser_sorted_keys={}
    for k in consider_only_keys :
        initialiser= ( complementary_dic_or_class[k][hd['initializer']], complementary_dic_or_class[k][hd['parallel/antiparallel']],complementary_dic_or_class[k][hd['binds_from_top']], complementary_dic_or_class[k][hd['first residue index']]+(complementary_dic_or_class[k][hd['Complementary_sequence']]).find(complementary_dic_or_class[k][hd['initializer']]) ) # we should also save the HB pattern of the initialiser but this is not currently saved in the class..
        if initialiser not in initialiser_sorted_keys : initialiser_sorted_keys[initialiser]=[]
        initialiser_sorted_keys[initialiser]+=[k]
    sorted_init_keys={}
    for init in initialiser_sorted_keys :
        sorted_keys=[ [k]+[complementary_dic_or_class[k][hd[sk]] for sk in sort_keys ] for k in initialiser_sorted_keys[init] ]
        if len(sorted_keys)<1 :
            continue    
        for j,sk in enumerate(sort_keys) :
            if reverse[j] : f=-1.
            else : f=1.
            sorted_keys.sort( key= lambda x : f*x[j+1] )
        sorted_init_keys[init]=list(zip(*sorted_keys))[0]
    return sorted_init_keys


def screen_complementary_peptides(protein, comp_len=[7,8],every=1,start_rescount=1,initializer_minimum_length=4, condition= lambda x,cd,hd : cd[x][hd['len_initializer']]>4 and cd[x][hd['Solubility_score']]>-0.5, allow_extra_res_in_range=True,assign_two_loop_pairs_offset=3,out_file_composition=None, copy=None,outpickle=None,debug=False ,save_coverage_results=None, score_candidates=True,logfile=sys.stdout,warn_file=None, just_save_coverage_results=False):
    '''
    given a protein sequence it generate all possible complementary fragmetnts of length comp_len
    starting at every every residue from the beginning of the protein. 
    copy can be used to init the Complement() class from another class, so that the database is linked with a pointer (super fast) rather then reloaded (super slow)
    it returns parallel_hb_0, parallel_hb_1 , anti_hb_0 , anti_hb_1
     which are four dictionaries, with keys that are indices along the protein sequence of where the complementary peptide begins, 
      note that some positions may be present as key, but the corresponding value may be []. This means that althought there is no complementary 
      in e.g. anti_hb_0 for that position, there is at least one complementary in at least one of the other 3 dictionaries (anti_hb_1,parallel_hb_0, parallel_hb_1).
      positions not present in keys are positions for which NO complementary has been found in EVERY of the four dictionaries. 
     values are lists of all the possible complements found at that position.
    the four dictionaries contain either parallel or antiparallel complementary peptides and 
    hb_0 means that the first residue of that peptide is NOT involved in Hydrogen bond, 1 means that it is.
    if not just_save_coverage_results
     inner list contain 'Complementary_sequence','parallel/antiparallel','len_initializer','initializer'
     and if score_candidates also 'Complementarity_score', 'Solubility_score'
    assign_two_loop_pairs_offset is the offset value to group together complementary peptides amenable for grafting onto a two-loop scaffold.
      set it to None if you do not wish to perform this calculation. offset is the N-terminal distance (in amino acids) between the N-termini of the complementary peptides
      it assigns an integer number to each peptide so that pairs with consecutive numbers (e.g. 1,2 or 5,6) can be used together for double loop design.
      At the moment it does not perform mixed parallel antiparallel design, only pepetides with the same direction will be used together.
      
parallel_hb_0, parallel_hb_1 , anti_hb_0 , anti_hb_1= beta_strands.protein_coverage(Nasyn,8,copy=C)
_=beta_strands.visualize_coverage( Nasyn,(parallel_hb_0, parallel_hb_1 , anti_hb_0 , anti_hb_1), outfile=sys.stdout)
    '''
    #if copy is not None : copy.logfile=None
    c=Complement(copy=copy,initializer_minimum_length=initializer_minimum_length,consider_extra_initializers=12)
    c.score_candidates=score_candidates
    c.initializer_minimum_length=initializer_minimum_length
    c.logfile=None
    c.warning_file=None #sys.stderr
    c.out_file_composition=None
    c.out_filename=None
    #print c.logfile
    close_out_file_composition=False
    if out_file_composition is not None :
        if type(out_file_composition) is str : 
            out_file_composition=open(out_file_composition,'w')
            close_out_file_composition=True     
        c.get_compositon_dic=True
    if not hasattr(comp_len, '__len__') : comp_len=[comp_len]
    tot_matching=0 # tot number of complementary peptides found in the protein that match the condition.
    parallel_hb_0={}
    parallel_hb_1={}
    anti_hb_0={}
    anti_hb_1={}
    if not just_save_coverage_results :
        for i in range(len(protein)-min(comp_len)+1) :
            parallel_hb_0[i]=[]
            parallel_hb_1[i]=[]
            anti_hb_0[i]=[]
            anti_hb_1[i]=[]
    if type(save_coverage_results) is bool and save_coverage_results==False : save_coverage_results=None
    if save_coverage_results is not None : 
        protein_coverage_res=csv_dict.Data()
        rhd=['first residue index','first residue position','len_complementary_seq' ,'Target Sequence', 'parallel/antiparallel','binds_from_top','FirstResHB','len_protein','Complementary_sequence','len_initializer','initializer','NumberOfPaths','Solubility_score','Solubility_bin']
        if score_candidates: rhd+=['Complementarity_score','Count_score','Promiscuity_score','Mean_enrichment']
        protein_coverage_res.hd=csv_dict.list_to_dictionary(rhd)
        tot=0
    interested_in_range=None # it changes after if required
    for c_len in comp_len :
        interval=list(range( 0 , len(protein)-c_len+1 , every))
        N=1.*len(interval)
        if logfile is not None : 
            logfile.write('\nlooking for complementary peptides of length %d\n' % (c_len))
        for j,i in  enumerate(interval): # interval is an index that runs from 0 (start of target sequence) till len(target_sequence)-c_len!
            if logfile is not None and j%10==0 :
                logfile.write(' pos:%d -->%5.2lf %% |' % (i,100.*j/(N)))
                logfile.flush()
            c.complementary_dic={}
            c.composition_dic={}
            if allow_extra_res_in_range  is not None :
                if i>0 and i+c_len < len(protein) : 
                    target=protein[i-1:i+c_len+1]
                    interested_in_range=[1,len(target)-1]
                    #print 'MID',i,interested_in_range,target,c_len,len(protein),i+c_len,j
                elif i>0 :
                    target=protein[i-1:i+c_len+1]
                    interested_in_range=[1,len(target)]
                    #print 'END',i,interested_in_range,target
                else :
                    target=protein[i:i+c_len+1]
                    interested_in_range=[0,len(target)-1]
                    #print '\nSTA',i,interested_in_range,target
            else : target=protein[i:i+c_len]
            c.interested_in_range=interested_in_range
            c.make(target ) # run the cascade method
            #d = beta_strands.Complement(target, copy=C, interested_in_range=interested_in_range,initializer_minimum_length=4,consider_extra_initializers=12, out_file_composition='complementary_composition_'+out_tag+'.txt',out_filename='complementary_peptides_'+out_tag+'.tsv')
            results=c.complementary_dic
            hd=c.complementary_hd
            if len(results)==0 :
                if warn_file is not None : warn_file.write( '****** NO COMPLEMENTS found at index=%d for segment %s\n' % (i,protein[i:i+c_len]))
                del results
                continue
            if c.interested_in_range is not None : target_seq=target[ c.interested_in_range[0]:c.interested_in_range[1]]
            else : target_seq=target
            keys_matching_condition = screen_complementary_dic(c, condition,debug=False)
            if debug : 
                print('\n',i,target,interested_in_range,target_seq,'len(keys_matching_condition)',len(keys_matching_condition),'len(results)',len(results),type(list(results.values())[0] ),list(results.values())[0],hd)
                #for jkb,kbk in enumerate(keys_matching_condition) :
                #    if 'YVKIIETSM' in kbk[0] : print '----------------->>',jkb,kbk
            if len(keys_matching_condition)==0 : continue
            tot_matching+=len(keys_matching_condition)
            
            for k in keys_matching_condition :
                if save_coverage_results is not None :
                    # hd should contain: 
                    # complementary_hd=['CandidateIndex', 'Complementary_sequence','parallel/antiparallel','len_initializer','initializer','FirstResHB','NumberOfPaths','Complementarity_score','Solubility_score','Count_score','Promiscuity_score' ]
                    tot+=1
                    binds_from_top=1 # silly index that assumes the first residue of the target protein has HB to the top, the second to the bottom and so on.
                    if results[k][hd['parallel/antiparallel']]=='+' : 
                        if i%2==0 and results[k][hd['FirstResHB']] : binds_from_top=0 #  the facing aa of a parallel complementary with HB=1 will NOT be involved in HB.
                        elif i%2==1 and not results[k][hd['FirstResHB']] : binds_from_top=0
                    else : # is antiparallel 
                        if i%2==1 and results[k][hd['FirstResHB']] : binds_from_top=0  # the facing aa of a antiparallel complementary with HB=1 will  be involved in HB.
                        elif i%2==0 and not results[k][hd['FirstResHB']] : binds_from_top=0
                    protein_coverage_res[tot]= [ i,i+start_rescount ,c_len, target_seq , results[k][hd['parallel/antiparallel']],binds_from_top,results[k][hd['FirstResHB']],len(protein),results[k][hd['Complementary_sequence']],results[k][hd['len_initializer']],results[k][hd['initializer']],results[k][hd['NumberOfPaths']],results[k][hd['Solubility_score']],results[k][hd['Solubility_bin']] ]
                    if score_candidates : protein_coverage_res[tot] += [ results[k][hd['Complementarity_score']],results[k][hd['Count_score']],results[k][hd['Promiscuity_score']],results[k][hd['Mean_enrichment']] ]
                if not just_save_coverage_results : # the order below is quite important (I think I don't remember exactly..)
                    to_save=[ results[k][hd['Complementary_sequence']],results[k][hd['parallel/antiparallel']],results[k][hd['len_initializer']],results[k][hd['initializer']] ]  #[ results[k][1:hd['FirstResHB']] ] # 'parallel/antiparallel','len_initializer','initializer'
                    if score_candidates : to_save += [ results[k][hd['Complementarity_score']],results[k][hd['Solubility_score']]  ]
                    if results[k][hd['parallel/antiparallel']]=='+' :
                        if results[k][hd['FirstResHB']]==0 : parallel_hb_0[i]+= [ to_save ]
                        else : parallel_hb_1[i]+= [ to_save ]
                    else :
                        if results[k][hd['FirstResHB']]==0 : anti_hb_0[i]+=[ to_save ]
                        else : anti_hb_1[i]+= [ to_save ]
                if out_file_composition is not None :
                    out_file_composition.write( ' index_of_first_seq_res: %d position %d len=%d\n%s\n' % (i,i+start_rescount,len(results[k][hd['Complementary_sequence']]),c.composition_dic[k]))
            del results
    if assign_two_loop_pairs_offset  is not None :
        protein_coverage_res= assign_two_loop_pairs(protein_coverage_res,offset=assign_two_loop_pairs_offset, index_key='first residue index',hb_key='FirstResHB',direction_key='parallel/antiparallel')
    print('\n %d complementary peptides passing the conditon found for the input protein' % (tot_matching))
    if logfile is not None : logfile.write('\n\n')
    if close_out_file_composition :
        out_file_composition.close()
    if save_coverage_results is not None :
        if type(save_coverage_results) is str : outname=save_coverage_results
        else : outname='protein_coverage_results'+protein[:8]+'etc.tsv'
        protein_coverage_res.Print(outname)
        print('coverage results saved in %s' % (outname))
    if not just_save_coverage_results :
        if outpickle  is not None :
            out=open(outpickle,'wb')
            pickle.dump( (parallel_hb_0, parallel_hb_1 , anti_hb_0 , anti_hb_1),out)
            out.close()
        return parallel_hb_0, parallel_hb_1 , anti_hb_0 , anti_hb_1





def assign_two_loop_pairs(coverage_res, offset=3,index_key='first residue index',hb_key='FirstResHB',peptide_sequence_key='Complementary_sequence',direction_key='parallel/antiparallel',offset_from_N_term_of_peptides=True):
    '''
CAREFUL:
1-check if index_key accounts for interest_in_range (peptide_sequence_key does).
    coverage_res is a Data() class printed for instance by screen_complementary_peptides()
    offset is the N-terminal distance (in amino acids) between the N-termini of the complementary peptides
    it assigns an integer number to each peptide so that pairs with consecutive numbers (e.g. 1,2 or 5,6) can be used together for double loop design.
    At the moment it does not perform mixed parallel antiparallel design, only pepetides with the same direction will be used together.
    if offset_from_N_term_of_peptides is False the offset is calculated from the N terminus of the target sequence (thus in the case of antiparallel peptides from their C terminus).
rhd=['first residue index','first residue position','len_complementary_seq' ,'Target Sequence', 'parallel/antiparallel','binds_from_top','FirstResHB','len_protein','Complementary_sequence','len_initializer','initializer','NumberOfPaths','Solubility_score']
if score_candidates: rhd+=['Complementarity_score','Count_score','Promiscuity_score','Mean_enrichment']
    '''
    if type(coverage_res) is str :
        fname=coverage_res
        coverage_res=csv_dict.Data()
        coverage_res.load(fname)
    if len(coverage_res)==0 : return coverage_res
    pairs_pks=[]
    connected_pks=[]
    pks_inds={}
    coverage_res = coverage_res.sort(index_key)
    candidate_pairs={}
    pair_count=0
    expected_len=len(list(coverage_res.values())[0])
    for k in coverage_res :
        r_ind=coverage_res[k][coverage_res.hd[index_key]]
        d=coverage_res[k][coverage_res.hd[direction_key]]
        r_hb=coverage_res[k][coverage_res.hd[hb_key]]
        if offset_from_N_term_of_peptides and d=='-' : # in this case we must save in the key the C terminus index of the target sequence (that would correspond to the N-terminus end of the peptide)
            l=len( coverage_res[k][coverage_res.hd[peptide_sequence_key]])
            if l%2==0 : r_hb= not r_hb
            r_ind= r_ind + l -1
        pk=( r_ind, d, r_hb )
        if pk not in candidate_pairs : candidate_pairs[pk]=[]
        candidate_pairs[pk]+=[ k ]
        
        p_r_ind = r_ind-offset # this is the index on the target sequence of peptides potentially suitable to match k in a two-loop design
        if p_r_ind>=0 :
            if offset%2==0 : p_hb= not r_hb # calculate the suitable hb pattern for potentially matching peptides
            else : p_hb= r_hb
            pk2=(p_r_ind, d ,p_hb)
            if pk2 in candidate_pairs :
                if pk2 in pks_inds : 
                    if pk not in connected_pks[pks_inds[pk2]] :
                        connected_pks[pks_inds[pk2]]+=[pk]
                        pks_inds[pk]=pks_inds[pk2]
                elif pk in pks_inds : # should actually never happen
                    if pk2 not in connected_pks[pks_inds[pk]] :
                        connected_pks[pks_inds[pk]]+=[pk2]
                        pks_inds[pk2]=pks_inds[pk]
                else :
                    pks_inds[pk2]=pks_inds[pk]=len(connected_pks)
                    connected_pks+=[  [pk2,pk] ]
                if [pk2,pk] not in pairs_pks : pairs_pks+=[ [pk2,pk] ]
    for group in connected_pks :
        for pk in group :
            for k in candidate_pairs[pk] :
                coverage_res[k]+=[pair_count]
            pair_count+=1
        pair_count+=1 # so that it gets disconnected
        '''
        if pk2 not in done_pks : 
            pair_count+=3
            done_pks[pk2]=pair_count
        for k2 in candidate_pairs[pk2] :
            #pairs+=[ [k2,k] ]
            if len(coverage_res[k2])!=expected_len+1 : coverage_res[k2]+=[ done_pks[pk2] ]
            elif coverage_res[k2][-1]!=done_pks[pk2] : sys.stderr.write("**WARNING** in assign_two_loop_pairs() when pairing %s we have pair_count already assinged at %s and trying to assign %s - leaving unchanged [len %d expected_len %d]\n" % ( str([k2,k]),str(coverage_res[k2][-1]),str(pair_count),len(coverage_res[k2]),expected_len))
        coverage_res[k]+=[ done_pks[pk2]+1 ]
        if pk not in done_pks :done_pks[pk]=done_pks[pk2]+1 
        '''
    # equalise to same length by adding an extra 'n' to non matched peptides (i.e. those that don't have a suitable two-loop partner)
    for k in coverage_res :
        if len(coverage_res[k])<expected_len+1 : coverage_res[k]+=['n']
    coverage_res._update_hd('two_loop_pairs_off_%d' % (offset))
    return coverage_res



# toDebug: epitope,offset_Nterm_target,peptide1,peptide2,initialiser1,initialiser2,first_hb1,first_hb2='DAEFRHDS',1,'DIRYVSRTLSG','DVVQAGS','DIRY','DVVQ',False,False
def _two_loops_pretty_print(epitope,offset_Nterm_target,peptide1,peptide2,initialiser1,initialiser2,first_hb1,first_hb2):
    j=peptide1.find(initialiser1)
    pretty1=''
    for i in range(len(peptide1)) :
        if i in range(j,j+len(initialiser1)) : pretty1+=peptide1[i]
        else : pretty1+=peptide1[i].lower()
    j=peptide2.find(initialiser2)
    pretty1+=' '*(len(epitope)-len(peptide1))
    pretty2=' '*offset_Nterm_target
    for i in range(len(peptide2)) :
        if i in range(j,j+len(initialiser2)) : pretty2+=peptide2[i]
        else : pretty2+=peptide2[i].lower()
    if first_hb1 is not None :
        if first_hb1 : hb1=(': '*len(epitope))[:len(epitope)]
        else : hb1=(' :'*len(epitope))[:len(epitope)]
        if first_hb2 : hb2=(' '*offset_Nterm_target)+(': '*len(epitope))[:len(epitope)-offset_Nterm_target]
        else : hb2=(' '*offset_Nterm_target)+(' :'*len(epitope))[:len(epitope)-offset_Nterm_target]
        line= '%s\n%s\n%s\n%s\n%s' % (pretty1,hb1,epitope,hb2,pretty2)
    else :
        line='%s\n%s\n%s' % (pretty1,epitope,pretty2)
    return line
    

def _merge_two_loop_peptides(coverage_res, k, k1 ,pair_reference='',merge_only_lengths_first_second_tup=None, add_pretty_line_csv=True, compostion_index_dict=None) :
    '''
    aux function used by merge_two_loop_pairs. it implements the conditions used to merge (average, min or max) two peptides and their scores.
    rhd=['first residue index','first residue position','len_complementary_seq' ,'Target Sequence', 'parallel/antiparallel','binds_from_top','FirstResHB','len_protein','Complementary_sequence','len_initializer','initializer','NumberOfPaths','Solubility_score']
    rhd+=['Complementarity_score','Count_score','Promiscuity_score','Mean_enrichment']
    '''
    if coverage_res[k][coverage_res.hd['binds_from_top']]== coverage_res[k1][coverage_res.hd['binds_from_top']] :
        sys.stderr.write('WARNING in _merge_two_loop_peptides peptides to merge both have binds_from_top= %s [k=%s, k1=%s]\n' % (str(coverage_res[k1][coverage_res.hd['binds_from_top']]),str(k),str(k1)))
    if coverage_res[k][coverage_res.hd['first residue index']] <= coverage_res[k1][coverage_res.hd['first residue index']] : 
        first_k,second_k=k,k1
    else : first_k , second_k = k1 , k
    first_peptide,second_peptide = coverage_res[first_k][coverage_res.hd['Complementary_sequence']],coverage_res[second_k][coverage_res.hd['Complementary_sequence']]
    if merge_only_lengths_first_second_tup is not None :
        if merge_only_lengths_first_second_tup[0] is not None and len(first_peptide)!= merge_only_lengths_first_second_tup[0] : 
            return None,None
        if merge_only_lengths_first_second_tup[1] is not None and len(second_peptide)!= merge_only_lengths_first_second_tup[1] :
            return None,None 
    offset_Nterm_target=(coverage_res[second_k][coverage_res.hd['first residue index']]-coverage_res[first_k][coverage_res.hd['first residue index']])
    target_sequence=  coverage_res[first_k][coverage_res.hd['Target Sequence']][ :offset_Nterm_target] + coverage_res[second_k][coverage_res.hd['Target Sequence']]
    if len(second_peptide)+offset_Nterm_target<len(first_peptide) :
        target_sequence +=  coverage_res[first_k][coverage_res.hd['Target Sequence']][len(target_sequence):]
    start_ind= min( coverage_res[k][coverage_res.hd['first residue index']], coverage_res[k1][coverage_res.hd['first residue index']] )
    start_pos= min( coverage_res[k][coverage_res.hd['first residue position']], coverage_res[k1][coverage_res.hd['first residue position']] )
    epitope_range='%03d;%03d' % (start_pos,start_pos+len(target_sequence)-1)
    
    first_peptide_d_hb,second_peptide_d_hb = coverage_res[first_k][coverage_res.hd['parallel/antiparallel']]+';'+str(coverage_res[first_k][coverage_res.hd['FirstResHB']]) , coverage_res[second_k][coverage_res.hd['parallel/antiparallel']]+';'+str(coverage_res[second_k][coverage_res.hd['FirstResHB']])
    peptide_lengths='%02d;%02d' % (len(first_peptide),len(second_peptide))
    init1,init2= coverage_res[first_k][coverage_res.hd['initializer']],coverage_res[second_k][coverage_res.hd['initializer']]
    
    # get the min for solubility
    solubility= min( coverage_res[k][coverage_res.hd['Solubility_score']], coverage_res[k1][coverage_res.hd['Solubility_score']] )
    # and the mean for the other scores
    mean_init_len=numpy.mean([ coverage_res[k][coverage_res.hd['len_initializer']], coverage_res[k1][coverage_res.hd['len_initializer']] ])
    mean_NumberOfPaths =numpy.mean([ coverage_res[k][coverage_res.hd['NumberOfPaths']], coverage_res[k1][coverage_res.hd['NumberOfPaths']] ])
    if 'Mean_enrichment' in coverage_res.hd : me = numpy.mean([ coverage_res[k][coverage_res.hd['Mean_enrichment']], coverage_res[k1][coverage_res.hd['Mean_enrichment']] ])
    if 'Count_score' in coverage_res.hd : 
        countscore = numpy.mean([ coverage_res[k][coverage_res.hd['Count_score']], coverage_res[k1][coverage_res.hd['Count_score']] ])
        promiscuity= numpy.mean([ coverage_res[k][coverage_res.hd['Promiscuity_score']], coverage_res[k1][coverage_res.hd['Promiscuity_score']] ])
    # retrieve and edit composition file
    twoloop_composition_lines=None
    if compostion_index_dict is not None :
        ind_k_1= ( coverage_res[first_k ][coverage_res.hd['first residue index']],   coverage_res[first_k ][coverage_res.hd['len_complementary_seq']] )
        ind_k_2= ( coverage_res[second_k][coverage_res.hd['first residue index']],   coverage_res[second_k][coverage_res.hd['len_complementary_seq']])
        key_pep_1=(first_peptide, coverage_res[first_k ][coverage_res.hd['parallel/antiparallel']] , coverage_res[first_k ][coverage_res.hd['FirstResHB']] )
        key_pep_2=(second_peptide,coverage_res[second_k][coverage_res.hd['parallel/antiparallel']] , coverage_res[second_k][coverage_res.hd['FirstResHB']] )
        try :
            str_1 = compostion_index_dict[ind_k_1][key_pep_1]
            str_2 = compostion_index_dict[ind_k_2][key_pep_2]
            remove_first_double1,firstbar1,remove_last_double1,lastbar1 = find_remove_range_from_composition_line(str_1[0])
            remove_first_double2,firstbar2,remove_last_double2,lastbar2 = find_remove_range_from_composition_line(str_2[0])
            twoloop_composition_lines=[]
            for l in str_1[::-1] : # we also flip these
                twoloop_composition_lines += [ remove_range_from_composition_line(l,remove_first_double1,firstbar1,remove_last_double1,lastbar1) ]
                
            for l in str_2 : # here we add offset
                if ind_k_1[0]==0 : twoloop_composition_lines += [ (' '*(offset_Nterm_target-int(remove_first_double2))*2)+ remove_range_from_composition_line(l,remove_first_double2,firstbar2,remove_last_double2,lastbar2) ]
                else : twoloop_composition_lines += [ (' '*(offset_Nterm_target)*2)+ remove_range_from_composition_line(l,remove_first_double2,firstbar2,remove_last_double2,lastbar2) ]
            
        except Exception :
            sys.stderr.write("\n")
            print_exc( file=sys.stderr)
            sys.stderr.write("*** raised in _merge_two_loop_peptides() for ind_k_1=%s  ind_k_2=%s key_pep_1=%s key_pep_2=%s\n\n" % (str(ind_k_1),str(ind_k_2), str(key_pep_1),str(key_pep_2)))
            pass
    if add_pretty_line_csv :
        pretty_line = _two_loops_pretty_print(target_sequence,offset_Nterm_target,first_peptide,second_peptide,init1,init2,coverage_res[first_k][coverage_res.hd['FirstResHB']],coverage_res[second_k][coverage_res.hd['FirstResHB']])
        return [ start_ind, start_pos, target_sequence,epitope_range, pair_reference,peptide_lengths ,pretty_line,mean_init_len, solubility,countscore,promiscuity,me,mean_NumberOfPaths, first_peptide, second_peptide,first_peptide_d_hb,second_peptide_d_hb, coverage_res[first_k][coverage_res.hd['initializer']],coverage_res[second_k][coverage_res.hd['initializer']] ],  twoloop_composition_lines
    return [ start_ind, start_pos, target_sequence,epitope_range, pair_reference,peptide_lengths, first_peptide, second_peptide,first_peptide_d_hb,second_peptide_d_hb ,mean_init_len, solubility,countscore,promiscuity,me,mean_NumberOfPaths, coverage_res[first_k][coverage_res.hd['initializer']],coverage_res[second_k][coverage_res.hd['initializer']] ],  twoloop_composition_lines


def find_remove_range_from_composition_line(comp_line):
    remove_last_double=False
    lastbar = comp_line.rfind('|')
    if comp_line[lastbar-1]=='|' : remove_last_double=True
    remove_first_double=False
    firstbar = comp_line.find('|')
    if comp_line[firstbar+1]=='|' : remove_first_double=True
    return remove_first_double,firstbar,remove_last_double,lastbar
def remove_range_from_composition_line(comp_line,remove_first_double,firstbar,remove_last_double,lastbar):
    if remove_last_double : comp_line=comp_line[:lastbar]+comp_line[lastbar+1:]
    if remove_first_double : comp_line=comp_line[:firstbar]+comp_line[firstbar+1:]
    if comp_line[0]=='|': comp_line=comp_line[1:]
    return comp_line

def merge_two_loop_pairs(coverage_res, composition_file=None,merge_only_lengths_first_second_tup=None,merge_only_best=2, sample_all_initialisers=True,   condition= lambda x,cd,hd : cd[x][hd['len_initializer']]>=4 and cd[x][hd['Solubility_score']]>=0.7, sort_keys=['Complementarity_score','Solubility_bin'],reverse=[True,True], pairs_key='two_loop_pairs_off_3',outfile_composition='two_loop_composition.txt',add_pretty_line_csv=True,outfile_coverage='two_loop_coverage.csv'):
    '''
    exploiting the tag that create_two_loop_pairs() has added ( pairs_key usually in the format 'two_loop_pairs_off_%d' % (offset) )
    this function creates a single table with average (or min/max) score reflecting the quality of both peptides.
      this is a combinatorial problem so if I have N complementaries with pairs_key 0 and M with pairs_key 1
      I have N*M different possible two-loop constructs
    to limit this issue one can give non None merge_only_best, which is the number of peptides to merge for each pairs_key.
      So if merge_only_best=5 you will generate <=25 two-loop constructs.
      To determine the best the entries are sorted according to sort_keys
       sort_keys are applied from first to last, so that entries are ultimately ordered according to sort_keys[-1] in the order reverse[-1]
    merge_lengths_first_second_tup can be used to select only peptides of the specified length, it can be a tuple with length of first peptide (the one closer to N terminus of epitope)
       and second peptide. Each of these can be None if only the lenght of one peptide should be constrained.
    Similarly condition can be given to neglect all peptides not matching the given conditon
    composition_file is an input file, if given it prints a corresponding output-composition file with the two-loop compositions.
rhd=['first residue index','first residue position','len_complementary_seq' ,'Target Sequence', 'parallel/antiparallel','binds_from_top','FirstResHB','len_protein','Complementary_sequence','len_initializer','initializer','NumberOfPaths','Solubility_score']
if score_candidates: rhd+=['Complementarity_score','Count_score','Promiscuity_score','Mean_enrichment']
    '''
    if type(coverage_res) is str :
        fname=coverage_res
        coverage_res=csv_dict.Data()
        coverage_res.load(fname)
    compostion_index_dict=None
    if composition_file is not None :
        if type(composition_file) is str : compostion_index_dict = read_composition_file(composition_file,remove_range_from_key=True)
        else : compostion_index_dict = composition_file # assumes it has been already loaded
        outfile_comp=open(outfile_composition,'w')
    uniq=[]
    map_k=OrderedDict()
    for key in coverage_res :
        p=coverage_res[key][coverage_res.hd[pairs_key]]
        if type(p) is not int : continue
        if p not in uniq : 
            uniq+=[p]
            map_k[p]=[]
        map_k[p]+=[key]
    twol=csv_dict.Data()
    if add_pretty_line_csv :
        rhd=['start residue index','start residue position','Target Sequence' ,'epitope_range','pair_reference','peptide_lengths','binding_mode','mean_len_initializer','Min_Solubility_score','mean_Count_score','mean_Promiscuity_score','Mean_enrichment','NumberOfPaths','first_peptide','second_peptide','first_pep_info','second_pep_info','first_initializer','second_initializer']
    else :
        rhd=['start residue index','start residue position','Target Sequence' ,'epitope_range','pair_reference','peptide_lengths','first_peptide','second_peptide','first_pep_info','second_pep_info','mean_len_initializer','Min_Solubility_score','mean_Count_score','mean_Promiscuity_score','Mean_enrichment','NumberOfPaths','first_initializer','second_initializer']
    twol.hd=csv_dict.list_to_dictionary(rhd)
    uniq.sort()
    n=0
    for j,p in enumerate(uniq) :
        if j+1<len(uniq) and uniq[j+1]==p+1 :
            keys_p=map_k[p]
            keys_p1=map_k[p+1]
            if condition is not None :
                keys_p = screen_complementary_dic(coverage_res, condition=condition ,consider_only_keys=keys_p, debug=False)
                keys_p1 = screen_complementary_dic(coverage_res, condition=condition ,consider_only_keys=keys_p1, debug=False)
            if merge_only_best is not None : # use keys from condition if given, otherwise use all
                if sample_all_initialisers :
                    init_sort_key = sort_complementary_dic_each_initialiser(coverage_res,consider_only_keys=keys_p, sort_keys=sort_keys,reverse=reverse)
                    #keys_p=[]
                    #for ik in  init_sort_key : keys_p += init_sort_key[ik]
                    init_sort_key1 = sort_complementary_dic_each_initialiser(coverage_res,consider_only_keys=keys_p1, sort_keys=sort_keys,reverse=reverse)
                    #keys_p1=[]
                    #for ik in  init_sort_key1 : keys_p1 += init_sort_key1[ik]
                else :
                    keys_p = sort_complementary_dic(coverage_res,consider_only_keys=keys_p, sort_keys=sort_keys,reverse=reverse)
                    keys_p1 = sort_complementary_dic(coverage_res,consider_only_keys=keys_p1, sort_keys=sort_keys,reverse=reverse)
            if not sample_all_initialisers :
                init_sort_key,init_sort_key1={},{}
                init_sort_key['a'],init_sort_key1['a']=keys_p,keys_p1
            pair_class=(p,p+1)
            for keys_p in list(init_sort_key.values()) :
                nbest=0
                for k in keys_p :
                    for keys_p1 in list(init_sort_key1.values()) :
                        nbest1=0
                        for k1 in keys_p1 :
                            new_entry,twoloop_composition_lines = _merge_two_loop_peptides(coverage_res, k, k1 ,pair_reference=str(pair_class) ,compostion_index_dict=compostion_index_dict,merge_only_lengths_first_second_tup=merge_only_lengths_first_second_tup,add_pretty_line_csv=add_pretty_line_csv)
                            if new_entry is None :
                                continue
                            twol[n] =new_entry
                            if composition_file is not None and twoloop_composition_lines is not None :
                                outfile_comp.write(' index_of_first_seq_res: %d position %d two_loop_key=%s actual_peps: %s %s\n' % ( twol[n][twol.hd['start residue index']], twol[n][twol.hd['start residue position']], str(n),twol[n][twol.hd['first_peptide']],twol[n][twol.hd['second_peptide']] ))
                                for line in  twoloop_composition_lines :
                                    outfile_comp.write(line)
                                outfile_comp.write('\n\n')
                            n+=1
                            nbest1+=1
                            if merge_only_best is not None and nbest1>=merge_only_best : 
                                break
                    nbest+=1
                    if merge_only_best is not None and nbest>=merge_only_best : 
                            break
    print("\nIdentified %d two-loop candidates" % len(twol))
    if composition_file is not None :
        print("Two loop outfile_composition printed in %s" % (outfile_composition)) 
        outfile_comp.close()
    if outfile_coverage is not None : 
        if add_pretty_line_csv :
            outfile_coverage=outfile_coverage.replace('.txt','.csv').replace('.tsv','.csv').replace('.dat','.csv')
        print("Two loop outfile_coverage printed in %s" % (outfile_coverage)) 
        twol.Print(outfile_coverage)
    return twol




class Interactome() :
    def __init__(self,protein_sequence=None,coverage_tuple=None, proteome_blastDB='/Users/francesco/Documents/proteomes/human_reference_proteome_filtered.fasta'):
        self.protein_sequence=protein_sequence # the protein of which we want to calculate possible beta-like interactions
        self.coverage_tuple = coverage_tuple # the tuple with the four different kinds of complementary peptides we obtain. These are p0,p1,a0,a1 corresponding to parallel without HB in first residue, parallel with, antiparlalle without and with respecively.
        self.proteome_blastDB=proteome_blastDB
        self.interesting_mistmatchN=[0,1]
        self.save_overall_up_to_mismatch=1
    def blast_coverage(self,coverage_tuple=None,**kwargs):
        if coverage_tuple is None: coverage_tuple=self.coverage_tuple
        if coverage_tuple is None : raise Exception("Coverage_tuple is None\n")
        self.overall, self.positions, self.number_of_seq_with_one_hit, self.tot = blast_coverage(coverage_tuple, blastDB=self.proteome_blastDB, interesting_mistmatchN=self.interesting_mistmatchN,save_overall_up_to_mismatch=self.save_overall_up_to_mismatch ,**kwargs)
        return





def find_ideal_pairs(Complement_class, antiparallel_pairs=True, target_length=9):
    '''
Problem - we need an ideal pair of longer length - need to rethink:
could start by joining fragments of maximum lenght to reach target length by requiring that the overlap region is as long as possible (as long as doesn't come from 
existing longer fragments).
    will start at maximum length pairs and will analyse subfragments to see if they belong to beta-pairs other than those orginating the longer fragments 
    '''
    C=Complement_class
    
    if antiparallel_pairs : frag_list=C.antiparallel_fragments
    else : frag_list=C.parallel_fragments
    # frag_list is list of lenght 8 where index 0 is fragment of length 1 and so on (length 8 in 2017 db has 43'347 entries)
    # thus length 3 is index 2 corresponding to -6
    n=-1
    while n>=-6 :
        for frag_seq in frag_list[n] : # frag_seq is amino acid sequence N to C terminus
            # we 
            cand_score=-1
            interesting_matches=[]
            for comp in frag_list[n][frag_seq] : # these comp belong to a list that should already be sorted according to count score.
                if comp.count_score>=cand_score :
                    interesting_matches+=[ comp ]
                    cand_score = comp.count_score
                elif cand_score > comp.count_score : 
                    break
            # now look for subfragments within interesting matches
            for comp in  interesting_matches :
                # should slide in components of length n-1 () and see if you have subfragments in addition to the originating one
                a=1
        n-=1
    return
    


# DIVIDE BY NUMBER OF COMPOSING FRAGMENTS?? Or maybe better by length of fragment so that fragments of different length may become more comparable? 
def score_candidate_from_len_dic(len_dic_score_pos_comp, relevant_ids, promiscuity_smooth=0.01):
    '''
    score the complementary_peptide according to the information in score_pos_comp,
      which should be already filtered for hb and parallel/antiparallel (basically should be returned by get_score_pos_comp_with_clustered_hb)
    using only the entries at position in relevant_indices, as returned by print_score_pos_comp
    '''
    count_contr=0
    prom_contr=0
    for length,ind in relevant_ids :
        c=len_dic_score_pos_comp[length][ind] # complementary fragment (composing the complementary peptide), comp is at index 2 in score_pos_comp
        count_contr+= (len( c[2].complementary_sequence)**2 ) * c[2].count_score 
        prom_contr-= promiscuity_smooth*len(c[2].complementary_sequence)*c[2].promiscuity_score
    return count_contr+prom_contr,count_contr,prom_contr


def get_solubility_bin(solubility_score):
    solbin= int(2*solubility_score+4)
    if solbin<=0 : return 0
    return solbin

def update_multi_paths_candidates(key,i,first_index,  relevant_ids_composition,len_dic_score_pos_comp,     multi_paths_candidates):
    add_one=False
    #if key[0]=='YMRIVVVVEYRL' : 
        #print "***=== DEBUG YMRIVVVVEYRL  i=",i,'first_index=',first_index,'key=',key,'\n relevant_ids_composition[first_index]=',relevant_ids_composition[first_index]
        #for tup_id in relevant_ids_composition[first_index] : print '   Xtup_id',tup_id,len_dic_score_pos_comp[tup_id[0]][tup_id[1]]
    for tup_id in relevant_ids_composition[i] : # i is newly found fragment, first_index first ever found
        if tup_id not in relevant_ids_composition[first_index] : 
            # check that it is not fully contained in another longer fragment (e.g. an initialiser)
            C=len_dic_score_pos_comp[tup_id[0]][tup_id[1]] # aka [length][ind] ; this is the new fragment we may add
            remove =[] #maybe some of the existing fragments are contained in C rather than the other way round
            add=True
            for length,ind in relevant_ids_composition[first_index] : # loop on existing fragments
                C1=len_dic_score_pos_comp[length][ind] # we are not going to check HB pattern as this function is only used among fragments filtered by HB pattern
                if C[2].complementary_sequence in C1[2].complementary_sequence and C[1]>=C1[1] and C[1]+tup_id[0]<=length+C1[1] : add=False
                elif C1[2].complementary_sequence in C[2].complementary_sequence and C1[1]>=C[1] and C[1]+tup_id[0]>=length+C1[1] :remove+=[(length,ind)]
            if add :
                #if key[0]=='YMRIVVVVEYRL' : print "***=== adding one tup_id=",tup_id, 'i=',i,len_dic_score_pos_comp[tup_id[0]][tup_id[1]]
                add_one=True # if we add at least one fragment we add this to multi_paths_candidates as an extra path
                relevant_ids_composition[first_index]+=[tup_id] # we do this so that the final scores and the composition we print take into account all possible cascade paths that lead to this candidate.
            if remove!=[] :
                #if key[0]=='YMRIVVVVEYRL' : print "***=== remove=",remove, 'i=',i
                for tup_id1 in remove : relevant_ids_composition[first_index].remove(tup_id1)
    if add_one :
        if key not in multi_paths_candidates : # handle identical peptides that have been built through different paths
            multi_paths_candidates[key]=[ first_index , 2 ] # save the first occurrence as the reference one, 2 is the count (complementary_dic[key][0] is the index i)
        else : multi_paths_candidates[key][1]+=1 # update the count
    #if key[0]=='YMRIVVVVEYRL' : print "***=== DEBUG YMRIVVVVEYRL  end multi_paths_candidates[key]=", multi_paths_candidates[key]
    return multi_paths_candidates, relevant_ids_composition

def score_candidate(score_pos_comp, relevant_indices, add_sub_frag_contribution=True,sub_frag_smooth=0.25,promiscuity_smooth=0.01) :
    '''
    score the complementary_peptide according to the information in score_pos_comp,
      which should be already filtered for hb and parallel/antiparallel (basically should be returned by get_score_pos_comp_with_clustered_hb)
    using only the entries at position in relevant_indices, as returned by print_score_pos_comp
    '''
    # first separate the fragment that were actually employed in the cascade from the one that are parts of other fragments (e.g. fragments of length 3 that are sub-parts of fragments of greater length used in the cascade).
    sub_frags=[]
    main_frags=[]
    done=[]
    
    for i in relevant_indices :
        c=score_pos_comp[i]
        covers=list(range( c[1],c[-1]))
        if any([ all( [ p in a_done for p in covers ]) for a_done in done ]) :
            sub_frags+=[i]
        else :
            main_frags+=[i]
        done+=[ covers ]
    
    count_contr=0
    prom_contr=0
    for i in main_frags :
        c=score_pos_comp[i]
        count_contr+= (len(c[2].complementary_sequence)**2)*c[2].count_score 
        prom_contr-= promiscuity_smooth*len(c[2].complementary_sequence)*c[2].promiscuity_score
    if add_sub_frag_contribution :
        for i in sub_frags :
            c=score_pos_comp[i]
            count_contr+= sub_frag_smooth*((len(c[2].complementary_sequence)**2)*c[2].count_score )
            prom_contr -= sub_frag_smooth*(promiscuity_smooth*len(c[2].complementary_sequence)*c[2].promiscuity_score)
    return count_contr+prom_contr,count_contr,prom_contr

"""
FUNCTIONS TO GENERATE THE COVERAGE ON A DATABASE or on a single protein
"""
def parse_disordered_regions_from_DisProt(disprot_fasta_file,outfilename=None,clean_outfilename='DisProtDB_clean.fasta',minimum_length=8 , clean_disprot_naming=True):
    '''
    since the disprot database contain also structured regions it parses it and saves a fasta file with sequences corresponding to the disordered parts of 
    the sequences in the DisProt DB.
    if clean_outfilename is given it saves a fasta file like the disprot DB but with the # disordered regions sorted and merged.
_=beta_strands.parse_disordered_regions_from_DisProt("DisProtDB.fasta",clean_outfilename='DisProtDB_clean.fasta', outfilename='DisProtDisorderedBits.fasta')
    '''
    from  Bio.SeqRecord import SeqRecord
    new_records=[]
    if clean_outfilename is not None : clean_records=[]
    records=mybio.get_SeqRecords(disprot_fasta_file)
    for i,rec in enumerate(records) :  
        fn=rec.id
        if clean_disprot_naming :
            try :
                fn=fn.split('|')[3] # should yield the Uniprot AC
            except Exception :
                fn=rec.id
                pass  
        ranges=get_annotated_ranges(rec.description, symbol='#', sid=str((rec.id,i)) )
        struranges=get_annotated_ranges(rec.description, symbol='&', sid=str((rec.id,i)), print_warning=False)
        # check if annotated structured and disordered regions overlap (this DB is a mess)
        tmpranges=ranges[:]+struranges[:]
        tmpranges.sort()
        skip_record=False
        for j,r in enumerate(tmpranges[:-1]) :
            if tmpranges[j+1][0]+3 <= r[1] : # the first 3 residues of a structurally annotated region could be disordered
                skip_record=True
                break
        if skip_record :
            print('WARNING SKIPPING record %s %d as annotated structured and disordered regions overlap (%s   %s) !!' % (rec.id,i,str(r),str(tmpranges[j+1])))
            continue
        l_ids=[]  
        for start,end in ranges :
            l_ids+=[  '#%d-%d' % (start,end)  ]
            if end+1-start < minimum_length :# the end is included (we decided)
                continue
            new_id=fn+'#%d-%d' % (start,end)
            new_seq=rec.seq[start-1:end]
            new_records+=[ SeqRecord(seq=new_seq,id=new_id, description='disordered residues from %d to %d' % (start,end))]
        if clean_outfilename is not None : 
            if struranges!=[] or ranges!=[] :
                s_ids=[]
                for start,end in struranges :
                    s_ids+=[  '&%d-%d' % (start,end)  ]
                clean_records+=[ SeqRecord(seq=rec.seq , id=rec.id, description=' '.join(l_ids+s_ids) ) ]
    if outfilename is not None :
        mybio.PrintSequenceRecords(new_records, outfilename)
    if clean_outfilename is not None :
        mybio.PrintSequenceRecords(clean_records, clean_outfilename)
    return new_records

def get_annotated_ranges(description,symbol='#', sid='sequence',print_warning=True):
    '''
    auxiliary function used to clean the DisProt DB
    '''
    ranges=None
    if symbol in description :
        des=description.split()
        ranges=[]
        for el in des :
            if symbol in el :
                ranges+= [ list(map(int,el.replace(symbol,'').split('-'))) ]
    else :
        if print_warning : print('NO %s Symbol in description of %s' % (symbol,str(sid)))
        return []
    # check for overlapping ranges and in case merge them (some entries like P14061 have overlaps)
    return misc.merge_ranges(ranges)


def fasta_to_coverage_results_folder(fasta_file, comp_len=8,every=1,outfolder=None,copy=None, clean_disprot_naming=True,logfile=sys.stdout,limit_new_candidates=10,fork_ncpu=None,score_candidates=True):
    '''
    creates a folder with one coverage results file for each sequence in the fasta file

# DisProt
beta_strands.fasta_to_coverage_results_folder('DisProtDisorderedBits.fasta',outfolder='DisProtDB_COVERAGE/',copy=C,comp_len=8,every=1,limit_new_candidates=10,score_candidates=True,fork_ncpu=6)
# Delta2D
beta_strands.fasta_to_coverage_results_folder('delta2d_proteins_disordered_bits20.fasta', comp_len=8,every=1,outfolder=None,copy=C, clean_disprot_naming=True,logfile=sys.stdout,limit_new_candidates=10,fork_ncpu=None,score_candidates=True)
    '''
    if outfolder is None :
        _,dbid,_  = misc.get_file_path_and_extension(fasta_file)     
        outfolder='COMPLEMENTARYDB_'+str(comp_len)+'_'+dbid+'/'
    if outfolder[-1]!='/': outfolder+='/'
    os.system('mkdir '+outfolder)
    records=mybio.get_SeqRecords(fasta_file)
    c=Complement(copy=copy)
    c.score_candidates=score_candidates
    c.logfile=None
    c.warning_file=None
    c.out_file_composition=None
    c.out_filename=None
    c.limit_new_candidates=limit_new_candidates
    fname_list=[]
    N=len(records)
    failed=0
    skipped_non_standard_aa=0
    if  logfile is not None : logfile.write('%d proteins to analyze\n\n' % (N))
    good_records=[]
    if fork_ncpu is not None and fork_ncpu>1 :
        misc.reset_forking()
    for j,rec in enumerate(records) : # remove bad ones, so that if we fork it should not be a problem
        if len(rec.seq)<comp_len :
            skipped_non_standard_aa+=1
            print('WARN in fasta_to_coverage_results_folder() Sequence %s %d in file %s too short!! (%d). Skipping!!' % (rec.id,j,fasta_file,len(rec.seq)))
            continue
        
        seq=str(rec.seq)
        skip=False # check that sequence is standard
        for aa in seq :
            if aa not in amino_list1 :
                print('WARN in fasta_to_coverage_results_folder() Sequence %s %d in file %s contains non standard amino acid (%s). Skipping!!' % (rec.id,j,fasta_file,aa))
                skip=True
                break
        if skip:
            skipped_non_standard_aa+=1 
            continue
        good_records+=[rec]
    
    for j,rec in enumerate(good_records) :
        seq=str(rec.seq)
        try :
            fn=rec.id
            if clean_disprot_naming :
                try :
                    fn=fn.split('|')[3] # should yield the Uniprot AC
                except Exception :
                    fn=rec.id
                    pass
            fname=fn
            i=0
            while fname in fname_list :
                fname=fn+'_%02d' % (i)
                i+=1
                
            # do the coverage for this protein and save the results
            if fork_ncpu is None or fork_ncpu<=1 :
                _ = protein_coverage(seq,comp_len=comp_len,every=every,copy=c,outpickle=None ,logfile=None,warn_file=None,score_candidates=score_candidates,debug=False , save_coverage_results=outfolder+fname+'.tsv' , just_save_coverage_results=True)
            else :
                misc.fork_function(None, fork_ncpu, len(good_records), j,   protein_coverage,       seq,comp_len=comp_len,every=every,copy=c,outpickle=None ,logfile=None,warn_file=None,score_candidates=score_candidates,debug=False , save_coverage_results=outfolder+fname+'.tsv' , just_save_coverage_results=True)
            
            if logfile is not None and j%10==0 :
                logfile.write('%d->%4.1lf %%\t' % (j,100.*j/N))
                logfile.flush() 
        except Exception :
            sys.stderr.write('\n\n***ERROR*** at sequence %d %s\n' % (j,fn))
            print_exc(file=sys.stderr)
            sys.stderr.flush()
            failed+=1
            continue
    if  logfile is not None : logfile.write('\n\n%d proteins processed, of these %d failed and %d were skipped because containing non standard residues\n\n' % (N,failed,skipped_non_standard_aa))
    return

def coverage_statistics_from_result_folder(result_folder,only_process_files_containing=['.tsv'],avoid_regions_shorter_than=None,min_initializer_len=4,min_complementarity_score=30.,min_Solubility_score=-1.5,overall_results_filename=None,warn_for_empty=True,save_all_results=True, return_all_complementary_scores=True):
    '''
    from a folder generated by fasta_to_coverage_results_folder it generates a summary with the coverage using the given thresholds (None ==> no condition on that variable)
result_folder='DisProtDB_COVERAGE_new'
overall=beta_strands.coverage_statistics_from_result_folder(result_folder,only_process_files_containing=['.tsv'],min_initializer_len=4,min_complementarity_score=30.,min_Solubility_score=-1.5)
    
    
# DELTA 2D results
overall,scores,all_comp_scores=beta_strands.coverage_statistics_from_result_folder('COMPLEMENTARYDB_8_delta2d_proteins_disordered_bits20/',avoid_regions_shorter_than=20,min_complementarity_score=0,min_Solubility_score=None,min_initializer_len=None)

    '''
    flist=os.listdir(result_folder)
    if result_folder[-1]!='/' : result_folder+='/'
    
    # do iterative for multiple condtions if asked to
    if type(min_complementarity_score) is list :
        run_res=csv_dict.Data()
        rhd=[]
        for k in ['positions_covered', 'antiparallel_positions_covered', 'parallel_positions_covered', 'a0_positions_covered', 'a1_positions_covered', 'p0_positions_covered', 'p1_positions_covered'] :
            rhd+=[ k, 'NoEmpties_'+k,'N_'+k]
        run_res.hd=csv_dict.list_to_dictionary(rhd+['tot_processed','tot_in_empties'])
        for minscore in min_complementarity_score :
            _,scores= coverage_statistics_from_result_folder(result_folder,overall_results_filename='Overall_scoreth_%d.tsv' % (minscore), only_process_files_containing=only_process_files_containing,avoid_regions_shorter_than=avoid_regions_shorter_than,warn_for_empty=False,min_initializer_len=min_initializer_len,min_complementarity_score=minscore,min_Solubility_score=min_Solubility_score,save_all_results=False, return_all_complementary_scores=False)
            run_res[minscore]=scores[:]
        return run_res
    
    flist=[ f for f in flist if f[0]!='.' and misc.loose_compare(f, only_process_files_containing)[0] ]
    print('%d files to anlayze in folder %s' % (len(flist),result_folder))
    all_results=csv_dict.Data()
    overall=csv_dict.Data()
    overall.hd=csv_dict.list_to_dictionary([ 'len_protein', 'positions_covered', 'antiparallel_positions_covered','a0_positions_covered', 'a1_positions_covered', 'parallel_positions_covered', 'p0_positions_covered', 'p1_positions_covered' \
                                            , 'relative_coverage', 'antiparallel_relative_coverage', 'a0_relative_coverage', 'a1_relative_coverage', 'parallel_relative_coverage', 'p0_relative_coverage', 'p1_relative_coverage'] )
    iii=0
    empties=0
    tot_skipped_residues=0
    total_residues_processed=0
    protein_cov={}
    peptide_considered=0
    peptide_discarded=0
    if return_all_complementary_scores : all_comp_scores=[]
    for f in flist :
        if 'overall' in f.lower() :
            print('SKIPPING ',f) 
            continue
        if 'Global_merged_results' in f :
            print('SKIPPING ',f) 
            continue
        if 'raw_list_protein_dic_coverage' in f :
            print('SKIPPING ',f) 
            continue
        _,pname,_=misc.get_file_path_and_extension(f)
        if pname in protein_cov :
            print('**ERROR** %s already in protein_cov, overwriting!!' % (pname))
        protein_cov[pname]={}
        protein_cov[pname]['tot']=[]
        d=csv_dict.Data()
        d.load(result_folder+f,verbose=False)
        if avoid_regions_shorter_than is not None :
            if '#' in pname :
                st,end=list(map(int, pname[pname.find('#')+1:].split('-') ))
                #print pname,st,end,end-st+1
                if end-st+1 < avoid_regions_shorter_than :
                    #print '   %s short' % (pname), 
                    continue
            else :
                print('**ERROR** cant assess lenght of entry %s to apply avoid_regions_shorter_than' % (f))
        if len(d)==0 :
            empties+=1
            if '#' in pname :
                st,end=list(map(int, pname[pname.find('#')+1:].split('-') ))
                tot_skipped_residues += (end-st+1)
            if warn_for_empty : print('**WARNING** empty file %25s , skipping completely!  %s %d' % (f,pname,(end-st+1)))
            continue
        protein_cov[pname]['len_protein']=list(d.values())[0][d.hd['len_protein']]
        total_residues_processed+=protein_cov[pname]['len_protein']
        for k in d :
            if return_all_complementary_scores : all_comp_scores+=[ (d[k][d.hd['Complementarity_score']], d[k][d.hd['len_initializer']],d[k][d.hd['Solubility_score']] ) ]
            if (min_initializer_len is None or d[k][d.hd['len_initializer']] >= min_initializer_len) \
                and (min_complementarity_score is None or d[k][d.hd['Complementarity_score']] >= min_complementarity_score) \
                and (min_Solubility_score is None or d[k][d.hd['Solubility_score']] >= min_Solubility_score) :
                
                peptide_considered+=1
                iii+=1
                if save_all_results :
                    if all_results.hd=={} :
                        all_results.hd=d.hd.copy()
                        all_results._update_hd('protein_name', 0)
                    all_results[iii]=[pname]+d[k]
                key=d[k][d.hd['parallel/antiparallel']]+str(d[k][d.hd['FirstResHB']])
                if key not in protein_cov[pname] : protein_cov[pname][key]=[]
                key2=d[k][d.hd['parallel/antiparallel']]
                if key2 not in protein_cov[pname] : protein_cov[pname][key2]=[]
                start_ind= d[k][d.hd['index_of_first_seq_res']]
                for i in range(start_ind , start_ind + len(d[k][d.hd['Complementary_sequence']])) :
                    if i not in protein_cov[pname][key] :
                        protein_cov[pname][key].append(i)
                    if i not in protein_cov[pname][key2] :
                        protein_cov[pname][key2].append(i)
                    if i not in protein_cov[pname]['tot'] :
                        protein_cov[pname]['tot'].append(i)
            else :
                peptide_discarded+=1
        for k in protein_cov[pname] :
            if k not in ['len_protein','tot','-','-0','-1','+','+0','+1'] :
                print('ERROR unrecognized key %s in protein_cov, not considering the values in it' % (k))
        count={}
        for k in ['tot','-','-0','-1','+','+0','+1'] : count[k]=0
        for i in range(0,protein_cov[pname]['len_protein']) :
            for k in ['tot','-','-0','-1','+','+0','+1'] :
                if k in protein_cov[pname] and i in protein_cov[pname][k] :
                    count[k]+=1
        overall[pname]=[ protein_cov[pname]['len_protein'] ]+ [ count[k] for k in ['tot','-','-0','-1','+','+0','+1'] ] + [ float(count[k])/float(protein_cov[pname]['len_protein']) for k in ['tot','-','-0','-1','+','+0','+1'] ]
    
    if overall_results_filename is None : overall_results_filename='Overall_coverage_results.tsv'
    overall.Print(result_folder+overall_results_filename)
    print(' overall results saved in ',result_folder+overall_results_filename)
    if save_all_results :
        all_results.Print(result_folder+'Global_merged_results.tsv')
        print(' merged results saved in ',result_folder+'Global_merged_results.tsv')
        out=open(result_folder+'raw_list_protein_dic_coverage.pkl','wb')
        pickle.dump(protein_cov ,out)
        out.close()
        print(' raw list protein coverage pickle saved in',result_folder+'raw_list_protein_dic_coverage.pkl')
    
    scores=[]
    print('peptide_considered',peptide_considered,'peptide_discarded',peptide_discarded,'considered:', 100.*peptide_considered/float(peptide_discarded+peptide_considered),'%% ; keep in mind this is heavily influenced by the parameter limit_new_candidates given to fasta_to_coverage_results_folder(), which selects only the best candidates a priori\n')
    print('%d EMPTY files have been skipped among %d ==> tot_skipped_residues= %d' % (empties,len(flist),tot_skipped_residues))
    print('total_residues_processed',total_residues_processed)
    for k in ['positions_covered', 'antiparallel_positions_covered', 'parallel_positions_covered', 'a0_positions_covered', 'a1_positions_covered', 'p0_positions_covered', 'p1_positions_covered'] :
        su= sum(overall.column(k))
        print('%s coverage: %d  %lf  and including skipped %lf' % (k,su, float(su)/total_residues_processed, float(su)/(total_residues_processed+tot_skipped_residues)))
        scores+= [ float(su)/(total_residues_processed+tot_skipped_residues), float(su)/total_residues_processed,su]
    scores+=[total_residues_processed,tot_skipped_residues]
    if return_all_complementary_scores : return overall,scores,all_comp_scores
    return overall,scores


def plot_minscore_screen(run_res, bootstrap=False, key_frac='positions_covered',key_count='N_positions_covered',key_tot=['tot_processed','tot_in_empties'], figsize=(14,10),label=None,ls='',marker='.',color=None,ms=10,set_to_percent=True,figure=None,vline=None,save=None):
    import plotter
    plotter.set_publish()
    if type(run_res) is list :
        fig=None
        if type(color) is list : cols=color 
        else :cols=plotter.iworkpalette*2
        if type(marker) is list : marks=marker
        else : marks=['o','D','s','^']*2
        for j,run in enumerate(run_res) :
            if j==0 :vl =vline
            else : vl=None
            fig,_=plot_minscore_screen(run,bootstrap=bootstrap,key_frac=key_frac,key_count=key_count,key_tot=key_tot,ls=ls, label=None,ms=ms,set_to_percent=set_to_percent, figsize=figsize,figure=fig, marker=marks[j] , color=cols[j],vline=vl ,save=None)
        if label is not None :
            plotter.add_custom_legend(label, cols, marker_types=marks, markersize=ms, figure=fig, frame=False, legend_location='upper right', legend_size=plotter.publication_small['legend_size'])
        if save  is not None and save != '':
            if '.' not in save: save += '.pdf'
            fig.savefig(save, dpi=300, bbox_inches="tight", bbox_extra_artists=[fig], transparent=True) # ,bbox_inches="tight" is such that everything is in the figure... even though margins are removed.
        return fig
    if bootstrap==False or bootstrap is None :
        profile=run_res.column(key_frac)
        if set_to_percent :
            profile=100.*numpy.array(profile)
        yerr=None
    else :
        ms=1
        if type(key_tot) is not list and type(key_tot) is not tuple : key_tot=[key_tot]
        if type(bootstrap) is int : cycles=bootstrap
        else : cycles=1000
        profile=[]
        yerr=[]
        low_distance_to_95CI=[]
        up_distance_to_95CI=[]
        if set_to_percent :m=100.
        else : m=1.
        for k in run_res :
            Num=run_res[k][run_res.hd[key_count]]
            tot=0
            for kt in key_tot : tot+=run_res[k][run_res.hd[kt]]
            outcomes=[1]*Num+[0]*(tot-Num)
            #un,bootstrap_mean,bootstrap_std=misc.bootstrap(outcomes, bootstrap_runs=bootstrap, return_only_count=True)
            count,un= misc.get_bootstrap_ensemble(outcomes, bootstrap_runs=bootstrap,return_only_count=True)
            bootstrap_mean=numpy.mean( count,axis=1)
            bootstrap_std=numpy.std(count,axis=1)
            ci25,ci975= numpy.percentile(count, [2.5,97.5], axis=1)
            
            ind=un.index(1)
            count_err=bootstrap_std[ind]
            count_ci25,count_ci975=ci25[ind],ci975[ind]
            if not misc.CompareFloats(Num,bootstrap_mean[ind] , sensibility=5) :
                print('*Warn bootstrap might not have converged for %s count is %d and corresponding bootstrap mean %lf' % (str(k),Num,bootstrap_mean[ind]))
            if not misc.CompareFloats(float(Num)/tot,run_res[k][run_res.hd[key_frac]] , sensibility=0.0001) :
                print('ERROR inconsistent count and fraction for %s count, tot is %d and corresponding fraction is %lf but calculated %lf' % (str(k),Num,tot,run_res[k][run_res.hd[key_frac]],float(Num)/tot))
            profile+=[ m*float(Num)/tot ]
            low_distance_to_95CI+=[ m*float(Num-count_ci25)/tot ]
            up_distance_to_95CI +=[ m*float(count_ci975-Num)/tot ] 
            yerr+=[ m* float(count_err)/tot ]
        run_res.add_column(up_distance_to_95CI,key_frac ,variable_name='DistTo95CI_up_'+key_frac)
        run_res.add_column(low_distance_to_95CI,key_frac ,variable_name='DistTo95CI_low_'+key_frac)
        run_res.add_column(yerr,key_frac ,variable_name='SDERR'+key_frac)

    xvals=list(map(float,list(run_res.keys())))
    if set_to_percent : 
        y_range=(0,100.)
        ym=10
        yM=25
    else : 
        y_range=(0,1)
        ym=0.1
        yM=0.25
    plotter.set_publish()
    if color is None : col=plotter.iworkpalette[0]
    else : col=color
    fig=plotter.profile(profile,xvals, yerr=yerr, label=label,x_minor_tick_every=10,x_major_tick_every=20,figure_size=figsize,frame=False,hgrid=False,ls=ls,y_range=y_range,y_major_tick_every=yM,y_minor_tick_every=ym,color=col,markerfacecolor=col,markeredgecolor='black',vline=vline,markersize=ms,marker=marker,figure=figure,save=save )
    return fig,run_res

def global_coverage_results_to_aa(fasta_file, global_results, threshold=10, add_skipped_from_fasta=True):
    '''
    as the analysis of the coverage folders skip empty files (that correspond to sequences for which no complementary exists) setting add_skipped_from_fasta
      recovers this sequences. However it works only if the fasta file is the same one that was used to generate the folder, and if all sequences where used to generate it.
    after a coverage folder has been analysed this can be used to extract amino acid info about
     sequences that have/don't have complementaries
    (see plot_aa_freqs for example usage)
    '''
    if type(global_results) is str :
        res=csv_dict.Data()
        res.load(global_results)
    else : 
        res=global_results
    prot_dic=mybio.get_SeqRecords(fasta_file, return_as_dictionary=True)
    print('Loaded %d proteins from %s' % (len(prot_dic),fasta_file))
    frag_dic={}
    tmp_frag_with_zero={}
    frag_with_th=None
    if threshold is not None :
        tmp_frag_th={}
        frag_with_th=[]
    frag_with_zero=[]
    frag_with_more1=[]
    for k in res :
        pname= res[k][res.hd['protein_name']] 
        if '.tsv' in pname : pname=pname.replace('.tsv','')
        
        if pname in prot_dic :
            if pname not in tmp_frag_with_zero :
                tmp_frag_with_zero[pname]= str(prot_dic[pname].seq)
                if threshold is not None :
                    tmp_frag_th[pname]={}
            st=res[k][res.hd['index_of_first_seq_res']]
            en=st+len(res[k][res.hd['Complementary_sequence']])
            protein_bit= str(prot_dic[pname].seq)[  st : en ]
            if protein_bit not in frag_dic: frag_dic[protein_bit]=0
            frag_dic[protein_bit]+=1
            if threshold is not None :
                if (st,en) not in tmp_frag_th[pname] : tmp_frag_th[pname][(st,en)]=0
                tmp_frag_th[pname][(st,en)]+=1
            tmp_frag_with_zero[pname]= tmp_frag_with_zero[pname][:st] + tmp_frag_with_zero[pname][st:en].lower() +tmp_frag_with_zero[pname][en:] # so that the upper case letters will be the ones without complementaries
        else :
            sys.stderr.write("**WARNING** %s not found in fasta file!!\n SKIPPED\n" % (pname))
            continue
    for pname in tmp_frag_with_zero :
        prev_upper=False
        for aa in tmp_frag_with_zero[pname] :
            if aa.isupper() :
                if not prev_upper : frag_with_zero+=[ aa ]
                else : frag_with_zero[-1]+=aa
                prev_upper=True
            else : 
                if prev_upper : frag_with_more1+=[ aa.upper() ]
                else :
                    if frag_with_more1==[] : frag_with_more1=['']
                    frag_with_more1[-1] += aa.upper() 
                prev_upper=False
    if add_skipped_from_fasta :
        for pname in prot_dic :
            if pname not in tmp_frag_with_zero :
                frag_with_zero+=[ str(prot_dic[pname].seq) ]
    if threshold is not None :
        for pname in tmp_frag_th :
            pseq=str(prot_dic[pname].seq)
            for t in tmp_frag_th[pname] :
                if tmp_frag_th[pname][t]>= threshold :
                    pseq=pseq[:t[0]]+pseq[t[0]:t[1]].lower()+pseq[t[1]:]
            prev_upper=False
            for aa in pseq :
                if aa.isupper() :
                    #if not prev_upper : frag_with_zero+=[ aa ]
                    #else : frag_with_zero[-1]+=aa
                    prev_upper=True
                else : 
                    if prev_upper : frag_with_th+=[ aa.upper() ]
                    else :
                        if frag_with_th==[] : frag_with_th=['']
                        frag_with_th[-1] += aa.upper() 
                    prev_upper=False
    
    return frag_dic, frag_with_zero,frag_with_more1, frag_with_th

def plot_aa_freqs(frag_dic, frag_with_zero,fasta_reference=None,threshold=10, bootstrap=False, plot_difference_from_fasta=True):
    '''
    plots the results of the function global_coverage_results_to_aa()

fasta='DisProtDisorderedBits.fasta'
frag_dic, frag_with_zero,frag_with_more1, frag_with_th=beta_strands.global_coverage_results_to_aa(fasta,'Global_merged_coverage_DisProt_noTH.tsv')
save='coverageDisProt_aa_composition.png'
# OR from delta2d 
fasta='delta2d_proteins_disordered_bits20.fasta'
frag_dic, frag_with_zero,frag_with_more1, frag_with_th=beta_strands.global_coverage_results_to_aa(fasta,'Global_merged_coverage_delta2D_20_noTH.tsv')
save='coverageD2D_aa_composition.png'
#  OR from the merged DBs
fasta='all_disordered_bits_merged.fasta'
frag_dic, frag_with_zero,frag_with_more1, frag_with_th=beta_strands.global_coverage_results_to_aa(fasta,'Global_coverage_results_bothDBs_noTH.tsv')
save='coverageBothDBs_aa_composition.png'

entries, yerr,aas=beta_strands.plot_aa_freqs(frag_with_more1, frag_with_zero,fasta,threshold=frag_with_th, bootstrap=1000,plot_difference_from_fasta=True)
plotter.set_publish(small_figure=False)
plotter.multibar(entries, yerr=yerr,xlabels=aas ,figure_size=(6,8),color=plotter.iworkpalette,hgrid=False,vgrid=False,y_minor_tick_every=0.05,y_major_tick_every=0.1,frame=False,hline=0, save=save)
    '''
    if fasta_reference is not None :
        if type(fasta_reference) is str :
            fasta_reference=mybio.get_SeqRecords(fasta_reference)
            ref_count,ref_tot=mybio.composition(fasta_reference,bootstrap=bootstrap, only_elements_in_list=amino_list1)
    noc_count,noc_tot=mybio.composition(frag_with_zero, bootstrap=bootstrap, only_elements_in_list=amino_list1)
    if type(frag_dic) is dict : comp_count,comp_tot=mybio.composition(list(frag_dic.keys()), bootstrap=bootstrap, only_elements_in_list=amino_list1)
    else : comp_count,comp_tot=mybio.composition(frag_dic, bootstrap=bootstrap, only_elements_in_list=amino_list1)
    
    if type(threshold) is list  :
        th_count,th_tot= mybio.composition( threshold, bootstrap=bootstrap, only_elements_in_list=amino_list1)
    elif threshold is not None and type(frag_dic) is dict : 
        ths=[k for k in frag_dic if frag_dic[k]>=threshold]
        th_count,th_tot= mybio.composition( ths, bootstrap=bootstrap, only_elements_in_list=amino_list1)
    else : threshold=None
    aas=sorted(amino_list1)
    entries=[]
    yerr=None
    if bootstrap :yerr=[]
    for aa in aas :
        entries+=[ [] ]
        if bootstrap :yerr+=[ [] ]
        if aa in noc_count : 
            noc_count[aa]/=float(noc_tot)
            if bootstrap :
                entries[-1]+=[ noc_count[aa][0] ]
                yerr[-1]+=[ noc_count[aa][1] ]
            else :
                entries[-1]+=[ noc_count[aa] ]
        if threshold is not None and aa in th_count : 
            th_count[aa]/=float(th_tot)
            if bootstrap :
                entries[-1]+=[ th_count[aa][0] ]
                yerr[-1]+=[ th_count[aa][1] ]
            else :
                entries[-1]+=[ th_count[aa] ]
        if aa in comp_count : 
            comp_count[aa]/=float(comp_tot)
            if bootstrap :
                entries[-1]+=[ comp_count[aa][0] ]
                yerr[-1]+=[ comp_count[aa][1] ]
            else :
                entries[-1]+=[ comp_count[aa] ]
        if fasta_reference is not None and aa in ref_count : 
            ref_count[aa]/=float(ref_tot)
            if bootstrap :
                entries[-1]+=[ ref_count[aa][0] ]
                yerr[-1]+=[ ref_count[aa][1] ]
            else :
                entries[-1]+=[ ref_count[aa] ]
    entries=numpy.array(entries).T
    if bootstrap : yerr=numpy.array(yerr).T
    # plotter.multibar(entries, yerr=yerr,xlabels=aas )
    if fasta_reference is not None and plot_difference_from_fasta :
        entries=(entries-entries[3])[:3:]
        yerr=numpy.sqrt((yerr**2+yerr[3]**2)[:3:])
    return entries, yerr,aas
    
    
def plot_complementarity_score_distribution(all_comp_scores, only_for_initializers_longer_than=3,x_range=(0,200), vline=30,figsize=(14,10),save='Complementarity_score_distribution.png' ):
    '''
    after a coverage folder has been analysed this can be used to plot the distribution of the complementarity score or other scores
    it needs the list all_comp_scores that is returned by coverage_statistics_from_result_folder() when return_all_complementary_scores=True
    '''
    import plotter
    res=[]
    if only_for_initializers_longer_than is None : only_for_initializers_longer_than=0
    for s in all_comp_scores :
        if s[1]> only_for_initializers_longer_than :
            res+=[s[0]]
    f,t,nf=numpy.percentile(res,[5,30,95])
    print('percentiles 5th, 30th, 95th :',f,t,nf)
    f=plotter.histogram(res,x_range=x_range,label=None,frame=False,hgrid=False,x_major_tick_every=40,x_minor_tick_every=20,color=plotter.iworkpalette[0],flag_labels=None,figure_size=figsize,show=False)
    if vline is not None :
        ax=f.gca()
        m,M=ax.get_ylim()
        ax.vlines(vline,m,M,color='red',linestyle='--',lw=3)
    plotter.plt.draw()
    plotter.plt.show(block=False)
    if save is not None :
        if '.png' in save:
            f.savefig(save, dpi=300, bbox_inches="tight", bbox_extra_artists=[ax, f],transparent=True) # ,bbox_inches="tight" is such that everything is in the figure... even though margins are removed.
        else:
            f.savefig(save, dpi=300, bbox_inches="tight", bbox_extra_artists=[ax, f]) # ,bbox_inches="tight" is such that everything is in the figure... even though margins are removed.
    return
    
    
    
# tot_coverage_stat,summary=beta_strands.build_coverage_database('DisprotDB.fasta',copy=c)
def build_coverage_database(fasta_file, comp_len=8,every=1,outfolder=None,copy=None, clean_disprot_naming=True,logfile=sys.stdout): 
    '''
    it generate the coverage database for a fasta file, basically it scans along each sequences every every position and
      generates all complementary peptide of length comp_len.
     copy is used to pass a Complement() class that has the BSn databases (fragments) already loaded..
        clean_disprot_naming is specific to the DisProt database and clean the protein names... (convert to Uniprot AC I think).
     outfolders saves all the pickles in a folder (if None the folder is authomatically created depending on the name of the fasta file)
    it also print a coverage_summary_ .tsv file in the outfolder
    Example :
    tot_coverage_stat,summary=beta_strands.build_coverage_database('DisprotDB.fasta',copy=c)
    '''
    if outfolder is None :
        _,dbid,_  = misc.get_file_path_and_extension(fasta_file)     
        outfolder='COMPLEMENTARYDB_'+str(comp_len)+'_'+dbid+'/'
    if outfolder[-1]!='/': outfolder+='/'
    os.system('mkdir '+outfolder)
    records=mybio.get_SeqRecords(fasta_file)
    c=Complement(copy=copy)
    c.score_candidates=False
    c.logfile=None
    c.warning_file=None
    c.limit_new_candidates=10
    summary=csv_dict.Data()
    hd=['coverage(%)','p0_coverage(%)','p1_coverage(%)','a0_coverage(%)','a1_coverage(%)','#different_sites_with_at_least_1_comlement','#different_sites_with_at_least_1_p0comlement','#different_sites_with_at_least_1_p1comlement','#different_sites_with_at_least_1_a0comlement','#different_sites_with_at_least_1_a1comlement']
    summary.hd=csv_dict.list_to_dictionary(hd)
    fname_list=[]
    tot_coverage_stat=structs.Stats('Coverage(%)')
    N=len(records)
    if  logfile is not None : logfile.write('%d proteins to analyze\n\n' % (N))
    for j,rec in enumerate(records) :
        try :
            seq=str(rec.seq)
            skip=False # check that sequence is standard
            for aa in seq :
                if aa not in amino_list1 :
                    print('WARN in build_coverage_database() Sequence %s %d in file %s contains non standard amino acid (%s). Skipping!!' % (rec.id,j,fasta_file,aa))
                    skip=True
                    break
            if skip: continue
            fn=rec.id
            if clean_disprot_naming :
                try :
                    fn=fn.split('|')[3] # should yield the Uniprot AC
                except Exception :
                    fn=rec.id
                    pass
            fname=fn
            i=0
            while fname in fname_list :
                fname=fn+'_%02d' % (i)
                i+=1
            p0,p1,a0,a1=protein_coverage(seq,comp_len=comp_len,every=every,copy=c,outpickle=None ,logfile=None)
            out=open(outfolder+fname+'_seq_p0p1a0a1.pkl','wb')
            pickle.dump((seq,p0,p1,a0,a1),out)
            out.close()
            plen=len(seq)
            total,p0C,p1C,a0C,a1C,totalM,p0M,p1M,a0M,a1M=coverage_results_to_percentages(plen, (p0,p1,a0,a1), comp_len)
            tot_coverage=100.*(plen-len(totalM))/plen
            p0_coverage =100.*(plen-len(p0M))/plen
            p1_coverage =100.*(plen-len(p1M))/plen
            a0_coverage =100.*(plen-len(a0M))/plen
            a1_coverage =100.*(plen-len(a1M))/plen
            summary[fname]=[tot_coverage,p0_coverage,p1_coverage,a0_coverage,a1_coverage,total,p0C,p1C,a0C,a1C]
            tot_coverage_stat.update(tot_coverage,fname)
            if logfile is not None and j%10==0 :
                logfile.write('%d->%4.1lf\t' % (j,100.*j/N))
                logfile.flush()
        except Exception :
            sys.stderr.write('\n\n***ERROR*** at sequence %d %s\n' % (j,fn))
            print_exc(file=sys.stderr)
            sys.stderr.flush()
            continue
    if  logfile is not None :logfile.write('\n\n')
    tot_coverage_stat.finalize()
    print('tot_coverage_stat.average',tot_coverage_stat.average)
    
    sname='coverage_summary_'+dbid+'.tsv'
    summary.Print(outfolder+sname)
    return tot_coverage_stat,summary

def coverage_results_to_percentages(protein_len, tuple_coverage_dict,comp_len): 
    if type(tuple_coverage_dict) is str :
        if 'seq' in tuple_coverage_dict :
            seq,p0,p1,a0,a1=pickle.load(open(tuple_coverage_dict,'rb'))
            protein_len=len(seq)
        else :
            p0,p1,a0,a1=pickle.load(open(tuple_coverage_dict,'rb'))
    else :
        p0,p1,a0,a1=tuple_coverage_dict
    p0C,p1C,a0C,a1C=0,0,0,0 # this will count the number of position where a complementary fragment begins!
    p0M,p1M,a0M,a1M=list(range(0,protein_len)),list(range(0,protein_len)),list(range(0,protein_len)),list(range(0,protein_len)) # this will count the number of positions that don't have any complementary fragment!
    total=0
    totalM=list(range(0,protein_len))
    for i in range(protein_len) :
        one_found=False
        one_missing=False
        if i in p0 and p0[i]!=[] :
            p0C+=1
            for j in range(i,i+comp_len) :
                if j in p0M :
                    p0M.remove(j)
                if j in totalM :
                    totalM.remove(j)
            one_found=True
        if i in p1 and p1[i]!=[] :
            p1C+=1
            for j in range(i,i+comp_len) :
                if j in p1M :
                    p1M.remove(j)
                if j in totalM :
                    totalM.remove(j)
            one_found=True
        if i in a0 and a0[i]!=[] :
            a0C+=1
            for j in range(i,i+comp_len) :
                if j in a0M :
                    a0M.remove(j)
                if j in totalM :
                    totalM.remove(j)
            one_found=True
        if i in a1 and a1[i]!=[] :
            a1C+=1
            for j in range(i,i+comp_len) :
                if j in a1M :
                    a1M.remove(j)
                if j in totalM :
                    totalM.remove(j)
            one_found=True
        if one_found :
            total+=1
    return total,p0C,p1C,a0C,a1C,totalM,p0M,p1M,a0M,a1M

# asyn='MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA'
def protein_coverage(protein, comp_len=8,every=1,copy=None,outpickle=None,debug=False ,save_coverage_results=None, score_candidates=True,logfile=sys.stdout,warn_file=sys.stderr, just_save_coverage_results=False):
    '''
used in paper, now screen_complementary_peptides() should be used instead.
    given a protein sequence it generate all possible complementary fragmetnts of length comp_len
    starting at every every residue from the beginning of the protein. 
    copy can be used to init the Complement() class from another class, so that the database is linked with a pointer (super fast) rather then reloaded (super slow)
    it returns parallel_hb_0, parallel_hb_1 , anti_hb_0 , anti_hb_1
     which are four dictionaries, with keys that are indices along the protein sequence of where the complementary peptide begins, 
      note that some positions may be present as key, but the corresponding value may be []. This means that althought there is no complementary 
      in e.g. anti_hb_0 for that position, there is at least one complementary in at least one of the other 3 dictionaries (anti_hb_1,parallel_hb_0, parallel_hb_1).
      positions not present in keys are positions for which NO complementary has been found in EVERY of the four dictionaries. 
     values are lists of all the possible complements found at that position.
    the four dictionaries contain either parallel or antiparallel complementary peptides and 
    hb_0 means that the first residue of that peptide is NOT involved in Hydrogen bond, 1 means that it is.
    if not just_save_coverage_results
     inner list contain 'Complementary_sequence','parallel/antiparallel','len_initializer','initializer'
     and if score_candidates also 'Complementarity_score', 'Solubility_score'
parallel_hb_0, parallel_hb_1 , anti_hb_0 , anti_hb_1= beta_strands.protein_coverage(Nasyn,8,copy=C)
_=beta_strands.visualize_coverage( Nasyn,(parallel_hb_0, parallel_hb_1 , anti_hb_0 , anti_hb_1), outfile=sys.stdout)
    '''
    #if copy is not None : copy.logfile=None
    c=Complement(copy=copy)
    c.score_candidates=score_candidates
    c.logfile=None
    c.warning_file=None #sys.stderr
    c.out_file_composition=None
    c.out_filename=None
    #print c.logfile
    parallel_hb_0={}
    parallel_hb_1={}
    anti_hb_0={}
    anti_hb_1={}
    if save_coverage_results==False : save_coverage_results=None
    if save_coverage_results is not None : 
        protein_coverage_res=csv_dict.Data()
        rhd=['index_of_first_seq_res' , 'parallel/antiparallel','FirstResHB','len_protein','Complementary_sequence','len_initializer','initializer','NumberOfPaths','Solubility_score']
        if score_candidates: rhd+=['Complementarity_score','Count_score','Promiscuity_score']
        protein_coverage_res.hd=csv_dict.list_to_dictionary(rhd)
        tot=0
    
    interval=list(range( 0 , len(protein)-comp_len , every))
    N=1.*len(interval)
    for j,i in  enumerate(interval):
        if logfile is not None and j%10==0 :
            logfile.write(' pos:%d -->%5.2lf %% |' % (i,100.*j/(N)))
            logfile.flush()
        c.complementary_dic={}
        c.make(protein[i:i+comp_len])
        results=c.complementary_dic
        hd=c.complementary_hd
        if len(results)==0 :
            if warn_file is not None : warn_file.write( '****** NO COMPLEMENTS found at index=%d for segment %s\n' % (i,protein[i:i+comp_len]))
            del results
            continue
        if debug : print(i,'len(results)',len(results),type(list(results.values())[0] ),list(results.values())[0],hd)
        if not just_save_coverage_results :
            parallel_hb_0[i]=[]
            parallel_hb_1[i]=[]
            anti_hb_0[i]=[]
            anti_hb_1[i]=[]
        for k in results :
            if save_coverage_results is not None :
                # hd should contain: 
                # complementary_hd=['CandidateIndex', 'Complementary_sequence','parallel/antiparallel','len_initializer','initializer','FirstResHB','NumberOfPaths','Complementarity_score','Solubility_score','Count_score','Promiscuity_score' ]
                tot+=1
                protein_coverage_res[tot]= [ i , results[k][hd['parallel/antiparallel']],results[k][hd['FirstResHB']],len(protein),results[k][hd['Complementary_sequence']],results[k][hd['len_initializer']],results[k][hd['initializer']],results[k][hd['NumberOfPaths']],results[k][hd['Solubility_score']] ]
                if score_candidates : protein_coverage_res[tot] += [ results[k][hd['Complementarity_score']],results[k][hd['Count_score']],results[k][hd['Promiscuity_score']] ]
            if not just_save_coverage_results :
                to_save=[ results[k][hd['Complementary_sequence']],results[k][hd['parallel/antiparallel']],results[k][hd['len_initializer']],results[k][hd['initializer']] ]  #[ results[k][1:hd['FirstResHB']] ] # 'parallel/antiparallel','len_initializer','initializer'
                if score_candidates : to_save += [ results[k][hd['Complementarity_score']],results[k][hd['Solubility_score']]  ]
                if results[k][hd['parallel/antiparallel']]=='+' :
                    if results[k][hd['FirstResHB']]==0 : parallel_hb_0[i]+= [ to_save ]
                    else : parallel_hb_1[i]+= [ to_save ]
                else :
                    if results[k][hd['FirstResHB']]==0 : anti_hb_0[i]+=[ to_save ]
                    else : anti_hb_1[i]+= [ to_save ]
        del results
    if logfile is not None : logfile.write('\n\n')
    if save_coverage_results is not None :
        if type(save_coverage_results) is str : outname=save_coverage_results
        else : outname='protein_coverage_results'+protein[:8]+'.tsv'
        protein_coverage_res.Print(outname)
    if not just_save_coverage_results :
        if outpickle  is not None :
            out=open(outpickle,'wb')
            pickle.dump( (parallel_hb_0, parallel_hb_1 , anti_hb_0 , anti_hb_1),out)
            out.close()
        return parallel_hb_0, parallel_hb_1 , anti_hb_0 , anti_hb_1


def blast_peptide(peptideNtoC, blastDB='/Users/francesco/Documents/proteomes/human_reference_proteome_filtered.fasta',ncpu=3,at_most_N_mismatch=3,update_this_dict=None,verbose=False):
    '''
    it does an ungapped blast search for the peptide in the decalred db
    returns a dictionary whose keys are the number of mismatching residues (calculated on the length of the peptide not of the alignmet)
    and values is a list of blast outputs:
      each blast output is a list with:
      [subject_id (name of the match), identity (percentage), alignment_length, mismatches (number within alignment length), gaps (number within alignemtn length), q_start (start index in the query sequence), q_end (end index), s_start (start index in the matching sequence), s_end (end index), evalue,nident (number of identical amino acids in aligned region), slen (sublect length) , sseq (subject sequence) ]
    '''
    if type(update_this_dict) is dict : results=update_this_dict
    else : results={}
    #HD of each match:
    #[subject_id (name of the match), identity (percentage), alignment_length, mismatches (number within alignment length), gaps (number within alignemtn length), q_start (start index in the query sequence), q_end (end index), s_start (start index in the matching sequence), s_end (end index), evalue,nident (number of identical amino acids in aligned region), slen (sublect length) , sseq (subject sequence) ]
    blast_results=mybio.blast(peptideNtoC,blastDB,ungapped=True,max_target_seqs=100,ncpu=ncpu,add_sseq=True)
    if blast_results is None : 
        if verbose : print("Warning blast_peptide returned empty file when blasting %s" % (peptideNtoC))
        return results
    for res in blast_results :
        non_matching_residues = len(peptideNtoC)-res[2]+res[3] # alignment_length + mismatches
        if at_most_N_mismatch is not None and  non_matching_residues > at_most_N_mismatch : 
            continue
        if non_matching_residues not in results : results[non_matching_residues]=[]
        results[non_matching_residues]+=[res]
    return results



    
def blast_coverage(tuple_coverage_dict,outputfile_overall=None, blastDB='/Users/francesco/Documents/proteomes/human_reference_proteome_filtered.fasta', interesting_mistmatchN=[0,1,2],save_overall_up_to_mismatch=1):
    '''
    tuple_coverage_dict should be (parallel_hb_0, parallel_hb_1 , anti_hb_0 , anti_hb_1 ) as returned by protein_coverage(),
    ALTERNATIVELY tuple_coverage_dict may be a dictionary as returned Complement class
overall, positions, number_of_seq_with_one_hit, tot=beta_strands.blast_coverage(asyn)

abeta_tot+iapp_tot+asyn_tot=15587
number_of_seq_with_one_hit: asyn+abeta+iapp: 22+12+1
P=0.00224546
number_of_seq_with_one_hit (1 mismatch): asyn+abeta+iapp: 1555+980+121
P=0.1703984
    '''
    tot=0
    overall={}
    number_of_seq_with_one_hit={} # like overall but adds one per peptide with at least a hit
    positions={}
    dones=[]
    add_final=[]
    at_most_N_mismatch=max(interesting_mistmatchN)
    for j in interesting_mistmatchN :
        number_of_seq_with_one_hit[j]=0
        if save_overall_up_to_mismatch is not None and save_overall_up_to_mismatch < j: continue
        overall[j]={}
    if isinstance(tuple_coverage_dict,Complement) or 'complementary_dic' in dir(tuple_coverage_dict) : #sometimes for no apparent reason the isinstance call returns False even when a Complement class is passed 
        dicval=list(tuple_coverage_dict.complementary_dic.values())
        hd=tuple_coverage_dict.complementary_hd
        add_final= ['FirstResHB','parallel/antiparallel','Promiscuity_score', 'Count_score', 'Complementarity_score', 'initializer','len_initializer', 'NumberOfPaths']
        for comp in dicval:
            # a comp here looks like [69, 'AITIEFY', 'AITIEFY', '-', 5, 'TIEFY', 1, 1, 39.1, 1.2951354941773072, 43, -3.9] provided that interest_in_range was given
            a=process_comp_and_blast(comp, number_of_seq_with_one_hit, overall, i_seq=hd['Complementary_sequence'],i_direction=hd['parallel/antiparallel'],query_info_list=[ hd[k] for k in add_final],done_list=dones,ncpu=2,at_most_N_mismatch=at_most_N_mismatch,interesting_mistmatchN=interesting_mistmatchN,save_overall_up_to_mismatch=save_overall_up_to_mismatch, blastDB=blastDB)
            if a is None : continue
            tot+=1
    else :
        comp=['']*6
        add_final= ['Position','FirstResHB','parallel/antiparallel','len_initializer','initializer','Complementarity_score', 'Solubility_score' ]
        for j,p in enumerate(tuple_coverage_dict) :
            for k in p : # keys are index in protein length (the protein from which the coverage is calculated)
                if p[k]==[] or p[k] is None : continue
                if k not in positions : positions[k]=numpy.zeros(len(interesting_mistmatchN))
                for comp in p[k] :
                    # a comp here looks like ['EQHAIVLT', '-', 4, 'EQHA', 15.310000000000002, 0.8386152082871657]
                    a=process_comp_and_blast(comp, number_of_seq_with_one_hit, overall, i_seq=0, i_direction=1,query_info_list=[k+1,j%2]+comp[1:],done_list=dones,ncpu=2,at_most_N_mismatch=at_most_N_mismatch,interesting_mistmatchN=interesting_mistmatchN,save_overall_up_to_mismatch=save_overall_up_to_mismatch, blastDB=blastDB)
                    if a is None : continue
                    tot+=1 
                    positions[k]+=a
                if len(comp)+1<len(add_final) : add_final=add_final[:len(comp)+1] # not score candidates
###
                sys.stdout.write( '%d %d;'%(j,k))
                sys.stdout.flush()
    if outputfile_overall is not None :
        out=open(outputfile_overall, 'w')
        out.write('complementary_peptideNtoC\tnumber_of_mismatches\tsubject_id\ts_start\ts_end\tsubject_seq\t'+'\t'.join(add_final)+'\n')
        for nm in overall :
            for c_pep in overall[nm] :
                for comp in overall[nm][c_pep] :
                    out.write('%s\t%d\t%s\n' % (c_pep,nm,'\t'.join(map(str,comp))))
        out.close()
    return overall, positions, number_of_seq_with_one_hit, tot

def process_comp_and_blast(comp, number_of_seq_with_one_hit,overall, i_seq=0,i_direction=1,query_info_list=[],ncpu=2,at_most_N_mismatch=3,positions_dict=None,done_list=None,interesting_mistmatchN=[0,1,2],save_overall_up_to_mismatch=1, blastDB='/Users/francesco/Documents/proteomes/human_reference_proteome_filtered.fasta'):
    '''
    auxiliary function used by blast_coverage
    '''
    if comp[i_direction]=='+' : seq=comp[i_seq]
    else : seq=comp[i_seq][::-1]
    if done_list is not None :
        if seq in done_list : return
        else : done_list+=[seq]
    try :
        r = blast_peptide(seq,blastDB=blastDB,ncpu=ncpu,at_most_N_mismatch=at_most_N_mismatch,update_this_dict=None)
    except Exception :
        sys.stderr.write('\nException raised when blasting seq=%s\n' % (seq))
        print_exc( file=sys.stderr )
        sys.stderr.flush()
        pass
    #    0                                1                    2                    3                                            4                                        5                                        6                    7                                            8                    9    10                                                            11                      12
    #[subject_id (name of the match), identity (percentage), alignment_length, mismatches (number within alignment length), gaps (number within alignemtn length), q_start (start index in the query sequence), q_end (end index), s_start (start index in the matching sequence), s_end (end index), evalue,nident (number of identical amino acids in aligned region), slen (sublect length) , sseq (subject sequence) ]
    a=numpy.zeros(len(interesting_mistmatchN))
    for j in range(len(interesting_mistmatchN)) :
        if j in r : 
            a[j]=len(r[j])
            if len(r[j])>0 : number_of_seq_with_one_hit[j]+=1
            if save_overall_up_to_mismatch is not None and save_overall_up_to_mismatch < j: continue
            # rj=['sp|P76205|ARPB_ECOLI', 100.0, 5, 0, 0, 1, 5, 173, 177, 3.3, 5, 632, 'HLSLQ']
            overall[j][seq] = [ [ rj[0], rj[7], rj[8], rj[12] ]+query_info_list for rj in r[j] ] # r[j] is a list of lists
    return a

def blast_many_peptides(peptides_recs,blastDB='/Users/francesco/Documents/proteomes/ecoli_reference_proteome.fasta',ncpu=4,at_most_N_mismatch=1,out=True):
    if os.path.isfile(peptides_recs) :
        if out==True : out=mybio.get_file_path_and_extension(peptides_recs)[1]+'_blast_results.txt'
        peptides_recs=mybio.get_SeqRecords(peptides_recs)
    data=csv_dict.Data()
    data.hd=csv_dict.list_to_dictionary(['Pep ID','Number of Mismatches','subject_id (name of the match)', 'identity (percentage)', 'alignment_length', 'mismatches (number within alignment length)', 'gaps (number within alignemtn length)', 'q_start (start index in the query sequence)', 'q_end (end index)', 's_start (start index in the matching sequence)', 's_end (end index)', 'evalue', 'nident (number of identical amino acids in aligned region)', 'slen (sublect length)','sseq (subject sequence)' ])
    n=0
    N=len(peptides_recs)
    for j,p in enumerate(peptides_recs) :
        if j%50==0 : 
            sys.stdout.write('| %5.2lf %% |' % (100.*j/N))
            sys.stdout.flush()
        seq=str(p.seq)
        r = blast_peptide(seq,blastDB=blastDB,ncpu=ncpu,at_most_N_mismatch=at_most_N_mismatch,update_this_dict=None)
        # r is a dict with number of mismatches as key and content the following:
        #    0                                1                    2                    3                                            4                                        5                                        6                    7                                            8                    9    10                                                            11                      12
        #[subject_id (name of the match), identity (percentage), alignment_length, mismatches (number within alignment length), gaps (number within alignemtn length), q_start (start index in the query sequence), q_end (end index), s_start (start index in the matching sequence), s_end (end index), evalue,nident (number of identical amino acids in aligned region), slen (sublect length) , sseq (subject sequence) ]
        for MM in r :
            for hit in r[MM] :
                data[n]=[p.id, MM]+hit[:]
                n+=1
    if type(out) is str :
        data.Print(out)
        print("\nOutput printed in %s\n" %(out))
    return data



def plot_coverage(protein, tuple_coverage_dict , plot=True,use_log_plus_one=False,start_rescount=1,second_has_hb_as_negative=False,negatives_according_to_top_bottom_of_seq=True,remove_empties=True,separate_parallel_antiparallel=True,stacked=True,just_count_start=False,min_initializer_len=4,min_complementarity_score=30.,min_Solubility_score=-1.5, insize_ind=2,compscore_ind=4, solubility_score_ind=5,figsize=(14,8),print_all_sequence=True,ym=None,yM=None,y_range=(0,None),save=None,**kwargs):   
    '''
    tuple_coverage_dict should be (parallel_hb_0, parallel_hb_1 , anti_hb_0 , anti_hb_1 ) as returned by protein_coverage() or screen_complementary_peptides(),
     it can also contain only one of more of the above if you are interested in plotting only that specific coverage
mat_count,pl=beta_strands.plot_coverage(mybio.asyn,tup,min_initializer_len=None,min_complementarity_score=None,min_Solubility_score=None,ym=0.5,yM=1.,separate_parallel_antiparallel=True,stacked=False,save='coverage_asyn.png')

# Plot Figure Coverage and get matrices for statistics Figure 2 PNAS
import beta_strands,plotter,mybio,cPickle as pickle
tupsyn,tupbeta,tupiapp = pickle.load(open('coverage_tup_syn_abeta_iapp.pkl','rb'))
mat_cous,pl=beta_strands.plot_coverage(mybio.asyn,tupsyn,min_initializer_len=None,min_complementarity_score=None,min_Solubility_score=None,ym=None,yM=None,separate_parallel_antiparallel=True,stacked=False,save='coverage_asyn.png',y_range=(0,10000),figsize=(22,5.5) )
mat_coua,pl=beta_strands.plot_coverage(mybio.abeta42,tupbeta,min_initializer_len=None,min_complementarity_score=None,min_Solubility_score=None,ym=None,yM=None,separate_parallel_antiparallel=True,stacked=False,save='coverage_abeta42.png',y_range=(0,10000),figsize=(9.3,5.5))
mat_coui,pl=beta_strands.plot_coverage(mybio.iapp,tupiapp,min_initializer_len=None,min_complementarity_score=None,min_Solubility_score=None,ym=None,yM=None,separate_parallel_antiparallel=True,stacked=False,save='coverage_amylin.png',y_range=(0,1000),figsize=(9,5.5) )

# Plot bars of db coverage
import beta_strands,plotter,mybio,csv_dict
plotter.set_publish(small_figure=False)
d20,d20Th,disP,disPTh=csv_dict.Data(),csv_dict.Data(),csv_dict.Data(),csv_dict.Data()
d20.load('complementary_scan_D2D_20noTH.tsv')
d20Th.load('complementary_scan_D2D_20.tsv')
disP.load('complementary_scan_DisProt_noTH.tsv')
disPTh.load('complementary_scan_DisProt.tsv')
M=100.
yerr=None 
if 'DistTo95CI_low_positions_covered' in d20.hd : 
    yerr= [ [ M*d20['0'][d20.hd['DistTo95CI_low_positions_covered']], M*disP['0'][disP.hd['DistTo95CI_low_positions_covered']] ],[ M*d20Th['0'][d20Th.hd['DistTo95CI_low_positions_covered']], M*disPTh['0'][disPTh.hd['DistTo95CI_low_positions_covered']]  ]   ]

plotter.multibar([ [M*d20['0'][0],M*disP['0'][0] ],[M*d20Th['0'][0],M*disPTh['0'][0]] ] ,color=plotter.threecols_palette,xlabels=['',''],hgrid=False,frame=False,y_major_tick_every=M*0.2,y_minor_tick_every=M*0.1,yerr=yerr,figure_size=(10,10))

# plot with separated parallel/antiparallel
M=100.
yerr=None
keys=['positions_covered','antiparallel_positions_covered','parallel_positions_covered']
d2d =[M*d20['0'][d20.hd[k]] for k in keys]
dP  =[M*disP['0'][disP.hd[k]] for k in keys]
cols=[plotter.threecols_palette[0],plotter.iworkpalette[1],plotter.iworkpalette[0]]
plotter.multibar(zip(d2d,dP) ,color=cols,xlabels=['',''], hgrid=False,frame=False,y_major_tick_every=M*0.2,y_minor_tick_every=M*0.1,yerr=yerr,figure_size=(10,8),save='global_DB_coverage_parAnti.png',linewidth=1.5,bar_sep=1.5,x_range=(0.2,7.8))


PLOT SCREEN with complementarity score
vlines=[51.96,50.05,33.17,54.59] # compelementarity score of peptide grafted in the DesAb-1-loop QYSVLIDA,GSVAEQFT,RQYVGLR,TIVSLKF respectively
_=beta_strands.plot_minscore_screen([disP,d20], bootstrap=False, figsize=(14,10),label=['                     \n','                     \n'],ls=':',marker='.',vline=vlines,color=None,ms=15,save='minScore_screen.png')
    '''
    if negatives_according_to_top_bottom_of_seq and second_has_hb_as_negative :
        sys.stderr.write("**WARNING** in plot_coverage() both negatives_according_to_top_bottom_of_seq and second_has_hb_as_negative set to True. Will use negatives_according_to_top_bottom_of_seq!\n")
    if separate_parallel_antiparallel==True and len(tuple_coverage_dict) != 4 :
        print('len(tuple_coverage_dict) != 4  cant separate_parallel_antiparallel, setting to False')
        separate_parallel_antiparallel=False
    mat_count=numpy.zeros( (len(tuple_coverage_dict),len(protein)))
    for i,covtype in enumerate(tuple_coverage_dict) : # covtype is a dictionary (coverage type parallel first HB=0, ''=1, Antiparallel first HB=0, ''HB=1 )
        for j in range(len(protein)) :
            if j in covtype :
                if covtype[j]==[] or covtype[j] is None : continue
                for comp in covtype[j] :
                    if (min_initializer_len is None or comp[insize_ind] >= min_initializer_len) \
                            and (min_complementarity_score is None or comp[compscore_ind] >= min_complementarity_score) \
                            and (min_Solubility_score is None or comp[solubility_score_ind] >= min_Solubility_score) :
                        if negatives_according_to_top_bottom_of_seq :
                            # in this case the 4 categories become parallel top, parallel bottom, antiparallel top, antiparallel bottom.
                            # the underlying assumpution is that the firs residue of the protein will make HB towards the top.
                            if i<2 : # parallel, i%2=1 means first HB=1
                                ii=abs(j%2-i%2) # the facing aa of a parallel complementary with HB=1 will NOT be involved in HB.
                            else : # antiparallel, i%2=1 means first HB=1
                                ii=2+abs(1-abs(j%2-i%2)) # the facing aa of a antiparallel complementary with HB=1 will  be involved in HB.
                            #if ii%2==0 : print comp
                            if just_count_start : mat_count[ii][j]+=1
                            else : mat_count[ii][j:j+len(comp[0])]+=1
                        else :
                            if just_count_start : mat_count[i][j]+=1
                            else : mat_count[i][j:j+len(comp[0])]+=1
    if plot :
        import plotter
        plotter.set_publish(small_figure=False,thick_ticks=False)
        #if use_log_plus_one :
        #    if  separate_parallel_antiparallel :
        #        if not stacked : pl=  [    numpy.log10(mat_count[:2,:].sum(axis=0) +1.) ,   numpy.log10(mat_count[2:,:].sum(axis=0) +1.) ]
        #        else : pl=  [    numpy.log10(mat_count[:2,:].sum(axis=0) +1.) , numpy.log10(mat_count.sum(axis=0) +1.)  - numpy.log10(mat_count[:2,:].sum(axis=0) +1.) ]
        #    else :pl= numpy.log10( mat_count.sum(axis=0) +1. )
        annotation_string=None
        if  separate_parallel_antiparallel :
            if negatives_according_to_top_bottom_of_seq :
                if y_range is not None and y_range[0]==0 :
                    if y_range[1] is None : y_range=None 
                    else : y_range=(-y_range[1],y_range[1])
                pl=[ mat_count[0,:],-1*mat_count[1,:], mat_count[2,:],-1*mat_count[3,:] ]
                legend=['Parallel top','Parallel bottom','Antiparallel top','Antiparallel bottom']
                annotation_string=(' :'*len(protein))[:len(protein)]
                color=[ plotter.iworkpalette_dict['blue'][0],plotter.iworkpalette_dict['blue'][1],plotter.iworkpalette_dict['green'][0],plotter.iworkpalette_dict['green'][1] ]
            elif second_has_hb_as_negative   :
                if y_range is not None and y_range[0]==0 :
                    if y_range[1] is None : y_range=None 
                    else : y_range=(-y_range[1],y_range[1])
                pl=[ mat_count[0,:],-1*mat_count[1,:], mat_count[2,:],-1*mat_count[3,:] ]
                legend=['Parallel HB0','Parallel HB1','Antiparallel HB0','Antiparallel HB1']
                color=[ plotter.iworkpalette_dict['blue'][0],plotter.iworkpalette_dict['blue'][1],plotter.iworkpalette_dict['green'][0],plotter.iworkpalette_dict['green'][1] ]
            else : 
                pl=[ mat_count[:2,:].sum(axis=0) , mat_count[2:,:].sum(axis=0) ]
                legend=['Parallel','Antiparallel']
                color=[ plotter.iworkpalette_dict['blue'][0],plotter.iworkpalette_dict['green'][0] ]
            if remove_empties :
                rem=[]
                for j,set in enumerate(pl) :
                    if all(set==0) : rem+=[j]
                if rem!=[] :
                    for j in sorted(rem,reverse=True) : del pl[j],color[j],legend[j]
        else : 
            pl= mat_count.sum(axis=0)
            legend=None
            color=plotter.iworkpalette[0]
        if use_log_plus_one : 
            log_scale=True
        else : 
            log_scale=False
        plotter.plot_seq_profile(protein, pl, bar=True,start_rescount=start_rescount,annotation_string=annotation_string,label=legend,log_scale=log_scale,stacked=stacked,bar_sep=0, xlabel=None, ylabel=None, title=None, zygg_like_lines=False, hline=False, vgrid=False, print_all_sequence=print_all_sequence,y_range=y_range, y_major_tick_every=yM, y_minor_tick_every=ym,figure_size=figsize,color=color, save=save,legend_location='lower right',legend_size=14,**kwargs)
    return mat_count,numpy.array(pl)




def visualize_coverage(protein, tuple_coverage_dict, outfile=sys.stdout,out_as_initializer_sizes=False, out_as_best_seq=True, merge_initializer_at_the_end=True,out_as_first_seq=False, split_hb_above_and_below=True, split_every=None,min_initializer_len=4,start_rescount=1, debug=False):
    '''
    Pretty prints a visualization of the protein coverage (in pure ascii so don't expect to get overwelmed by beauty!)
    tuple_coverage_dict should be (parallel_hb_0, parallel_hb_1 , anti_hb_0 , anti_hb_1 ) as returned by protein_coverage() or a pickle containg them
     it writes to outfile a "not too ugly" representation of the coverage
    KNOWN BUGS: the representation of parallel fragments with split_hb_above_and_below is probably opposite as one would expect (aka wrong)
_=beta_strands.visualize_coverage( Nasyn,(parallel_hb_0, parallel_hb_1 , anti_hb_0 , anti_hb_1), outfile=sys.stdout)
    '''
    
    if out_as_first_seq or out_as_best_seq: 
        out_as_initializer_sizes=False
        merge_initializer_at_the_end=True
    if type(tuple_coverage_dict) is str :
        p0,p1,a0,a1=pickle.load(open(tuple_coverage_dict,'rb'))
    else :
        p0,p1,a0,a1=tuple_coverage_dict
    p0_str={}
    a0_str={}
    p1_str={}
    a1_str={}
    if out_as_initializer_sizes : 
        in_sizesp0=[]
        in_sizesp1=[]
        in_sizesa0=[]
        in_sizesa1=[]
    plen=len(protein)
    if debug : print(' |Deb visualise_coverage: about to fill 4 str[i], i ranges up to plen:',plen)
    for i in range(plen) :
        if debug : print(' |i=',i)
        if debug and i in p0 : print('    p0=',len(p0[i]))
        p0_str[i] = repr_coverage_dic_entry(p0,i, is_parallel=True ,out_as_initializer_sizes=out_as_initializer_sizes,out_as_best_seq=out_as_best_seq,out_as_first_seq=out_as_first_seq,initializer_size_index=2, min_initializer_len=min_initializer_len,debug=False)
        if debug and i in p1 : print('    p1=',len(p1[i]))
        p1_str[i] = repr_coverage_dic_entry(p1,i, is_parallel=True ,out_as_initializer_sizes=out_as_initializer_sizes,out_as_best_seq=out_as_best_seq,out_as_first_seq=out_as_first_seq,initializer_size_index=2, min_initializer_len=min_initializer_len,debug=False)
        if debug and i in a0 : print('    a0=',len(a0[i]))
        a0_str[i] = repr_coverage_dic_entry(a0,i, is_parallel=False,out_as_initializer_sizes=out_as_initializer_sizes,out_as_best_seq=out_as_best_seq,out_as_first_seq=out_as_first_seq,initializer_size_index=2, min_initializer_len=min_initializer_len,debug=False)
        if debug and i in a1 : print('    a1=',len(a1[i]))
        a1_str[i] = repr_coverage_dic_entry(a1,i, is_parallel=False,out_as_initializer_sizes=out_as_initializer_sizes,out_as_best_seq=out_as_best_seq,out_as_first_seq=out_as_first_seq,initializer_size_index=2, min_initializer_len=min_initializer_len,debug=False)
        
        if out_as_initializer_sizes :
            if type(p0_str[i]) is dict : in_sizesp0+=list(p0_str[i].keys()) # add it only if there was at least one complementary for this dict at this position
            if type(p1_str[i]) is dict : in_sizesp1+=list(p1_str[i].keys())
            if type(a0_str[i]) is dict : in_sizesa0+=list(a0_str[i].keys())
            if type(a1_str[i]) is dict : in_sizesa1+=list(a1_str[i].keys())
    if debug : print('  Deb visualise_coverage: constructed 4 a0_str[i]')
    if out_as_initializer_sizes : 
        in_sizesp0=misc.uniq(in_sizesp0)
        in_sizesp1=misc.uniq(in_sizesp1)
        in_sizesa0=misc.uniq(in_sizesa0)
        in_sizesa1=misc.uniq(in_sizesa1)
        in_sizesp0.sort(reverse=True)
        in_sizesp1.sort(reverse=True)
        in_sizesa0.sort(reverse=True)
        in_sizesa1.sort(reverse=True)
        if not merge_initializer_at_the_end :
            # in this case we revert the dictionaries of initializers key, so that initializer sizes becomes the outer key (now it is the inner).
            p0_str_ini={}
            for i_s in in_sizesp0 :p0_str_ini[i_s]={}
            for i in range(plen) :
                for i_s in in_sizesp0 :
                    if type(p0_str[i]) is dict and i_s in p0_str[i]:
                        p0_str_ini[i_s][i]=p0_str[i][i_s]
                    else : p0_str_ini[i_s][i]=' '
            p1_str_ini={}
            for i_s in in_sizesp1 :p1_str_ini[i_s]={}
            for i in range(plen) :
                for i_s in in_sizesp1 :
                    if type(p1_str[i]) is dict and i_s in p1_str[i]:
                        p1_str_ini[i_s][i]=p1_str[i][i_s]
                    else : p1_str_ini[i_s][i]=' '
            a0_str_ini={}
            for i_s in in_sizesa0 :a0_str_ini[i_s]={}
            for i in range(plen) :
                for i_s in in_sizesa0 :
                    if type(a0_str[i]) is dict and i_s in a0_str[i]:
                        a0_str_ini[i_s][i]=a0_str[i][i_s]
                    else : a0_str_ini[i_s][i]=' '
            a1_str_ini={}
            for i_s in in_sizesa1 :a1_str_ini[i_s]={}
            for i in range(plen) :
                for i_s in in_sizesa1 :
                    if type(a1_str[i]) is dict and i_s in a1_str[i]:
                        a1_str_ini[i_s][i]=a1_str[i][i_s]
                    else : a1_str_ini[i_s][i]=' '
            p0_lines={}
            for i_s in in_sizesp0 :
                p0_lines[i_s]=generate_coverage_line(p0_str_ini[i_s],plen)
            p1_lines={}
            for i_s in in_sizesp1 :
                p1_lines[i_s]=generate_coverage_line(p1_str_ini[i_s],plen)
            a0_lines={}
            for i_s in in_sizesa0 :
                a0_lines[i_s]=generate_coverage_line(a0_str_ini[i_s],plen)
            a1_lines={}
            for i_s in in_sizesa1 :
                a1_lines[i_s]=generate_coverage_line(a1_str_ini[i_s],plen)
        else :
            for i in range(plen) :
                if type(p0_str[i]) is dict :
                    p0_str[i]=list(p0_str[i].values())[0]
                if type(p1_str[i]) is dict :
                    p1_str[i]=list(p1_str[i].values())[0]
                if type(a0_str[i]) is dict :
                    a0_str[i]=list(a0_str[i].values())[0]
                if type(a1_str[i]) is dict :
                    a1_str[i]=list(a1_str[i].values())[0]
            
            
            #p0_lines=generate_coverage_line(p0_str,plen)
            #p1_lines=generate_coverage_line(p1_str,plen)
            #a0_lines=generate_coverage_line(a0_str,plen)
            #a1_lines=generate_coverage_line(a1_str,plen)
    #else :
    if split_hb_above_and_below :
        if debug : print('  Deb visualise_coverage: about to construct 4 a0_lines[i] with split_hb_above_and_below')
        p_up={}
        p_down={}
        a_up={}
        a_down={}
        line_up=''
        line_down=''
        for i in range(0,plen,2) :
            line_up+=': '
            a_up[i]=a1_str[i]
            p_up[i]=p1_str[i]
            line_down+=' :'
            a_down[i]=a0_str[i]
            p_down[i]=p0_str[i]
            if i+1<plen :
                a_up[i+1]=a0_str[i+1]
                p_up[i+1]=p0_str[i+1]
                a_down[i+1]=a1_str[i+1]
                p_down[i+1]=p1_str[i+1]
        if len(line_up)>plen : line_up=line_up[:-1]
        if len(line_down)>plen : line_down=line_down[:-1]
        p0_lines=generate_coverage_line(p_up,plen)
        a0_lines=generate_coverage_line(a_up,plen)
        a1_lines=generate_coverage_line(a_down,plen)
        p1_lines=generate_coverage_line(p_down,plen)
        if debug : print('  Deb visualise_coverage: constructed 4 a0_lines[i] with split_hb_above_and_below')
    else :
        a1_lines=generate_coverage_line(a1_str,plen)
        a0_lines=generate_coverage_line(a0_str,plen)
        p0_lines=generate_coverage_line(p0_str,plen)
        p1_lines=generate_coverage_line(p1_str,plen)

    #if split_every>0:
    if debug : print('  Deb visualise_coverage: about to fill wholefile') 
    wholefile='' 
    close_out=False
    if type(outfile) is str :
        outfile=open(outfile,'w')
        close_out=True
    if not out_as_initializer_sizes or merge_initializer_at_the_end :
        ts='        '+'%-9d' % (start_rescount) # line with numbers to count positions along the sequence
        for i in range(start_rescount+9,start_rescount+plen,10) :ts+='%-10d' % (i)
        ts+='\n'
        wholefile += ts
        sep='-------<%s|\n' % ('-~-~-~-~-|'*(len(protein)/10)+'-~-~-~-~-|'[:len(protein)%10] ) # separator line with ticks represented by - or ~
        #sep+='\n'
        wholefile+=sep
        for l in p0_lines :
            if split_hb_above_and_below : wholefile+='paralle<'+l+'|\n' #wholefile+=('paralle<%s|' % (l) ).strip()+'\n' # if the sting l is formatted empty spaces will be added afterwards, unfortunately strip does not solve the problem..
            else : wholefile+=('paral_0<%s|' % (l)).strip()+'\n'
        for l in a0_lines :
            if split_hb_above_and_below : wholefile+='antipar<'+l+'|\n'# wholefile+=('antipar<%s|' % (l)).strip()+'\n'
            else : wholefile+=('antip_0<%s|' % (l)).strip()+'\n'

        if split_hb_above_and_below : wholefile+=('hb_patt<%s|' % (line_up)).strip()+'\n' # hb pattern 
        wholefile+='protein>%s|\n' % (protein)  # target protein sequence
        if split_hb_above_and_below : wholefile+=('hb_patt<%s|' % (line_down)).strip()+'\n' # hb pattern
        del l
        for l in a1_lines :
            if split_hb_above_and_below : wholefile+='antipar<'+l+'|\n' # wholefile+=('antipar<%s|' % (l)).strip()+'\n'
            else : wholefile+=('antip_1<%s|' % (l)).strip()+'\n'
        for l in p1_lines :
            if split_hb_above_and_below : wholefile+='paralle<'+l+'|\n'#  wholefile+=('paralle<%s|' % (l)).strip()+'\n'
            else : wholefile+=('paral_1<%s|' % (l)).strip()+'\n'
        wholefile+=sep
        wholefile += ts
        if split_every is not None and split_every>0 : # if Fake_str are use to print colors the split_every will have many problems and won't work
            wholefile=misc.split_file_in_blocks(wholefile,split_every)
        outfile.write(wholefile)
    else :
        print('**** WARN**** in visualize_coverage() print of (out_as_initializer_sizes or not merge_initializer_at_the_end) NOT IMPLEMENTED.\n ')
    if close_out:
        outfile.close()
    
    
    return p0_lines,p1_lines,a0_lines,a1_lines

        
class Color:
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'         

class Fake_seq :
    # mimic of the Biopython SeqRecord object
    def __init__(self,seq='',seq_id='',seq_name='',description='') :
        self.seq=seq
        self.id=seq_id
        self.name=seq_name
        self.description=description
    def __len__(self) :
        return len(self.seq)
    def __str__(self):
        return self.seq
    def __getslice__(self,i,j):
        return self.seq[i:j]
    def __getitem__(self,y):
        return self.seq[y]
    def __add__(self,y):
        return self.seq+y
    def __iadd__(self,y):
        self.seq+=y
    
class Fake_str(str) :
    def __init__(self, *argv, **argd):
        self.length=None
        super(Fake_str, self).__init__(*argv, **argd) # this inits the str
    def __len__(self):
        if self.length is not None :
            return self.length
        else : return len(str(self))


    
def generate_coverage_lineOLD(position_dic, protein_len):
    '''
    generates a line of the coverage representation (for a given protein)
    '''
    lines=[Fake_str('')]
    lines[-1].length=0
    j=0
    #print position_dic
    available={j:True} # this dictionary will tell for each row (initially only one row in line) whether the position under scrutiny (i.e. i index in the for below) is available to add the complementary sequence. 
    #                    If no row has the position available then a new row is created
    for i in range(protein_len) :
        j=0
        #done_on=None
        while True :
            if j not in available :
                #print len(lines),i,len(position_dic[i])
                lines+=[Fake_str( (i)*' ' )] # this fills with spaces up to index i excluded
                lines[-1].length = i
                available[j]=True
            if available[j] :
                oldlen=len(lines[j])
                lines[j]=Fake_str(lines[j] + position_dic[i])
                lines[j].length =  oldlen+len(position_dic[i]) # which might be a Fake_str
                if len(position_dic[i])>1 : available[j]=False
                #done_on=j
                break
            j+=1
        # reset the availables
        for k in available:
            if len(lines[k])<=i+1 :
                available[k]=True
                oldlen=len(lines[k])
                lines[k] = Fake_str( lines[k] + (i+1-oldlen)*' ' )
                lines[k].length= i+1
    return lines


def generate_coverage_line(position_dic, protein_len, allow_empty_space=True):
    '''
    generates a line of the coverage representation (for a given protein)
    position_dic can have as value a list in cases where more than one sequence should be printed at a given position
    '''
    lines=[Fake_str('')]
    lines[-1].length=0
    j=0
    #print position_dic

    available_lines=[0] #  this list will tell for each row (row index is the key of the dictionary, initially only one row in line) 
    #                    whether the position under scrutiny (i.e. i index in the for below) is available to add the complementary sequence. 
    #                    If no row has the position available then a new row is created
    for i in range(protein_len) :
        if i  in position_dic and position_dic[i]!=' ' and position_dic[i]!='' : # add a space to available one and update availability of other lines
            #done_on=None
            if type(position_dic[i]) is not list : pos_i= [ position_dic[i] ] # so that position_dic does not get changed
            else : pos_i= position_dic[i]
            for seq_i in pos_i :
                if len(available_lines)==0 :
                    lines+=[Fake_str( (i)*' ' )] # this adds a new line as none of the existing ones is available (it fills with spaces up to index i excluded)
                    lines[-1].length = i
                    available_lines+=[ len(lines)-1 ]
                    
                j=available_lines[0]
                oldlen=len(lines[j])
                if i+len(seq_i)==protein_len : # no extra empty space
                    lines[j]=Fake_str(lines[j] + seq_i)
                    lines[j].length =  oldlen+len(seq_i)
                else :
                    lines[j]=Fake_str(lines[j] + seq_i+' '*int(allow_empty_space))
                    lines[j].length =  oldlen+len(seq_i)+int(allow_empty_space) # seq_i may be a Fake_str
                available_lines.remove(j) # no longer available at this postion
                
        # add a space to available one and update availability of other lines
        for j in range(len(lines)) :
            if j in available_lines :
                oldlen=len(lines[j])
                lines[j] = Fake_str( lines[j] +' ')
                lines[j].length=oldlen+1
            elif lines[j].length<=i+1 : 
                available_lines+=[j]
    
    #print "DEB: ",len(lines),[len(l) for l in lines]
    return lines


def repr_coverage_dic_entry(cov_dict,index,is_parallel=False,out_as_initializer_sizes=True, out_as_best_seq=False,out_as_first_seq=False,return_one_per_initializer_of_at_least=5,bold_initializers_greater_than=4,initializer_size_index=2,min_initializer_len=4,min_Solubility_score=-1.2,debug=False):   
    '''
    auxiliary function that represents (convert to string) a single value of the coverage dictionary
    '''
    if index in cov_dict :
        if cov_dict[index]==[] : ch=' '
        else :
            chosen_seq=None
            complementary_len=len(cov_dict[index][0][0])
            #if len(cov_dict[index])<4 : print cov_dict[index]
            if out_as_best_seq : # this will have an underlined initializer and a bold start font
                #comps=sorted( cov_dict[index] , key=lambda x : x[4], reverse=True) # sort according to complementarity score
                ch, chosen = _aux_out_as_best_seq(cov_dict, index, initializer_size_index=initializer_size_index , min_initializer_len=min_initializer_len, min_Solubility_score=min_Solubility_score,bold_initializers_greater_than=bold_initializers_greater_than,already_processed_initialisers=None)
                if debug : print('         _aux_out_as_best_seq first returned')
                if return_one_per_initializer_of_at_least > 0 :
                    if return_one_per_initializer_of_at_least < min_initializer_len :
                        sys.stderr.write('**Warning in repr_coverage_dic_entry() return_one_per_initializer_of_at_least < min_initializer_len %d %d\n' % (return_one_per_initializer_of_at_least,min_initializer_len))
                        return_one_per_initializer_of_at_least=min_initializer_len
                    dones, ch_res = [chosen[initializer_size_index+1]],[ch]
                    while ch!='' :
                        ch,chosen= _aux_out_as_best_seq(cov_dict, index, initializer_size_index=initializer_size_index , min_initializer_len=return_one_per_initializer_of_at_least, min_Solubility_score=min_Solubility_score,bold_initializers_greater_than=bold_initializers_greater_than,already_processed_initialisers=dones)
                        if ch=='' or chosen[initializer_size_index+1] in dones: break
                        dones+=[ chosen[initializer_size_index+1].upper()]
                        if debug : print('           ch=',ch,chosen,chosen[initializer_size_index+1],chosen[initializer_size_index+1] in dones)
                        ch_res+=[ch]
                    if ch_res==[''] : return ' '
                    else : return ch_res
            elif out_as_initializer_sizes :
                iniz={}
                for comp in cov_dict[index] :
                    if comp[initializer_size_index] not in iniz : iniz[comp[initializer_size_index]]=0
                    iniz[comp[initializer_size_index]]+=1
                comp_str={}
                for k in iniz :
                    num=str(iniz[k])
                    in_size=str(k)
                    off=len(num)+len(in_size)
                    if complementary_len-off>= 4: # we add 4 char  [] and -i
                        c_str=num+'-i'+in_size
                        q=(complementary_len-off-4)/2
                        r=(complementary_len-off-4)%2
                        c_str='['+'-'*q+c_str+'-'*(q+r)+']'
                    else :
                        c_str='?'*complementary_len
                        print('****WARN**** NOT IMPLEMENTED in repr_coverage_dic() complementary_len-off< 4')
                    comp_str[k]=c_str
                ch=comp_str
            elif out_as_first_seq : # this will have an underlined initializer and a bold start font
                comp= cov_dict[index][0] # get first sequence
                if chosen_seq is None : chosen_seq=comp[:]
                isize=comp[initializer_size_index]
                if chosen_seq[initializer_size_index+1] not in chosen_seq[0] : # not found, hopefully means that interest_in_range was used and initializers falls outside the range.
                    if chosen_seq[initializer_size_index+1][1:] in chosen_seq[0] :
                        chosen_seq[initializer_size_index+1]=chosen_seq[initializer_size_index+1][1:]
                        if isize > len(chosen_seq[initializer_size_index+1]): # with most option the size would have already been adjusted when generating the complementary petpides using interest_in_range
                            chosen_seq[initializer_size_index]=isize-1
                            isize-=1
                    elif chosen_seq[initializer_size_index+1][:-1] in chosen_seq[0] :
                        chosen_seq[initializer_size_index+1]=chosen_seq[initializer_size_index+1][:-1]
                        if isize > len(chosen_seq[initializer_size_index+1]): # with most option the size would have already been adjusted when generating the complementary petpides using interest_in_range
                            chosen_seq[initializer_size_index]=isize-1
                            isize-=1
                    else :
                        sys.stderr.write("**WARNING** in repr_coverage_dic_entry() intializer %s not found within comp_seq %s. NOT REPRESENTING CORRESPONDING REGION (index %d) IN COVERAGE\n" % (chosen_seq[initializer_size_index+1], chosen_seq[0],index))
                        sys.stderr.flush()
                        return ''
                istart= comp[0].find(comp[initializer_size_index+1]) # initializer_size_index+1  is the initializer sequence. comp[0] the full sequence
                if istart==0 :
                    ch = Fake_str( Color.RED+Color.UNDERLINE+comp[initializer_size_index+1][0]+Color.END+Color.UNDERLINE+comp[initializer_size_index+1][1:]+Color.END+comp[0][isize:] )
                    ch.length = len(comp[0]) 
                else : 
                    init=Color.UNDERLINE+comp[initializer_size_index+1]+Color.END
                    ch = Fake_str( Color.RED+comp[0][0]+Color.END+comp[0][1:istart]+init+comp[0][istart+isize:] )
                    ch.length = len(comp[0])
            else :            
                num=str(len(cov_dict[index]))
                off=len(num)
                if complementary_len-off>= 2: # we add 2 char  []
                    c_str=num
                    q=(complementary_len-off-2)/2
                    r=(complementary_len-off-2)%2
                    c_str='['+'-'*q+c_str+'-'*(q+r)+']'
                else :
                    c_str='?'*complementary_len
                    print('****WARN**** NOT IMPLEMENTED in repr_coverage_dic() complementary_len-off< 4')
                ch=c_str
    else : ch=' '
    return ch
  
def _aux_out_as_best_seq(cov_dict, index, initializer_size_index , min_initializer_len=None, min_Solubility_score=None,bold_initializers_greater_than=4, already_processed_initialisers=None):
    '''
    auxiliary_function used by repr_coverage_dict
    '''
    #comps=sorted( cov_dict[index] , key=lambda x : x[4], reverse=True) # sort according to complementarity score
    comps=sorted( cov_dict[index] , key=lambda x : x[5], reverse=True) # sort according to solubility score
    comps=sorted( comps , key=lambda x : x[initializer_size_index], reverse=True) # sort according to initializer size
    chosen_seq=None
    for c in comps :
        if ( min_initializer_len is None or ((type(min_initializer_len) is dict and min_initializer_len[c[1]] <= c[initializer_size_index]) or (type(min_initializer_len) is int and min_initializer_len <= c[initializer_size_index]) ) )\
        and (min_Solubility_score is None or c[5]>=min_Solubility_score) : # c[1] is +/-
            #print c[initializer_size_index+1] , already_processed_initialisers
            if already_processed_initialisers is None or already_processed_initialisers==[] or c[initializer_size_index+1].upper() not in already_processed_initialisers :
                chosen_seq=c[:]
                break
    if chosen_seq is None :
        if already_processed_initialisers is None or already_processed_initialisers==[] : 
            chosen_seq=comps[0][:]
            chosen_seq[0]=chosen_seq[0].lower() # make the complementary sequence lower case as it does not satisfy the thresholds...
            chosen_seq[initializer_size_index+1]=chosen_seq[initializer_size_index+1].lower() # initializer_size_index+1  is the initializer sequence. comp[0] the full sequence
        else :
            return '',None # we have finished the initializers of suitable quality
        
    isize=chosen_seq[initializer_size_index]
    if chosen_seq[initializer_size_index+1] not in chosen_seq[0] : # not found, hopefully means that interest_in_range was used and initializers falls outside the range.
        if chosen_seq[initializer_size_index+1][1:] in chosen_seq[0] :
            chosen_seq[initializer_size_index+1]=chosen_seq[initializer_size_index+1][1:]
            if isize > len(chosen_seq[initializer_size_index+1]): # with most option the size would have already been adjusted when generating the complementary petpides using interest_in_range
                chosen_seq[initializer_size_index]=isize-1
                isize-=1
        elif chosen_seq[initializer_size_index+1][:-1] in chosen_seq[0] :
            chosen_seq[initializer_size_index+1]=chosen_seq[initializer_size_index+1][:-1]
            if isize > len(chosen_seq[initializer_size_index+1]): # with most option the size would have already been adjusted when generating the complementary petpides using interest_in_range
                chosen_seq[initializer_size_index]=isize-1
                isize-=1
        else :
            sys.stderr.write("**WARNING** in repr_coverage_dic_entry() intializer %s not found within comp_seq %s. NOT REPRESENTING CORRESPONDING REGION (index %d) IN COVERAGE\n" % (chosen_seq[initializer_size_index+1], chosen_seq[0],index))
            sys.stderr.flush()
            return '',None
    istart= chosen_seq[0].find(chosen_seq[initializer_size_index+1]) # initializer_size_index+1  is the initializer sequence. comp[0] the full sequence
    #if chosen_seq[0][:2]=='LF' :
    #    print '##  ',chosen_seq[0], chosen_seq[initializer_size_index+1],isize,chosen_seq[initializer_size_index]
    if istart==0 :
        if bold_initializers_greater_than is not None and isize>bold_initializers_greater_than : init=Color.RED+Color.BOLD+Color.UNDERLINE+chosen_seq[initializer_size_index+1][0]+Color.END+Color.UNDERLINE+Color.BOLD+chosen_seq[initializer_size_index+1][1:]+Color.END
        else : init=Color.RED+Color.UNDERLINE+chosen_seq[initializer_size_index+1][0]+Color.END+Color.UNDERLINE+chosen_seq[initializer_size_index+1][1:]+Color.END
        ch = Fake_str( init +chosen_seq[0][isize:] )
        ch.length = len(chosen_seq[0])
    else : 
        if bold_initializers_greater_than is not None and isize>bold_initializers_greater_than : init=Color.UNDERLINE+Color.BOLD+chosen_seq[initializer_size_index+1]+Color.END 
        else : init=Color.UNDERLINE+chosen_seq[initializer_size_index+1]+Color.END
        ch = Fake_str( Color.RED+chosen_seq[0][0]+Color.END+chosen_seq[0][1:istart]+init+chosen_seq[0][istart+isize:] )
        ch.length = len(chosen_seq[0])
    return ch, chosen_seq
  






class Strand :
    def __init__(self,sequence=None,direction=None,hb_pattern=None):
        self.sequence=sequence
        self.hb_pattern=hb_pattern
        self.direction=direction
    def __str__(self):
        return str(self.sequence)+' '+str(self.direction)+' HB: '+str(self.hb_pattern)
    def __repr__(self):
        return self.__str__()

class BuildSheet :
    '''
    this class builds a strand between two given other strands (here top_strand and bottom_strand) in 
     order to construct a beta sheet of three strands. if bottom_strands is none it just builds complementary peptides...
B=beta_strands.BuildSheet( ('','',() ), ('','',()  ), C)
B.make()
    '''
    def __init__(self,top_strand,bottom_strand=None,complement_class=None,interested_in_range_top=None,interested_in_range_bottom=None):
        if type(top_strand) is tuple or type(top_strand) is list : self.top_strand=Strand(*top_strand)
        else : self.top_strand=top_strand
        if type(bottom_strand) is tuple or type(bottom_strand) is list : self.bottom_strand=Strand(*bottom_strand)
        else : self.bottom_strand=bottom_strand
        if complement_class is None : self.C=Complement()
        else : self.C=complement_class
        self.interested_in_range_top=interested_in_range_top
        self.interested_in_range_bottom=interested_in_range_bottom
        self.logfile=sys.stdout
        self.warn_file=sys.stderr
        self.out_file_composition='sheet_matching_compostion.txt'
        self.out_table_filename='sheet_matching_peptides.txt'
    def make(self,debug=False,score_candidates=True):
        '''
         NB at the moment to save on coding we build all possible peptides - this is highly suboptimal as only those of relevant direction and pattern should be built.
        '''
        self.top_peptides_satisfying_condition,self.top_peptides_composition_dic =  self.screen_one_strand(self.top_strand,score_candidates=score_candidates,debug=debug)
        if self.bottom_strand is not None : 
            self.bottom_peptides_satisfying_condition,self.bottom_peptides_composition_dic =  self.screen_one_strand(self.bottom_strand,score_candidates=score_candidates,debug=debug)
            matches=self.get_merging_keys(self.top_peptides_satisfying_condition, self.bottom_peptides_satisfying_condition)
            if matches!=[] :
                #if type(self.out_file_composition) is str : 
                outc=open(self.out_file_composition,'w')
                #else : outc=self.out_file_composition
                self.peptides_satisfying_condition=csv_dict.Data()
                self.peptides_satisfying_condition.hd=self.top_peptides_satisfying_condition.hd.copy()
                self.peptides_satisfying_condition._update_hd('top/bottom', 0)
                self.peptides_satisfying_condition._update_hd('Id top/bottom', 0)
                j=0
                for kt,kb in matches :
                    self.peptides_satisfying_condition[j]=[kt,'top']+self.top_peptides_satisfying_condition[kt]
                    outc.write('> %d top\n%s\n' % (j, self.top_peptides_composition_dic[kt]))
                    j+=1
                    self.peptides_satisfying_condition[j]=[kb,'bottom']+self.bottom_peptides_satisfying_condition[kb]
                    outc.write('> %d bottom\n%s\n' % (j, self.bottom_peptides_composition_dic[kb]))
                    j+=1
                self.peptides_satisfying_condition.Print(self.out_table_filename)
                if self.logfile is not None : self.logfile.write('found %d sheet peptides, compostion saved in %s and output table in %s' % (j/2, self.out_file_composition,self.out_table_filename))
        return
    def get_merging_keys(self, group1,group2):
        matches=[]
        for k1 in group1 :
            for k2 in group2 :
                if 1 <= self.merging_condition(group1[k1][group1.hd['Complementary_sequence']],group1[k1][group1.hd['Complementary_sequence']]) :
                    matches+=[ (k1,k2) ]
        return matches
    def merging_condition(self,sequence1,sequence2):
        '''
        at the moment we require these to be equal, but a blosum-based seq similarity would be more appropriate (or instead of blosum a purely phys-chem matrix).
        '''
        return sequence1==sequence2
    def screen_one_strand(self,target_strand, debug=False,score_candidates=True):
        peptides_satisfying_condition=csv_dict.Data() # keys will be index
        rhd=['Target Sequence', 'parallel/antiparallel','FirstResHB','Complementary_sequence','len_initializer','initializer','NumberOfPaths','Solubility_score']
        if score_candidates: rhd+=['Complementarity_score','Count_score','Promiscuity_score','Mean_enrichment']
        peptides_satisfying_condition.hd=csv_dict.list_to_dictionary(rhd)
        peptides_composition_dic=OrderedDict()
        tot=0
        
        tot_matching=0
        # get complementaries of top_strand
        target=target_strand.sequence
        if self.logfile is not None : 
            self.logfile.write('\nlooking for complementary peptides of length %d\n' % (len(target)))
        
        self.C.complementary_dic={}
        self.C.composition_dic={}
        self.C.interested_in_range=self.interested_in_range_top
        self.C.make(target ) # run the cascade method
        results=self.C.complementary_dic
        hd=self.C.complementary_hd
        # filter out only the direction and hb pattern requested
        if len(results)==0 :
            if self.warn_file is not None : self.warn_file.write( '****** NO COMPLEMENTS for top_strand\n')
            del results
        else :
            if target_strand.hb_pattern is not None : 
                if target_strand.hb_pattern[0]==(None,None) or target_strand.hb_pattern==[None,None]: FirstResHB=False
                else : FirstResHB=True
                condition = lambda x,cd,hd : cd[x][hd['direction']]==target_strand.direction and cd[x][hd['FirstResHB']]==FirstResHB
            else : condition = lambda x,cd,hd : cd[x][hd['direction']]==target_strand.direction
            keys_matching_condition = screen_complementary_dic(self.C, condition,debug=False)
            if debug : print('top','len(keys_matching_condition)',len(keys_matching_condition),'len(results)',len(results),type(list(results.values())[0] ),list(results.values())[0],hd)
            if len(keys_matching_condition)>0 : 
                tot_matching+=len(keys_matching_condition)
                
                for k in keys_matching_condition :
                    tot+=1
                    # hd should contain: 
                    # complementary_hd=['Target Sequence', 'parallel/antiparallel','FirstResHB','Complementary_sequence','len_initializer','initializer','NumberOfPaths','Solubility_score']
                    peptides_satisfying_condition[tot]= [ target , results[k][hd['parallel/antiparallel']],results[k][hd['FirstResHB']],results[k][hd['Complementary_sequence']],results[k][hd['len_initializer']],results[k][hd['initializer']],results[k][hd['NumberOfPaths']],results[k][hd['Solubility_score']] ]
                    if score_candidates : peptides_satisfying_condition[tot] += [ results[k][hd['Complementarity_score']],results[k][hd['Count_score']],results[k][hd['Promiscuity_score']],results[k][hd['Mean_enrichment']] ]
                    peptides_composition_dic[tot]=self.C.composition_dic[k]
                del results
        if self.logfile is not None : 
            self.logfile.write(' %d strand matching condition imposed by %s\n' % (tot_matching,str(target_strand)))
        return peptides_satisfying_condition,peptides_composition_dic
        
        
        
        
        

#compare two float numbers at a specified sensibliity
def CompareFloats(float1, float2, sensibility=0.0001) :
    if float1-sensibility <= float2 <= float1+sensibility :
        return True
    else :
        return False  
    
    
    
    
    
    
    
    
    
# when executed and not imported, just run analyze_dssp_folder() to create the BSn database:
if __name__ == '__main__':
    if sys.argv[1]=='make_database' : # old way of making database (it runs dssp on the pdb)
        sys.stdout.write("making dssp database. Hopefully you have dssp in /usr/bin otherwise kill this!\n  looking for pdb files in %s, writing dssp_files to %s\n\n" %(sys.argv[2],sys.argv[3]))
        make_dssp_database_from_pdb_folder(sys.argv[2],sys.argv[3], check_for_already_present=True, only_these_files=None)
        sys.exit(0)
    if '--all' in sys.argv[1:] :# this makes BS and BSn libraries from folder with DSSP file and pdb90 chain list. If the latter is NOT GIVEN it will aslo rsync the DSSP folder...
                                #   output folder name can also be given as output:my_folder_name (a single argv with : that separate the keyword from the name)
        list_chain=None 
        dssp_folder=None
        pickle_folder=None
        pickle_file=None
        debug=False
        pickle_file_strange=None
        output_folder='BSn_libraries/' 
        for arg in sys.argv[1:] :
            if 'list' in arg : list_chain=arg
            elif 'debug' in arg : debug=True
            elif os.path.isdir(arg) :
                file_list=os.listdir(arg)
                if '.dssp' in file_list[0] or '.dssp' in file_list[3] :
                    dssp_folder=arg
                elif '.pkl' in file_list[0] or '.pkl' in file_list[1] :
                    pickle_folder=arg
            elif '.pkl' in arg :
                if 'strange' in arg :
                    if pickle_file_strange is None : pickle_file_strange=arg
                    else :
                        sys.stderr.write("WARNING pickle_file_strange = %s already given..\n" % pickle_file_strange)
                else :
                    if pickle_file is None : pickle_file=arg
                    else :
                        sys.stderr.write("WARNING pickle_file = %s already given..\n" % pickle_file)
            elif 'output:' in arg :
                output_folder=arg.replace('output:','').replace('"','').replace("'","")
        sys.stdout.write("\ninput: list_chain=%s dssp_folder=%s output_folder=%s pickle_file=%s pickle_file_strange=%s\n" % (str(list_chain),str(dssp_folder),str(output_folder),str(pickle_file),str(pickle_file_strange) ) )
        if debug : sys.stdout.write('using debug\n')
        sys.stdout.flush()
        if dssp_folder is not None :
            if list_chain is None : chain_list_file_name,_=update_dssp_database(dssp_folder, update_repr_chain_list=True, chain_list_identity=90,chain_list_name=list_chain, overwrite_dssp=False)
            else : chain_list_file_name=list_chain
            strand_dict, strange_strand_dict, get_aa_abundance = analyze_dssp_folder(dssp_folder,chain_list_file=chain_list_file_name,get_aa_abundance=True,debug=debug)
            generate_frag_dic_folder(strand_dict, length_range=list(range(1,9)), only_actual_n_plets=False, strange_strand_dic_or_pickle=strange_strand_dict, output_folder=output_folder )
            sys.stdout.write("OUTPUT in : \n BSn in: %s  BS in: strand_dict.pkl, strange_strand_dict.pkl \n" % (output_folder))
        elif pickle_file is not None :
            generate_frag_dic_folder(pickle_file, length_range=list(range(1,9)), only_actual_n_plets=False, strange_strand_dic_or_pickle=pickle_file_strange, output_folder=output_folder )
            sys.stdout.write("OUTPUT in : \n %s \n" % (output_folder))
        sys.exit(0) 
            
    if len(sys.argv)==2 :
        if '.pkl' in sys.argv[1] :
            fragments_from_database(sys.argv[1])
        else :
            analyze_dssp_folder(sys.argv[1])
    elif len(sys.argv)==3 :
        analyze_dssp_folder(sys.argv[1],chain_list_file=sys.argv[2],get_aa_abundance=True)
    else :
        sys.stderr.write("\nUSAGE: python beta_strand.py folder_with.dssp_files [optional: chain_list_file that contains 2 column, pdb_id and semicolon-separated chain ids to consider ]\n\n") 
        sys.exit(1)
    
    
    
