 # -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 18:46:29 2020

@author: elbangue
"""
import os

import sys
import pandas as pd
import glob
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import db_to_db_blaster_unisaar_2021 as db
#
#First Arm: Blast protens query agaisnst reference proteomes
#

#Second Arm organiztion for PPIFT
#
#   First Part: Generation of single table with blast results 
#



print('\n                   ========================================================================')
print('| Second_arm_strike       | Preparing for PPIFT analysis | Please wait! | All is good so far! |')
print('                      =========================================================================')


#Progreso rutinas to screen
def progres(count, total, status='', bar_len=60):
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    fmt = '[%s] %s%s ...%s' % (bar, percents, '%', status)
    #print('\b' * len(fmt), end='\r')  # clears the line
    sys.stdout.write(fmt)
    sys.stdout.flush()
    

def mafft_trimal (var_mi, var_out): # is only mafft, trimAl was removed from here due to Alphfold option
   #call mafft
    all_together = str ('mafft --auto '+ var_mi +'>'+var_out)
    os.system(all_together)
   #call trimal 
   #call_trimal ='/home/elbangue/Documents/Bioinfo/Programs/trimAl/source/trimal'
   #all_together_trimal_1 = str (call_trimal+ ' -in ' + var_out +' -out ' + var_out_trimal_1 +' -automated1' )
   #os.system(all_together_trimal_1)
  
   #delete gaps from alignments
    #with open(var_out_trimal_2, "w") as o:
     #   for record in AlignIO.read(var_out_trimal_1, "fasta"):
      #      record.seq = record.seq.ungap("-")
       #     SeqIO.write(record, o, "fasta") 
    
   #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 25 16:13:29 2021

@author: elbangue
"""
  
first_row = ''     
internal_path =db.nextscripts
all_data = pd.DataFrame([])
x = 0
text_files = [f for f in sorted(os.listdir(internal_path)) if f.endswith('.txt')]

os.chdir(internal_path+'/')  
print ('Init..')  
for n in text_files:
    temp= pd.read_csv(n, sep='\t', header=[3], engine = 'python')
    all_data.loc[:,x] = temp.iloc [:,3]
    progres (x, len (text_files), status = '#')
    x = x+1
if not os.path.exists(internal_path + r'/Results_BRUT'):    
    os.mkdir(internal_path+r'/Results_BRUT')     
all_data.to_csv (internal_path+r'/Results_BRUT/blast_BRUT.csv', index = False, header = False)    


print ('keep going..')
os.chdir(internal_path+r'/Results_BRUT')

with open(internal_path+r'/Results_BRUT/blast_BRUT.csv', 'r') as file:
    data = file.read().split('\n')
    
for row in range(0, len(data)):
    progres (row, len (data), status = '#')
    third_col= data[row].split(',')
    for x in third_col:
      #  A.thaliana, must be in the analysis  
        if x.find ('ARATH')!=-1: #CORCC   ARATH YEAST ECO HUMAN
            first_row = x
            #A0A2T2NIJ9_CORCC
           
    i = first_row [first_row.rfind('|')+1: len (first_row)] 
    if 'ARATH' in first_row: #CORCC ARATH YEAST ECO HUMAN
        with open(i + '.txt', 'w') as output: #str(row)
            
            for x in third_col:
                
                os.chdir(internal_path+r'/Results_BRUT')
                output.write(x + os.linesep)
file.close()

# Second Arm organiztion for PPIFT
#
# Second Part: Generation of files to feed PPIFT
#

selected_seqs = list()

internal_path_genomes =db.mypath_genomes

# Initial adjustment. Need to made a big fasta file with all proteomes

print ('still working on second arm, please wait...') 

os.chdir(internal_path_genomes)

oh = open('one_fasta_file.fasta', 'w')
for f in db.floating_var:#os.listdir(internal_path_genomes):
    
    if f.endswith('.fasta'):
        fh = open(os.path.join(internal_path_genomes, f))
        for line in fh:
            oh.write(line)
        fh.close()
     
oh.close()


if not os.path.exists(internal_path+r'/Results_BRUT/Results_II'):
    os.makedirs(internal_path+r'/Results_BRUT/Results_II')
 
x=0
for filename in glob.glob(os.path.join(internal_path+r'/Results_BRUT/', '*.txt')):
   totalfiles =   len([name for name in sorted(os.listdir('.')) if os.path.isfile(name)])
   
   with open(os.path.join(os.getcwd(), filename), 'r') as f: 
       progres(x, totalfiles, status='#', bar_len=60) 
       accessionids = [line.rstrip("\n") for line in f]
       selected_seqs = list()
       
   for seq_record in SeqIO.parse(internal_path_genomes+'/'+'one_fasta_file.fasta', "fasta"):
          header = seq_record.name # (or .id or so)
                            
          for accession_id in accessionids:
              if accession_id == header:
                  selected_seqs.append(seq_record)
                  os.chdir (internal_path+r'/Results_BRUT/Results_II')
                   
                  SeqIO.write(selected_seqs, filename[filename.rfind('/')+1:filename.rfind('.')]+".fasta", "fasta")
   x=x+1  
       
#align_with mafft,  edit alignment with trimal and delete '-' gap character from sequences for next analysis

if not os.path.exists(internal_path+r'/Results_BRUT/Results_III'):
    os.makedirs(internal_path+r'/Results_BRUT/Results_III')
if not os.path.exists(internal_path+r'/Results_BRUT/Results_III_notAl'):
    os.makedirs(internal_path+r'/Results_BRUT/Results_III_notAl')
if not os.path.exists(internal_path+r'/Results_BRUT/Results_IV'):
    os.makedirs(internal_path+r'/Results_BRUT/Results_IV')
if not os.path.exists(internal_path+r'/Results_BRUT/Results_V'):
    os.makedirs(internal_path+r'/Results_BRUT/Results_V')

     
#Mafft align, trimAl edition and asa_Alphafold2 selection
  
for filename in glob.glob(os.path.join(internal_path+r'/Results_BRUT/Results_II', '*.fasta')):
   totalfiles =   len([name for name in sorted(os.listdir('.')) if os.path.isfile(name)])
 
   with open(os.path.join(os.getcwd(), filename), 'r') as f: 
       progres(x, totalfiles, status='#', bar_len=60) 
       accessionids = [line.rstrip("\n") for line in f]
       selected_seqs = list()
   var_out_p  = internal_path+r'/Results_BRUT/Results_III'+filename[filename.rfind('/'):filename.rfind('.')]+".msa"
   #var_out_trimal_1_p  = internal_path+r'/Results_BRUT/Results_III'+filename[filename.rfind('/'):filename.rfind('.')]+".fasta"
   mafft_trimal(filename, var_out_p)
   
   
#Comienzo asa aminoacidos!!!!!!!!!! Gracias Aplphafold2!!!!!!!!!
 
#función para extraer nombre con extension de ficheros en directorio
def files(path):  
    for file in sorted(os.listdir(path)):
        if os.path.isfile(os.path.join(path, file)):
            yield file
            
#abrir fichero con los id y valores de exposicion a solvente per aminoacido
#Para YEAST
#file1 = open(r'/home/elbangue/Downloads/UP000002311_559292_YEAST/SASA/new_sequences_scerevisiae_goodheader.sa', 'r') 

#Para Arabidospsis
#file1 = open(r'/home/elbangue/Downloads/UP000006548_3702_ARATH/SASA/athaliana_refineassit_3dfrom_alphafold2.as', 'r')

# Lines = file1.readlines()
asa_line = str()
var_out_trimal_1_all_files_fasta = internal_path+r'/Results_BRUT/Results_III/*.msa'
var_out_asa_all_files = internal_path+ r'/Results_BRUT/Results_IV/'
var_out_asa_all_files_notAl_first = internal_path+ r'/Results_BRUT/Results_III/'
var_out_trimal_2_all_files = internal_path+r'/Results_BRUT/Results_V/'
var_out_asa_all_files_notAl = internal_path+ r'/Results_BRUT/Results_III_notAl/'#para no procesar alphafold se analiza el folder III directamente : r'/Results_BRUT/Results_III_notAl'
otro_contador = 0



#Esta seccion debe moverse de lugar y es para alphafold

#abrir fichero con los id y valores de exposicion a solvente per aminoacido
#Para YEAST
#file1 = open(r'/home/elbangue/Downloads/UP000002311_559292_YEAST/SASA/new_sequences_scerevisiae_goodheader.sa', 'r') 

#Para Arabidospsis
#file1 = open(r'/home/elbangue/Downloads/UP000006548_3702_ARATH/SASA/athaliana_refineassit_3dfrom_alphafold2.as', 'r')

#Lines = file1.readlines()
#asa_line = str()
#recorrer carpeta con alineamientos mafft  
#for filename in glob.glob(var_out_trimal_1_all_files_fasta):
 #buscar en fichero con exposición a solvente el id correspodiente y la secuencia correspondiente de la exposicion de los aminoacidos al solvente
 #   aux_file_var = str (filename [filename.rfind('/')+1:filename.rfind('_')])#cambie . por _
  #  var_out_asa_p  = internal_path+r'/Results_BRUT/Results_IV/'+aux_file_var+".fasta"
   # var_out_trimal_2_p =internal_path+r'/Results_BRUT/Results_V/'+aux_file_var+".fasta"
   # count = 0
    # Strips the newline character
    #print(filename)
    #print(filename [filename.rfind('/')+1:filename.rfind('.')]
    
    #for line in Lines:
     #   if (line.find (aux_file_var) !=-1) and (count < (len(Lines))):
      #      asa_line = Lines[count +1] 
       #     print (line)
        #    print ('_------------_')
         #   print (count)
          #  print (asa_line)
            
           # len (Lines)
           # break
       # count = count + 1 
   
 #buscar alineamiento correspondiente mediante id del fichero, se debe incluir chequeo de secuencia   
   # align = AlignIO.read(filename, "fasta")
    
    #for record in align:
     #   aux =str(record.id)
      #  if (aux.find(aux_file_var) != -1):
       #     aux2 = record.id
 #moverse por secuencia de alineamiento correspondiente y colocar gap donde los aminoacidos no esten expuestos a solvente y adiciona secuencia guía
            
        #    position = 0            
         #   my_aux_fake_prot = str()
          #  for aminoacid in record.seq:
           #     if position <len(asa_line): 
                    
            #        if  asa_line [position] == 'B' or asa_line [position] == 'C':#asa_line [position] == 'B' or asa_line [position] == 'C': B-buried S-exposed C-core
             #           my_aux_fake_prot = my_aux_fake_prot + '-'
                        
                        
                        #print (aminoacid,record.seq[position])
                        #print (Lines[count])
                        #print (aux2)
              #      else:
               #         my_aux_fake_prot = str (my_aux_fake_prot + aminoacid)
                        
                #    if aminoacid !='-':
                 #       position = position +1
            #remiendo por problemas que no entiendo en algunos resultados finales  
           # while len (my_aux_fake_prot)<len (record.seq):
            #    my_aux_fake_prot = my_aux_fake_prot + '-'
    #seq =Seq(my_aux_fake_prot)
    #dummy = SeqRecord(seq, id=aux2)
    #print(len(dummy.seq), len (record))
    #align.append (dummy)
    #dummy2 = []
    #print ('aminoacids relative to alignment and Alphafold2 3D structure prediction: ', end = '\r')
    #progres(otro_contador, totalfiles/2, status='#', bar_len=60)
#remove similar sequences
     
    #for var1 in align:
     #   booleana = False  
        #print (var1)
      #  i =0 
       # dummy3 = []
        #for var2 in  align:
         #  if var1.id == var2.id: 
          #     i =i +1
           #    if i==2:
            #       booleana = True
             #      dummy3 = SeqRecord (var1.seq, var1.id, str(len (var1.seq)))
                   #break
        #if booleana == False: 
            #print('%%%%')            
         #   dummy2.append(var1)
    #dummy2.append(dummy3)
    #fd = os.open( var_out_asa_p, os.O_RDWR|os.O_CREAT)
    #SeqIO.write(dummy2, var_out_asa_p, "fasta") 
    ##fd =os.open (var_out_asa_p, os.O_WRONLY)
    #os.close (fd)
    #otro_contador = otro_contador +1
    

#trimal for final gap removoal 
y = 0
alphafold = False
if alphafold:
    
    #Activate for Alphafold, desactivar en caso contrario:
    for file in files(var_out_asa_all_files):  
    
#Activate for NOT Alphafold use, desactivar en caso contrario:
#for file in files(var_out_asa_all_files_notAl_first):    
    
        call_trimal =r'/home/elbangue/Documents/Bioinfo/Programs/trimAl/source/trimal'
    #var_out_trimal_1 = 'filename'
    #var_out_trimal_2 =
    #all_together_trimal_2 = str (call_trimal+ ' -in ' +internal_path+r'/Results_BRUT/Results_IV/' + file + ' -out ' + var_out_asa_all_files + file [file.rfind('/')+1:file.rfind('.')]+'.fasta' +' -automated1' ) 
   
#activar para procesar datos con Alphafold, desactivar en caso contrario:
        all_together_trimal_3 = str (call_trimal+ ' -in ' +internal_path+r'/Results_BRUT/Results_IV/' + file + ' -out ' + var_out_trimal_2_all_files+file [file.rfind('/')+1:file.rfind('.')]+'.fasta' +' -nogaps' )
    
#activar para NO procesar datos de Alphafold desactivar en caso contrario. Edit after Mafft alignment with trimAl -automated1 and -nogaps, also include call to trimAl: 
    #all_together_trimal_2 = str (call_trimal +' -in ' +internal_path+r'/Results_BRUT/Results_III/' + file + ' -out ' + var_out_asa_all_files_notAl +file [file.rfind('/')+1:file.rfind('.')]+'.fasta' +' -automated1' )
    #all_together_trimal_3 = str (call_trimal +' -in ' +internal_path+ r'/Results_BRUT/Results_III_notAl/' + file [file.rfind('/')+1:file.rfind('.')]+'.fasta' + ' -out ' + var_out_trimal_2_all_files+file [file.rfind('/')+1:file.rfind('.')]+'.fasta' +' -nogaps' )
    #os.system(all_together_trimal_2)
    
        os.system(all_together_trimal_3)
    #s = str(file) + '%'                        # string for output
    #print (s, end ='\r')                # just print and flush
        sys.stdout.write("trim_asa_progress: ")
        sys.stdout.flush()
    #time.sleep(0.05)
        progres(y, totalfiles/2, status='#', bar_len=60)
        y = y+1
##fin de asa aminoacidos!!!!!!!!!! Gracias Aplphafold2!!!!!!!!!

else:
    #Activate for Alphafold, desactivar en caso contrario:
    #for file in files(var_out_asa_all_files):  
    
#Activate for NOT Alphafold use, desactivar en caso contrario:
    for file in files(var_out_asa_all_files_notAl_first):    
    
        call_trimal =r'/home/elbangue/Documents/Bioinfo/Programs/trimAl/source/trimal'
        #var_out_trimal_1 = 'filename'
        #var_out_trimal_2 =
        #all_together_trimal_2 = str (call_trimal+ ' -in ' +internal_path+r'/Results_BRUT/Results_IV/' + file + ' -out ' + var_out_asa_all_files + file [file.rfind('/')+1:file.rfind('.')]+'.fasta' +' -automated1' ) 
   
#activar para procesar datos con Alphafold, desactivar en caso contrario:
        #all_together_trimal_3 = str (call_trimal+ ' -in ' +internal_path+r'/Results_BRUT/Results_IV/' + file + ' -out ' + var_out_trimal_2_all_files+file [file.rfind('/')+1:file.rfind('.')]+'.fasta' +' -nogaps' )
    
#activar para NO procesar datos de Alphafold desactivar en caso contrario. Edit after Mafft alignment with trimAl -automated1 and -nogaps, also include call to trimAl: 
        all_together_trimal_2 = str (call_trimal +' -in ' +internal_path+r'/Results_BRUT/Results_III/' + file + ' -out ' + var_out_asa_all_files_notAl +file [file.rfind('/')+1:file.rfind('.')]+'.fasta' +' -automated1' )
        all_together_trimal_3 = str (call_trimal +' -in ' +internal_path+ r'/Results_BRUT/Results_III_notAl/' + file [file.rfind('/')+1:file.rfind('.')]+'.fasta' + ' -out ' + var_out_trimal_2_all_files+file [file.rfind('/')+1:file.rfind('.')]+'.fasta' +' -nogaps' )
        os.system(all_together_trimal_2)
    
        os.system(all_together_trimal_3)
    #s = str(file) + '%'                        # string for output
    #print (s, end ='\r')                # just print and flush
        sys.stdout.write("trim_asa_progress: ")
        sys.stdout.flush()
    #time.sleep(0.05)
        progres(y, totalfiles/2, status='#', bar_len=60)
        y = y+1
##fin de asa aminoacidos!!!!!!!!!! Gracias Aplphafold2!!!!!!!!!
  
   

os.chdir(internal_path+r'/Results_BRUT/Results_V') 
eraser = []
for name in sorted(os.listdir('.')):   
   fd = os.open( name, os.O_RDWR|os.O_CREAT)
   booleana = False  
   for seq_record in SeqIO.parse (name, "fasta"):
          
       for seq_record_2 in  SeqIO.parse (name, "fasta"):
           if len (seq_record.seq)>(len(seq_record_2.seq)*2)  or len (seq_record.seq)==0: 
               booleana = True
   if booleana == True: 
       eraser.append(name)
   os.close (fd)
print (len (eraser)+1, 'proteins were taken out of the analysis due to sequences problems..')    
for real_eraser in eraser:
   print (real_eraser.partition('.')[0], end=",")
   os.remove (real_eraser)
print ('--------------------------->')
print ('Second arm strike completed')

os.system ('ntfy -b telegram send second_arm_PlantInt_finished')
#with open("out.fa", "w") as o:
  #  for record in AlignIO.read("test.fa", "fasta"):
   #     record.seq = record.seq.ungap("-")
    #    SeqIO.write(record, o, "fasta")

