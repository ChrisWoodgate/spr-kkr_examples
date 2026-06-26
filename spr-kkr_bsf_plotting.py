'''
sprkkr_bsf_plotting.py

Routine(s) for handling SPR-KKR BSF output for pretty plotting in
Python, e.g. using matplotlib.

Original code written by Dr. Aki Pulkkinen, University of West Bohemia

C. D. Woodgate, University of Bristol, 2026
'''

import sprkkr_dos_plotting
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.constants as c
from matplotlib import rc
from matplotlib import cm
import matplotlib.ticker as ticker
import matplotlib.colors as mcolors
import sys
from scipy.interpolate import griddata

Ry = 0.5*c.value('Hartree energy in eV')

class BSFdata():
    def __init__(self):
        self.KEYWORD = ''
        self.NQeff = 0
        self.NTeff = 0
        self.EFERMI = 0
        self.mode = ''
        self.NE = 0
        self.NK = 0
        self.NK1 = 0            
        self.NK2 = 0    
        self.VECK1 = [0,0,0]
        self.VECK2 = [0,0,0]
        self.ERYD = [0,0]
        self.EMIN = 0
        self.EMAX = 0
        self.NKDIR = 0
        self.INDKDIR = []   
        self.IT = []
        self.k = [] 
        self.E = []
        self.k1 = []
        self.k2 = []
        self.Iup = []
        self.Idn = []     
        self.Ix = []
        self.Iy = []
        self.Iz = []
        self.I = []   



def parse_bsfdata(filename):
    data = BSFdata()

    with open(filename,'r') as f:
        line = f.readline()
        data.KEYWORD = line.split()[1]
        
        
        if data.KEYWORD == 'BSF':
            while line[0:6] != 'NQ_eff':
                line = f.readline()
            data.NQeff = int(line.split()[1])
            line = f.readline()
            data.NTeff = int(line.split()[1])  
            line = f.readline()
            line = f.readline()
            line = f.readline()
            line = f.readline()
            data.EFERMI = float(line.split()[1])  

            while line[0:5] != '   IT':
                line = f.readline()   
            for i in range(data.NTeff):    
                line = f.readline()
                IQAT = int(line.split()[3])
                if(IQAT>10): # 2 lines. Need to genearlize for more?
                    for j in range(10):
                        index = line.split()[1]
                        data.IT.insert(int(line.split()[4+j]),index)
                    line = f.readline()
                    for j in range(IQAT-10):
                        data.IT.insert(int(line.split()[j]),index)
                else:
                    for j in range(IQAT):
                        data.IT.insert(int(line.split()[4+j]),line.split()[1])

            line = f.readline()            
            line = f.readline()
            line = f.readline()        
            data.mode = line.split()[1]
        
###############
#   k-E map   #
###############        
            if data.mode == 'EK-REL':
                line = f.readline()
                data.NE = int(line.split()[1])
                line = f.readline() 
                data.NK = int(line.split()[1])        
                line = f.readline()    
                line = f.readline() 
                data.EMIN = float(line.split()[2])            
                data.EMAX = float(line.split()[3])    
                line = f.readline()            
                line = f.readline()            
                data.NKDIR = int(line.split()[1])
            
                while line[0:5] != '#####':
                    line = f.readline()      
            
                for i in range(data.NKDIR):           
                    line = f.readline()            
                    data.INDKDIR.append(int(line.split()[1]))
            
                line = f.readline()             
                line = f.readline()   
                     
                for i in range(data.NK):           
                    line = f.readline()            
                    data.k.append(float(line))   
                                             
                line = f.readline()  
                data.I = np.zeros((data.NQeff,data.NE,data.NK))        
                data.Iup = np.zeros((data.NQeff,data.NE,data.NK))            
                data.Idn = np.zeros((data.NQeff,data.NE,data.NK))            

                for i in range(data.NE):
                    for j in range(data.NQeff):
                        for k in range(data.NK):
                            line = f.readline()
                            data.Iup[j,i,k] = float(line)
                                                   
                    for j in range(data.NQeff):
                        for k in range(data.NK):
                            line = f.readline()
                            data.Idn[j,i,k] = float(line)	                
                        
                for i in range(data.NE):                        
                    for j in range(data.NQeff):
                        for k in range(data.NK):
                            line = f.readline()
                            data.I[j,i,k] = float(line)	      
          
                data.k,data.E = np.meshgrid(np.linspace(0,np.max(data.k),data.NK),np.linspace(data.EMIN,data.EMAX,data.NE))
          
###########################
#   constant energy map   #
########################### 
            if data.mode == 'CONST-E':
                line = f.readline()
                data.NK1 = int(line.split()[1])
                line = f.readline()            
                data.NK2 = int(line.split()[1])        
                line = f.readline()   
                data.ERYD[0] = float(line.split()[1])
                data.ERYD[1] = float(line.split()[2])            
                line = f.readline() 
                line = f.readline() 
                line = f.readline() 
                data.VECK1[0] = float(line.split()[1])
                data.VECK1[1] = float(line.split()[2])
                data.VECK1[2] = float(line.split()[3])
                line = f.readline() 
                line = f.readline()             
                data.VECK2[0] = float(line.split()[1])
                data.VECK2[1] = float(line.split()[2])
                data.VECK2[2] = float(line.split()[3])
                line = f.readline()             
            
            
                data.I = np.zeros((data.NQeff,data.NK1,data.NK2))
                data.Iup = np.zeros((data.NQeff,data.NK1,data.NK2))            
                data.Idn = np.zeros((data.NQeff,data.NK1,data.NK2))            

                for i in range(data.NK1):
                    for j in range(data.NQeff):
                        for k in range(data.NK2):
                            line = f.readline()
                            data.Iup[j,i,k] = float(line)
                                               
                    for j in range(data.NQeff):
                        for k in range(data.NK2):
                            line = f.readline()
                            data.Idn[j,i,k] = float(line)	                
                        
                for i in range(data.NK1):                        
                    for j in range(data.NQeff):
                        for k in range(data.NK2):
                            line = f.readline()
                            data.I[j,i,k] = float(line)	               
          
# Here the K1 and K2 vectors are constructed from zero to the length of the vector. They should start from the coordinates of the vector KA specified in the input file, but these coordinates are not available in the .bsf file.
                #data.k1,data.k2 = np.meshgrid(np.linspace(0,data.VECK2[1],data.NK2),np.linspace(0,data.VECK1[0],data.NK1))
                data.k1,data.k2 = np.meshgrid(np.linspace(0,2,data.NK2),np.linspace(0,1,data.NK1))
                    
        elif data.KEYWORD == 'BSF-SPOL':
                   
            while line[0:6] != 'NQ_eff':
                line = f.readline()
            data.NQeff = int(line.split()[1])
            line = f.readline()
            data.NTeff = int(line.split()[1])  
            line = f.readline()
            line = f.readline()
            line = f.readline()
            line = f.readline()
            data.EFERMI = float(line.split()[1])  

            while line[0:5] != '   IT':
                line = f.readline()   
            for i in range(data.NTeff):    
                line = f.readline()
                IQAT = int(line.split()[3])
                if(IQAT>10): # 2 lines. Need to genearlize for more?
                    for j in range(10):
                        index = line.split()[1]
                        data.IT.insert(int(line.split()[4+j]),index)
                    line = f.readline()
                    for j in range(IQAT-10):
                        data.IT.insert(int(line.split()[j]),index)
                else:
                    for j in range(IQAT):
                        data.IT.insert(int(line.split()[4+j]),line.split()[1])


            line = f.readline()            
            line = f.readline()
            line = f.readline()        
            data.mode = line.split()[1]
        
###############
#   k-E map   #
###############        
            if data.mode == 'EK-REL':
                line = f.readline()
                data.NE = int(line.split()[1])
                line = f.readline() 
                data.NK = int(line.split()[1])        
                line = f.readline()    
                line = f.readline() 
                data.EMIN = float(line.split()[2])            
                data.EMAX = float(line.split()[3])    
                line = f.readline()            
                line = f.readline()            
                data.NKDIR = int(line.split()[1])
            
                while line[0:5] != '#####':
                    line = f.readline()      
            
                for i in range(data.NKDIR):           
                    line = f.readline()            
                    data.INDKDIR.append(int(line.split()[1]))
            
                line = f.readline()             
                line = f.readline()   
                     
                for i in range(data.NK):           
                    line = f.readline()            
                    data.k.append(float(line))   
                                             
                line = f.readline()  
                data.I = np.zeros((data.NQeff,data.NE,data.NK))        
            
                data.Ix = np.zeros((data.NQeff,data.NE,data.NK))            
                data.Iy = np.zeros((data.NQeff,data.NE,data.NK))            
                data.Iz = np.zeros((data.NQeff,data.NE,data.NK))
            
                for i in range(data.NE):                        
                    for j in range(data.NQeff):
                        for k in range(data.NK):
                            line = f.readline()
                            data.I[j,i,k] = float(line)            

                    for j in range(data.NQeff):
                        for k in range(data.NK):
                            line = f.readline()
                            data.Ix[j,i,k] = float(line)
                                               
                    for j in range(data.NQeff):
                        for k in range(data.NK):
                            line = f.readline()
                            data.Iy[j,i,k] = float(line)	                

                    for j in range(data.NQeff):
                        for k in range(data.NK):
                            line = f.readline()
                            data.Iz[j,i,k] = float(line)	
                        
	      
          
                data.k,data.E = np.meshgrid(np.linspace(0,np.max(data.k),data.NK),np.linspace(data.EMIN,data.EMAX,data.NE))
          
###########################
#   constant energy map   #
########################### 
            if data.mode == 'CONST-E':
                line = f.readline()
                data.NK1 = int(line.split()[1])
                line = f.readline()            
                data.NK2 = int(line.split()[1])        
                line = f.readline()   
                data.ERYD[0] = float(line.split()[1])
                data.ERYD[1] = float(line.split()[2])            
                line = f.readline() 
                line = f.readline() 
                line = f.readline() 
                data.VECK1[0] = float(line.split()[1])
                data.VECK1[1] = float(line.split()[2])
                data.VECK1[2] = float(line.split()[3])
                line = f.readline() 
                line = f.readline()             
                data.VECK2[0] = float(line.split()[1])
                data.VECK2[1] = float(line.split()[2])
                data.VECK2[2] = float(line.split()[3])
                line = f.readline()             
            
            
                data.I = np.zeros((data.NQeff,data.NK1,data.NK2))
                data.Ix = np.zeros((data.NQeff,data.NK1,data.NK2))            
                data.Iy = np.zeros((data.NQeff,data.NK1,data.NK2))            
                data.Iz = np.zeros((data.NQeff,data.NK1,data.NK2))  
            
            
                for i in range(data.NK1):      
                    for j in range(data.NQeff):
                        for k in range(data.NK2):
                            line = f.readline()
                            data.I[j,i,k] = float(line)
                for i in range(data.NK1):            
                    for j in range(data.NQeff):
                        for k in range(data.NK2):
                            line = f.readline()
                            data.Ix[j,i,k] = float(line)
                for i in range(data.NK1):                                               
                    for j in range(data.NQeff):
                        for k in range(data.NK2):
                            line = f.readline()
                            data.Iy[j,i,k] = float(line)	                
                for i in range(data.NK1):
                    for j in range(data.NQeff):
                        for k in range(data.NK2):
                            line = f.readline()
                            data.Iz[j,i,k] = float(line)	
                                               
	               
# Same problem with K1 and K2 as in the .bsf file parser (see above).          
                data.k1,data.k2 = np.meshgrid(np.linspace(0,np.linalg.norm(data.VECK2),data.NK2),np.linspace(0,np.linalg.norm(data.VECK1),data.NK1))
          
        else:
            print("Datafile is not a valid BSF or BSF-SPOL file.")
            sys.exit(1)  
          
                    
    return data
