#----------------------------------------------------------------------
#script to run gravitional potential field heatmap, analytic and simulational methods
#----------------------------------------------------------------------
#usage: python3 runmap.py
#----------------------------------------------------------------------
#                         0.3 version
#                           10/2021
#Bruno Morini
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#packages  
#----------------------------------------------------------------------
import os #terminal 
import numpy as npy #math solver
import PyGnuplot as gp #GNUplot module package

#----------------------------------------------------------------------
#constants
#----------------------------------------------------------------------
global G
G = 6.67408E-11 #gravitataional constant

global k
k = 8.98755E9 #Coulomb constant

global pi
pi = 3.141592653589793238 #pi 









#----------------------------------------------------------------------
#read data.info file information  
#----------------------------------------------------------------------
print('Current data file:\n')

with open('/home/bruno/ic/simulations/gmap/lmp/data/data.info', 'r') as rdatainfo: #open for read lammps data.info file

  text = rdatainfo.read()
  print(text+'\n')

rdatainfo.close() #close file

with open('/home/bruno/ic/simulations/gmap/lmp/data/data.info', 'r') as rdatainfo: #open for read lammps data.info file

  for i in range(0,1): #read file line 1 
    global a
    a = rdatainfo.readline() #save line text
    a = a.replace('Semi-axis x: ', '') #save semi-axis x value 
    a = float(a)
  
  for i in range(0,1): #read file line 2
    global b
    b = rdatainfo.readline() #save line text
    b = b.replace('Semi-axis y: ', '') #save semi-axis y value
    b = float(b)

  for i in range(0,1): #read file line 3 
    global c
    c = rdatainfo.readline() #save line text
    c = c.replace('Semi-axis z: ', '') #save semi-axis z value
    c = float(c)
  
  for i in range(0,1): #read file line 4
    global cnt
    cnt = rdatainfo.readline() #save line text
    cnt = cnt.replace('Number of particles: ', '') #save number of particles value
    cnt = float(cnt)
    cnt = int(cnt)
  
  for i in range(0,1): #read line file 5 
    global M
    M = rdatainfo.readline() #save line text
    M = M.replace('Body mass: ', '') #save body mass value 
    M = float(M)
    m = M/cnt #mass of individual particles
    q = m*(npy.sqrt(G/k)) #Converting mass to charge 

  
  for i in range(0,1): #read line file 6 
    global de
    de = rdatainfo.readline() #save line text
    de = de.replace('Density: ', '') #save body mass value 
    de = float(de)

rdatainfo.close() #close file










#----------------------------------------------------------------------
#Select or create data file
#----------------------------------------------------------------------
def selectdata(): #function to select or create data file

  global newdata #global variable
  newdata = str(input('New data file ? (y/n)\n')) 

  if (newdata == 'n'): #test variable text
    return #end function

  elif (newdata == 'y'): #test variable text
    global m
    m = float(input('Particle masses (q=m): '))
    global q
    q = m*(npy.sqrt(G/k)) #Converting mass to charge 
    global d
    d = float(input('Particle Diameter: '))
    global de
    de = m/((4.0/3.0)*pi*((d/2.0)**3.0)) #particle density 

    def semi_axis_a(): #function to verify semi-axis x > and multiple particle diameter)
      global a
      a = float(input('Semi-axis x (> and multiple of Particle Diameter): ')) 
  
      if (npy.abs(a/d).is_integer()): #verify if a is bigger and multiple than d
        return
      else:
        print('Error: semi-axis x is not bigger and multiple of the particle diameter')
        semi_axis_a() #restart function

    semi_axis_a() #start function     

    def semi_axis_b(): #verify semi-axis y > and multiple particle diameter)
      global b
      b = float(input('Semi-axis y (> and multiple of Particle Diameter): ')) 
  
      if (npy.abs(b/d).is_integer()): #verify if b is bigger and multiple than d
        return
      else:
        print('Error: semi-axis y is not bigger and multiple of the particle diameter')
        semi_axis_b() #restart function

    semi_axis_b() #start function   

    def semi_axis_c(): #verify semi-axis c > and multiple particle diameter)
      global c
      c = float(input('Semi-axis z (> and multiple of Particle Diameter): ')) 
  
      if (npy.abs(c/d).is_integer()): #verify if c is bigger and multiple than d
        return  
      else:
        print('Error: semi-axis z is not bigger and multiple of the particle diameter')
        semi_axis_c() #restart function

    semi_axis_c() #start function     

    #------------------------------------------------------------------
    #golden ratio and angle
    #------------------------------------------------------------------
    gr=(npy.sqrt(5.0)+1.0)/2.0 #golden ratio
    ga=(2.0-gr)*(2.0*3.1415) #golden angle

    #------------------------------------------------------------------
    #set major segment of ellipsoid
    #------------------------------------------------------------------
    if ((a > b) and (a > c)): #a as major segment
      r = a #sphere radius
      nl = a/d #number of layers
    elif ((b > a) and (b > c)): #b as major segment
      r = b #sphere radius
      nl = b/d #number of layers
    elif ((c > a) and (c > b)): #c as major segment
      r = c #sphere radius
      nl = c/d #number of layers  
    else:
      r = a #sphere radius  
      nl = a/d #number of layers 

    #------------------------------------------------------------------
    #set particles position for each layer 
    #------------------------------------------------------------------
    global cnt
    cnt = 1
    pinfo = str()
    velocities = str()
    xyz = str()


    for i in range(1,int(nl)): #for number of particle each layer  
      nps = ((4*pi*((d*i)**2))/(pi*(d**2))) #number of particles in each layer 
  
      for j in range(1,int(nps)): #for position of particles
        lat = npy.arcsin(-1+2.0*j/(nps+1)) #latitude particle position
        lon = ga*j #longitude particle position

        x = d*i*npy.cos(lon)*npy.cos(lat) #x particle coordenates 
        y = d*i*npy.sin(lon)*npy.cos(lat) #y particle coordenates 
        z = d*i*npy.sin(lat) #z particle coordenates 

        if ((((x**2)/(a**2))+((y**2)/(b**2))+((z**2)/(c**2)))<1): #cut sphere to ellispdoid format
          cnt = cnt+1 #update counter
          
          pinfo = str(pinfo)+str(int(cnt))+' '+'1'+' '+str(x)+' '+str(y)+' '+str(z)+' '+str(d)+' '+str(de)+' '+str(q)+' '+'0 0 0'+'\n' #write lammps data

          velocities = str(velocities)+str(int(cnt))+' 0.000000000 0.000000000 0.000000000 0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00'+'\n' #write lammps data

          xyz = str(xyz)+'c '+str(x)+' '+str(y)+' '+str(z)+'\n' #write vmd xyz

    #----------------------------------------------------------------
    #write lammps data file
    #----------------------------------------------------------------
    preconfig = '#LAMMPS data file via lattice_fibonacci_ellipsoid_generator.py, version 0.1 Oct 2020, by Bruno Morini\n' + str(int(cnt)) + ' atoms\n' + '1 atom types\n\n' + str(d*-1000) + ' ' + str(d*1000) + ' ' + 'xlo xhi\n' + str(d*-1000) + ' ' + str(d*1000) + ' ' + 'ylo yhi\n' + str(d*-1000) + ' ' + str(d*1000) + ' ' + 'zlo zhi\n\n' + 'Masses\n\n' + '1' + ' ' + str(m) + '\n\n' + 'Atoms # hybrid' #origin particle

    with open('/home/bruno/ic/simulations/gmap/lmp/data/body.data', 'w+') as data: #open for read lammps output file
      data.write(str(preconfig) + '\n\n' + '1 1 0.0000000000000000 0.0000000000000000 0.0000000000000000 ' + str(d) + ' ' + str(de) + ' ' + str(q) + ' 0 0 0\n' + str(pinfo) + '\n' + 'Velocities\n\n' + '1 0.000000000 0.000000000 0.000000000 0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00\n' + str(velocities)) #write text in file

    data.close() #close file

    with open('/home/bruno/ic/simulations/gmap/lmp/data/body.xyz', 'w+') as wxyz: #open for read lammps output file
      wxyz.write(str(cnt) + '\n' + 'generated by VMD\n' + 'c 0.0 0.0 0.0\n' + str(xyz)) #write text in file
  
    wxyz.close() #close file

    #----------------------------------------------------------------
    #write data.info file
    #----------------------------------------------------------------
    datainfo = 'Semi-axis x: ' + str(a) + '\nSemi-axis y: ' + str(b) + '\nSemi-axis z: ' + str(c) + '\nNumber of particles: ' + str(cnt) + '\nBody mass: ' + str(cnt*m) + '\nDensity: ' + str(de)  #data.info text

    with open('/home/bruno/ic/simulations/gmap/lmp/data/data.info', 'w+') as wdatainfo: #open for read lammps output file
      wdatainfo.write(datainfo) #write text in file
  
    wdatainfo.close() #close file

    return #end function

  else:
    print('Invalid command!') #test variable text
    selectdata() #restart function

selectdata() #start function









#----------------------------------------------------------------------
#select type and creating gravitational map 
#--------------------------------------------------------------------- 
def selectmap(): #function to select type of gravitational map
  global typemap
  typemap = str(input('Simulational map ? (s), Analytic map ? (a), Both ? (b)\n'))

  global  npx
  npx = int(input('Map size ?\n')) #set map size

  global  pxl_den
  pxl_den = int(input('Pixels density ?\n')) #set map size


  #------------------------------------------------------------------
  #setting lammps input for simulation
  #------------------------------------------------------------------
  z = 0

  if ((typemap == 's') or (typemap == 'b')): #
    particles = str() #particle lammps input
    ctr = 1 #counter

    for x in range(-int(npx*pxl_den),int(npx*pxl_den)): #x coordenates proof particle
      x = x/pxl_den

      for y in range(-int(npx*pxl_den),int(npx*pxl_den)): #y coordenates proof particle
        y = y/pxl_den
        
        particles = str(particles)+'create_atoms       1 single '+str(x)+' '+str(y)+' '+'0.0\n'+'set                atom '+str(cnt+ctr)+' charge '+str(q)+'\n' #create first particle

        ctr = ctr+1 #counter of particles

    particles = str(particles)+'\n'+'group              gpe id '+str(cnt+1)+':'+str(cnt+ctr-1)+'\n'   

    particles = str(particles)+'neigh_modify       exclude group gpe gpe check no'

    #------------------------------------------------------------------
    #defining neighbor list and particles in lammps input file
    #----------------------------------------------------------------- 
    with open('/home/bruno/ic/simulations/gmap/lmp/in.gfield', 'r') as rinp: #open for read lammps map input file 
      with open('/home/bruno/ic/simulations/gmap/lmp/temp/in.map', 'w') as winp: #create lammps input file for proof particles
      
        for line in rinp: #search line in lammps input file
          winp.write(line.replace('datfileneighbor', 'neigh_modify       one '+str(ctr*2)+' page '+str(ctr*20)).replace('createparticles', str(particles))) #replace line in lammps input file
     
      winp.close() #close lammps particles file
    rinp.close() #close lammps input file

    #------------------------------------------------------------------
    #running simulation in lammps  
    #----------------------------------------------------------------- 
    cmd = 'lmp_gserial < /home/bruno/ic/simulations/gmap/lmp/temp/in.map' #run the simulation in lammps in serial mode
    #cmd = 'mpirun -np 4 lmp_gserial < /home/bruno/ic/simulations/gmap/lmp/in.wmap' #run the simulation in lammps parallel mode
    os.system(cmd) #execute cmd in terminal 

    #------------------------------------------------------------------
    #rewriting dump file 
    #------------------------------------------------------------------
    penergy = str() 

    with open('/home/bruno/ic/simulations/gmap/lmp/out/dump.pe', 'r') as pe: #open for read lammps output file
  
      for i in range(1,ctr+9): #search line in output lammps file
    
        if (i >= 10): #read from 9 line
          penergy = penergy+pe.readline()

        else: 
          pe.readline()

    pe.close() #close datafile

    with open('/home/bruno/ic/simulations/gmap/lmp/out/dump.pe', 'w') as pe: #open for read lammps output file
      pe.write(penergy) #write datafile

    pe.close() #close datafile
  
  global sdataout
  sdataout = str()
  global adataout
  adataout = str()
  global edataout
  edataout = str()
  
  with open('/home/bruno/ic/simulations/gmap/lmp/out/dump.pe', 'r') as pe: #open for read lammps output file     
         
    for x in range(-int(npx*pxl_den),int(npx*pxl_den)): #x coordenates proof particle
      x = x/pxl_den
      sdataout = str(sdataout)+'\n'
      adataout = str(adataout)+'\n'

      for y in range(-int(npx*pxl_den),int(npx*pxl_den)): #y coordenates proof particle  
        y = y/pxl_den

        #-----------------------------------------------------------------
        #simulational map
        #-----------------------------------------------------------------
        if ((typemap == 's') or (typemap == 'b')):
          cenergy = pe.readline() #coulomb potential energy dump file
          
          if (((((x**2)/(a**2))+((y**2)/(b**2))+((z**2)/(c**2)))>1) and (cenergy != 'inf') and (cenergy != '-inf') and (cenergy != '0.0') and (cenergy != '')):
            cenergy = float(cenergy)
            genergy = ((2*G*m)*cenergy)/((q**2)*k) #conversion coulomb energy to gravitational energy
            sdataout = sdataout+str(x)+' '+str(y)+' '+str(genergy)+'\n' #simulational data file position xy and genergy

          else:
            sdataout = sdataout+str(x)+' '+str(y)+' 0.0\n'  #body interior postion xy and 0.0 gravitatinal energy

        #-----------------------------------------------------------------
        #analytical map
        #-----------------------------------------------------------------
        if ((typemap == 'a') or (typemap == 'b')):
          z = 0
          
          Ixy = (1/5)*M*((a**2)+(b**2)) #moment of inertia
          Ixz = (1/5)*M*((a**2)+(c**2)) #moment of inertia
          Iyz = (1/5)*M*((b**2)+(c**2)) #moment of inertia
    
          R = npy.sqrt((x**2)+(y**2)+(z**2))
    
          if ((((x**2)/(a**2))+((y**2)/(b**2))+((z**2)/(c**2)))>1): #body exterior particle data only 
            Ug = -((G*M)/R)-(((G/2)/(R**5))*(((Ixy+Ixz-(2*Iyz))*(x**2))+((Ixy+Iyz-(2*Ixz))*(y**2))+((Ixz+Iyz-(2*Ixy))*(z**2))))

            adataout = str(adataout) + str(x) + ' ' + str(y) + ' ' + str(Ug) + '\n'

          else:  
            adataout = str(adataout) + str(x) + ' ' + str(y) + ' ' + '0.0' + '\n'

        #-----------------------------------------------------------------
        #error
        #-----------------------------------------------------------------
        if (typemap == 'b'):
          
          if ((((x**2)/(a**2))+((y**2)/(b**2))+((z**2)/(c**2)))>1):
            mp = npy.absolute(1 - genergy/Ug)
            edataout = str(edataout) + str(R) + ' ' + str(mp) + '\n'

        #-----------------------------------------------------------------
        #comand verify
        #-----------------------------------------------------------------
        if ((typemap != 's') and (typemap != 'a') and (typemap != 'b')):
          print('Invalid command!')	#test newdata variable text
          selectmap() #restart function
 
  pe.close() #close lammps output file
  return #end function
  
selectmap() #start function
    
#-------------------------------------------------------------
#creating data file and plotting graph  
#------------------------------------------------------------- 
gp.c('set terminal png size 800,800') #set map size
gp.c('set border 0') #remove border
gp.c('unset title') #remove title
gp.c('set tics out') #set tics
gp.c('set tics nomirror') #set tics
gp.c('set key font ",15"') #set label font
gp.c('set tics font ",12"') #set tics font
gp.c('set xlabel font ",15"') #set label font
gp.c('set ylabel font ",15"') #set label font
gp.c('set cblabel font ",15"') #set label font
gp.c('set size ratio 1') #set figure ratio to 1
gp.c('set xlabel "Distancia (m)"')
gp.c('set ylabel "Distancia (m)"')
gp.c('set cblabel "Potencial (Kgm/s^2)"')
gp.c('set rmargin at screen 0.78')
gp.c('set cblabel offset 3.6')
gp.c('set palette defined (0 "#fcffa4", 1 "#f8c932", 3 "#f9950a", 4 "#e35933", 5 "#a83655", 6 "#88226a", 7 "#550f6d", 8 "#1f0c48", 9"#000004")') # set palette colors

gp.c('set xrange ['+str(-int(npx))+':'+str(int(npx))+']') #map x range
gp.c('set yrange ['+str(-int(npx))+':'+str(int(npx))+']') #map y range

if ((typemap == 's') or (typemap == 'b')):   
  with open('/home/bruno/ic/simulations/gmap/py/data/s'+str(a)+'-'+str(b)+'-'+str(c)+'p'+str(cnt)+'mapdata.csv', 'w+') as smapfile: #create data file 
    smapfile.write(sdataout) #write datafile
    
  smapfile.close() #close datafile

  #-------------------------------------------------------------
  #plotting map
  #-------------------------------------------------------------  
  gp.c('set output "/home/bruno/ic/simulations/gmap/py/map/s'+str(a)+'-'+str(b)+'-'+str(c)+'p'+str(cnt)+'map.png"') #set map location
  gp.c('plot "/home/bruno/ic/simulations/gmap/py/data/s'+str(a)+'-'+str(b)+'-'+str(c)+'p'+str(cnt)+'mapdata.csv" using 1:2:3 with image notitle') #plot graph
  
if ((typemap == 'a') or (typemap == 'b')): 
  with open('/home/bruno/ic/simulations/gmap/py/data/a'+str(a)+'-'+str(b)+'-'+str(c)+'p'+str(cnt)+'mapdata.csv', 'w+') as amapfile: #create data file 
    amapfile.write(adataout) #write datafile

  amapfile.close() #close datafile

  #-------------------------------------------------------------
  #plotting map
  #------------------------------------------------------------- 
  gp.c('set output "/home/bruno/ic/simulations/gmap/py/map/a'+str(a)+'-'+str(b)+'-'+str(c)+'p'+str(cnt)+'map.png"') #set map location
  gp.c('plot "/home/bruno/ic/simulations/gmap/py/data/a'+str(a)+'-'+str(b)+'-'+str(c)+'p'+str(cnt)+'mapdata.csv" using 1:2:3 with image notitle') #plot graph

if (typemap == 'b'): 
  with open('/home/bruno/ic/simulations/gmap/py/data/e'+str(a)+'-'+str(b)+'-'+str(c)+'p'+str(cnt)+'mapdata.csv', 'w+') as emapfile: #create data file 
    emapfile.write(edataout) #write datafile

  emapfile.close() #close datafile

  #-------------------------------------------------------------
  #plotting map
  #------------------------------------------------------------- 
  gp.c('set ylabel "Erro %"')
  gp.c('unset xrange')
  gp.c('unset yrange')
  gp.c('set output "/home/bruno/ic/simulations/gmap/py/map/e'+str(a)+'-'+str(b)+'-'+str(c)+'p'+str(cnt)+'map.png"') #set map location
  gp.c('f(x) = a+b*x')
  gp.c('fit f(x) "/home/bruno/ic/simulations/gmap/py/data/e'+str(a)+'-'+str(b)+'-'+str(c)+'p'+str(cnt)+'mapdata.csv" via a,b')
  gp.c('f(x) = a*exp(-b*x)')
  gp.c('fit f(x) "/home/bruno/ic/simulations/gmap/py/data/e'+str(a)+'-'+str(b)+'-'+str(c)+'p'+str(cnt)+'mapdata.csv" via a,b')
  gp.c('plot "/home/bruno/ic/simulations/gmap/py/data/e'+str(a)+'-'+str(b)+'-'+str(c)+'p'+str(cnt)+'mapdata.csv" with points pt 5 lc rgb "#ffe6e6" title "erro", f(x) lw 2 lc rgb "red" title "fit"') #plot graph
