import os
import numpy
import PyGnuplot as gp
import cv2
import glob

line 1

  outline=''
  
  with open('/home/bruno/ic/simulations/gmap/impact/0.2/lammps/in.map', 'rt') as r_lammpsdata:
    line 2
    
      for data in r_lammpsdata:
        line 3

    w_lammpsdata.close()
  r_lammpsdata.close()

  movex = 0.00
  movey = 0.00

  for x in range(-20,80):

    for y in range(-50,50):

      line 4 
        line 5

          for line in r_lammps:
            w_lammps.write(line.replace('create_atoms       x', 'create_atoms       ' 
             +'1 single '+str((x)+(frames*movex))+' '+str((y)+(frames*movey))+' '+'0.0').replace(
             line 6

        w_lammps.close()
      r_lammps.close()

      line 7
      #cmd = 'mpirun -np 4 lmp_gravity < in.ellipsoid_map_data_pc'
      os.system(cmd)

      line 8
           
        for i in range(0,10):
          dump = r_dump.readline()

        outline = (str(outline)+str(dump))
         
      r_dump.close()

  line 9
    w_dump.write(outline)
    
  w_dump.close()

  gp.c('set terminal png size 800,800')
  line x1
  gp.c('unset xtics')
  gp.c('unset ytics')
  gp.c('set border 0')
  gp.c('unset title')
  gp.c('set size ratio 1')
  gp.c('set cbrange[0:-2.5e3]')
  gp.c('set palette defined (0 "#fcffa4", 1 "#f8c932", 2 "#f8c932", 3 "#f9950a", 4 "#e35933", 5 "#a83655", 6 "#88226a", 7 "#550f6d", 8 "#1f0c48", 9"#000004")')
  line x2
    


