import os

path = '/home/bruno/ic/simulations/gmap/impact/0.2'

#cmd = "lmp_gserial < "+path+"/lammps/in.impact"
cmd = "mpirun -np 2 lmp_gmpi < "+path+"/lammps/in.impact"
os.system(cmd)

nframes = 320
nprocess = 10
runlines = ''

for mapfiles in range(0,nprocess):

  with open('gmap.py', 'rt') as rgm:
    with open('out_map/py/gmap_run'+str(mapfiles)+'.py', 'wt') as wgm:

  	  for gmap in rgm:

  	  	line = gmap.replace("line 1", "for frames in range(1,"+str(int(nframes/nprocess)+1)+"):").replace(

  	  	"line 2", "with open('"+path+"/out_map/lammps/in/in.map"+str(mapfiles)+"', 'wt') as w_lammpsdata:").replace(
  	  	
  	  	"line 3", "w_lammpsdata.write(data.replace('read_data          data', 'read_data          "+path+"/lammps/data/'+str(int("+str(((nframes/nprocess)*1000)*mapfiles)+"+frames*1000))+'.data group g_1'))").replace(

        "line 4", "with open('"+path+"/out_map/lammps/in/in.map"+str(mapfiles)+"', 'rt') as r_lammps:").replace(

        "line 5", "with open('"+path+"/out_map/lammps/in/in.map_p"+str(mapfiles)+"', 'wt') as w_lammps:").replace(

        "line 6", "'dump               2', 'dump               2 2 custom 1 "+path+"/out_map/lammps/dump/dump_map"+str(mapfiles)+".out x y c_1[*] v_amx v_amy'))").replace(

        "line 7", "cmd = 'lmp_gserial < "+path+"/out_map/lammps/in/in.map_p"+str(mapfiles)+"'").replace(

        "line 8", "with open('"+path+"/out_map/lammps/dump/dump_map"+str(mapfiles)+".out', 'rt') as r_dump:").replace(

        "line 9", "with open('"+path+"/out_map/map_data/map_data"+str(mapfiles)+".csv', 'wt') as w_dump:").replace(

        "line x1", "gp.c('set output \""+path+"/frames/map/'+str(int("+str(((nframes/nprocess)*1000)*mapfiles)+"+(frames*1000)))+'.png\"')").replace(

        "line x2", "gp.c('plot \""+path+"/out_map/map_data/map_data"+str(mapfiles)+".csv\" with image notitle')")

  	  	wgm.write(line)

    wgm.close()
  rgm.close()  
  
  with open('runprocess.py', 'wt') as run:   	

  	runlines = (runlines + "cmd = 'nohup python3  < "+path+"/out_map/py/gmap_run"+str(mapfiles)+".py > gmap_run.out & disown'\nos.system(cmd) \n")

  run.close()

with open('runprocess.py', 'wt') as run:   

  run.write("import os \n\n")
  run.write(runlines)

run.close()

cmd = 'nohup python3 < runprocess.py > ruprocess.out &  disown'
os.system(cmd)
