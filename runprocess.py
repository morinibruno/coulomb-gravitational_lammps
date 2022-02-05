import os 

cmd = 'nohup python3  < /home/bruno/ic/simulations/gmap/impact/0.2/out_map/py/gmap_run0.py > gmap_run.out & disown'
os.system(cmd) 
cmd = 'nohup python3  < /home/bruno/ic/simulations/gmap/impact/0.2/out_map/py/gmap_run1.py > gmap_run.out & disown'
os.system(cmd) 
