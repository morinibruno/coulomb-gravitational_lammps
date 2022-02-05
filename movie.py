import os
import numpy
import cv2
import glob
import re

frames = glob.glob('frames/simulation/*.png')

def tryint(s):
    try:
        return int(s)
    except:
        return s

def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(frames):
    """ Sort the given list in the way that humans expect.
    """
frames.sort(key=alphanum_key)
    
img_array = []
for filename in frames:
    img = cv2.imread(filename)
    height, width, layers = img.shape
    size = (width,height)
    img_array.append(img)
 
out = cv2.VideoWriter('frames/impact_s.avi',cv2.VideoWriter_fourcc(*'DIVX'), 21, size)
 
for i in range(len(img_array)):
    out.write(img_array[i])
out.release()


