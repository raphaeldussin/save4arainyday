import PyRaf
import sys

myfile=str(sys.argv[1])
myvar=str(sys.argv[2])

var = PyRaf.readfull(myfile,myvar)

var = PyRaf.mask_value(var,-9999.)

#if (myvar.find('precip')) or (myvar.find('snow')):
if (myvar.find('precip')) :
	print var.min() * 86400 , var.max() * 86400
else:
	print var.min() , var.max()
