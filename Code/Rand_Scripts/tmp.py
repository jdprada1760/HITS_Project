import numpy as np
import ctypes

indata = np.ones((10,3), dtype = np.float32)
outdata = np.zeros((10,3), dtype = np.float32)

print("in: %s" %indata)
print("out: %s" %outdata)
print("--------------------------------------------------------")

lib = ctypes.cdll.LoadLibrary('./C_Libs/inertia.so')

# Calls C function
print(len(indata))
lib.cfun(ctypes.c_void_p(indata.ctypes.data),ctypes.c_int(3*len(indata)), ctypes.c_void_p(outdata.ctypes.data))
print("in: %s" %indata)
print("out: %s" %outdata)
