import numpy as np
from scipy import weave
from scipy.weave import converters
import os

A = np.zeros([10,10], dtype=float);
B = np.zeros([10,10], dtype=float);
teststr = """
int i,j;
for (i = 0; i < NA[0]; i++) {
for (j = 0; j < NA[1]; j++) {
A(i,j) = i*j;
B(i,j) = A(i,j)/2;
}
}
"""
thevars=set(dir()) - set(dir(__builtins__))
print thevars
weave.inline(teststr, ['A', 'B'], 
             type_converters = converters.blitz, 
             compiler='gcc')
print A
print B

x = -.1
teststr = """
return_val = infHN(-0.036, -0.0085, x);
"""
support_code = '''
extern "C"{
double infHN( double A, double B, double V );
double tau( double A, double B, double C, double V );
double phi( double x, double kNa );
}
'''
tmp = weave.inline(teststr, ['x'], support_code=support_code,
                   type_converters = converters.blitz, 
                   compiler='gcc',
                   library_dirs = [os.getcwd()],
                   libraries=['helpers'],
                   headers=['<math.h>'],
                   runtime_library_dirs = [os.getcwd()])
print tmp


