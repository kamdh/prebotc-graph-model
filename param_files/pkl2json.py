#!/usr/bin/env python
# Kameron Decker Harris

import re
import pickle
import json
import glob

for f in glob.glob('./*.pkl'):
    fh = open(f, 'r')
    params = pickle.load(fh)
    fh.close()
    fnew = re.sub(r'(pkl$)', 'json', f)
    fh = open(fnew, 'w')
    json.dump(params, fh)
    fh.close()
