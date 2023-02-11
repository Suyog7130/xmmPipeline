
"""
For conversion of the RA and DEC.
"""

import argparse
import numpy as np

if __name__=="__main__":
    
    parser = argparse.ArgumentParser(description='Convert HH:MM:SS coordinates to decimal.')
    
    subparsers = parser.add_subparsers()
    method1 = subparsers.add_parser('method')
    method1.add_argument('--obsIDs', nargs='+', action='store')
    #method1.set_defaults(func=hello)
    
    parser.add_argument('--ra', action='store', default='12:48:15', \
                        help='RA in HH:MM:SS\n default is the RA of ASASSN-14li, %(default)s')
    parser.add_argument('--dec', action='store', default='+17:46:26.20', \
                        help='DEC in DD:MM:SS\n default is the DEC of ASASSN-14li, %(default)s')
    args = parser.parse_args()
    
    print(args)
    if 'func' in args:
        args.func(args.obsIDs)
    
    if ':' in list(args.ra):
        h, m, s = np.array(args.ra.split(':')).astype(float)
        ra = (h + m/60 + s/3600 )*360/24
    else:
        ra = float(args.ra)
        
    if ':' in list(args.dec):
        d, m, s = np.array(args.dec.split(':')).astype(float)
        dec = np.abs(d) + m/60 + s/3600
        if d<0:
            dec = float('-'+str(dec))
    else:
        dec = float(args.dec)
    
    print('RA:{}, DEC:{}'.format(ra, dec))
    
    
