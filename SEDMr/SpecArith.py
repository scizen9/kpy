
import argparse
import numpy as np
import os
from scipy.interpolate import interp1d
import sys


def specDiv(A, B, out):
    """ Divide spectra in A, B and store to out """

    s_A = np.load(A)[0]
    s_B = np.load(B)[0]

    t_A, t_B = s_A['exptime'], s_B['exptime']
    if t_A != t_B:
        print("Exposure times do not match (%s v %s). This may be a problem." %
                (t_A, t_B))

    # Combine keys in rational way
    result = {}
    for k,v in s_A.items():
        # Skip these
        if k in ['meta', 'Extinction Correction', 'doc']: continue
        # Are there discrepant dictionaries?
        if k not in s_B:
            raise Exception("A has key %s but B does not. Quitting" % k)

        # Take average of these
        if "dlam" in k or "extinction_corr" in k:
            result[k] = (s_A[k] + s_B[k])/2.0
        # Combine these lists
        elif "spectra" in k or "object_spaxel_ids" in k:
            result[k] = np.concatenate((s_A[k], s_B[k]))
        # Otherwise just add them together
        else:
            result[k] = s_A[k] + s_B[k]

    # Get interpolated flux
    f2 = interp1d(s_B['nm'], s_B['ph_10m_nm'], bounds_error=0,
                  fill_value=np.nan)

    # Use first for reference wavelengths
    result['nm'] = s_A['nm']
    # Calculate division
    result['ph_10m_nm'] = (s_A['ph_10m_nm'] / f2(s_A['nm']))

    # Average sky
    if s_A.has_key('skyph') and s_B.has_key('skyph'):
        s2 = interp1d(s_B['nm'], s_B['skyph'], bounds_error=0,
                      fill_value=np.nan)
        result['skyph'] = (s_A['skyph'] + s2(s_A['nm']))/2.0

    # Average variance
    if s_A.has_key('var') and s_B.has_key('var'):
        v2 = interp1d(s_B['nm'], s_B['var'], bounds_error=0, fill_value=np.nan)
        result['var'] = (s_A['var'] + v2(s_A['nm']))/2.0

    # Use one of the docs
    result['doc'] = s_A['doc']
    result['operation'] = "specDiv"

    # Record individual meta data
    result['meta_1'] = s_A['meta']
    result['meta_2'] = s_B['meta']
    # Need this for compatibility with other programs
    result['meta'] = s_A['meta']

    np.save(out, [result])


def specAvg(A, B, out):
    """ Add spectra in A, B and store to out """

    s_A = np.load(A)[0]
    s_B = np.load(B)[0]

    t_A, t_B = s_A['exptime'], s_B['exptime']
    if t_A != t_B:
        print("Exposure times do not match (%s v %s). This may be a problem." %
                (t_A, t_B))

    # Combine keys in a rational way
    result = {}
    for k,v in s_A.items():
        # Skip these
        if k in ['operation', 'meta_1', 'meta_2', 'meta', 
                 'Extinction Correction', 'doc']:
            continue
        # Are there discrepant dictionaries?
        if k not in s_B:
            raise Exception("A has key %s but B does not. Quitting" % k)

        # Take average of these
        if "dlam" in k or "extinction_corr" in k:
            result[k] = (s_A[k] + s_B[k])/2.0
        # Combine these lists
        elif "spectra" in k or "object_spaxel_ids" in k:
            result[k] = np.concatenate((s_A[k], s_B[k]))
        # Otherwise just add them together
        elif "quality" in k:
            result[k] = np.max((s_A[k], s_B[k]))
        else:
            result[k] = s_A[k] + s_B[k]

    # Get interpolated flux and sky
    f2 = interp1d(s_B['nm'], s_B['ph_10m_nm'], bounds_error=0,
                  fill_value=np.nan)
    s2 = interp1d(s_B['nm'], s_B['skyph'], bounds_error=0, fill_value=np.nan)

    # Use first for reference wavelengths
    result['nm'] = s_A['nm']
    # Get average flux
    result['ph_10m_nm'] = (s_A['ph_10m_nm'] + f2(s_A['nm']))/2.0
    # Get average sky
    result['skyph'] = (s_A['skyph'] + s2(s_A['nm']))/2.0

    # Get average variance
    if s_A.has_key('var'):
        v2 = interp1d(s_B['nm'], s_B['var'], bounds_error=0, fill_value=np.nan)
        result['var'] = (s_A['var'] + v2(s_A['nm']))/2.0

    # Use one of the docs
    result['doc'] = s_A['doc']
    result['operation'] = "specAdd"

    # Record individual meta data
    result['meta_1'] = s_A['meta']
    result['meta_2'] = s_B['meta']
    # Need this for compatibility with other programs
    result['meta'] = s_A['meta']

    np.save(out, [result])


def specAdd(A, B, out):
    """ Add spectra in A, B and store to out """

    s_A = np.load(A)[0]
    s_B = np.load(B)[0]

    t_A, t_B = s_A['exptime'], s_B['exptime']
    if t_A != t_B:
        print("Exposure times do not match (%s v %s). This may be a problem." %
                (t_A, t_B))

    # Combine keys in a rational way
    result = {}
    for k,v in s_A.items():
        # Skip these
        if k in ['operation', 'meta_1', 'meta_2', 'meta',
                 'Extinction Correction', 'doc']:
            continue
        # Are there discrepant dictionaries?
        if k not in s_B:
            raise Exception("A has key %s but B does not. Quitting" % k)

        # Take average of these
        if "dlam" in k or "extinction_corr" in k:
            result[k] = (s_A[k] + s_B[k])/2.0
        # Combine these lists
        elif "spectra" in k or "object_spaxel_ids" in k:
            result[k] = np.concatenate((s_A[k], s_B[k]))
        # Otherwise just add them together
        elif "quality" in k:
            result[k] = np.max((s_A[k], s_B[k]))
        else:
            result[k] = s_A[k] + s_B[k]

    # Get interpolated flux and sky
    f2 = interp1d(s_B['nm'], s_B['ph_10m_nm'], bounds_error=0,
                  fill_value=np.nan)
    s2 = interp1d(s_B['nm'], s_B['skyph'], bounds_error=0, fill_value=np.nan)

    # Use first for reference wavelengths
    result['nm'] = s_A['nm']
    # Get flux sum
    result['ph_10m_nm'] = s_A['ph_10m_nm'] + f2(s_A['nm'])
    # Get sky sum
    result['skyph'] = s_A['skyph'] + s2(s_A['nm'])

    # Get variance sum
    if s_A.has_key('var'):
        v2 = interp1d(s_B['nm'], s_B['var'], bounds_error=0, fill_value=np.nan)
        result['var'] = s_A['var'] + v2(s_A['nm'])

    # Use A docs
    result['doc'] = s_A['doc']
    result['operation'] = "specAdd"

    # Record individual meta data
    result['meta_1'] = s_A['meta']
    result['meta_2'] = s_B['meta']
    # Need this for compatibility with other programs
    result['meta'] = s_A['meta']

    np.save(out, [result])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        """Perform arithmetic operations on extracted spectra sp_*.npy.
        Operations are:
        'a': average (default)
        '+': add
        '/': divide
        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--operation', type=str, help='Operation to perform',
                        default='a')
    parser.add_argument('A', type=str, help='First term', default="must set A")
    parser.add_argument('B', type=str, help='Second term', default="must set B")
    parser.add_argument('outname', type=str, help='Output name',
                        default="Must set outname")

    args = parser.parse_args()
    err = ""
    if not os.path.isfile(args.A):
        err += "File '%s' does not exist.\n" % args.A
    if not os.path.isfile(args.B):
        err += "File '%s' does not exist.\n" % args.B
    if os.path.isfile(args.outname):
        err += "File '%s' exists. Use a different file name.\n" % args.outname

    if err != "":
        print(err)
        sys.exit(1)

    print "%s %s %s > %s" % (args.A, args.operation, args.B, args.outname)
    if args.operation == '+':
        specAdd(args.A, args.B, args.outname)
    elif args.operation == 'a':
        specAvg(args.A, args.B, args.outname)
    elif args.operation == '/':
        specDiv(args.A, args.B, args.outname)
    else:
        print "%s not recognized as an op." % args.operation
