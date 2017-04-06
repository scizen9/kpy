import os
import numpy as np
import math
import argparse
import pylab as pl
import SEDMr.Version as Version

drp_ver = Version.ifu_drp_version()


def plot_drp_ver():
    ax = pl.gca()
    ax.annotate('DRP: ' + drp_ver, xy=(0.0, 0.01), xytext=(0, 0),
                xycoords=('axes fraction', 'figure fraction'),
                textcoords='offset points', size=6,
                ha='center', va='bottom')


def sig2noise(norm_exptime=3600., mag_file='mags.dat', ref_mag='R',
              only_exptime=False, no_stamp=False):
    """Calculate a normalized Signal-to-Noise for a set of sp_*.npy files
    Args:
        norm_exptime (float): exposure time to normalize magnitudes to (s)
        mag_file (str): filename for table of magnitudes (see Note below)
        ref_mag (str): reference filter for magnitudes (defaults to 'R')
        only_exptime (bool): set to use norm_exptime as a filter for exptime
        no_stamp (bool): set to prevent DRP stamp from printing
    Notes:
        mag_file should contain one line per object file with the following
        whitespace delimited columns:
            1. Object name (e.g. PTF16evs, or NGC_869-1_obs1 [NB: no spaces!])
            2. B mag
            3. V mag
            4. R mag
            5. I mag
        This procedure will use the mag_file to find the sp_[Object name].npy
        files in the current directory.  It will then normalize the magnitude
        based on the exposure time and the norm_exptime value.  It then plots
        the effective signal-to-noise ratio at that exposure time for bins
        in wavelength of 100 nm starting at 400 nm and ending at 1000 nm.
        """
    # Magnitudes list
    magnames = ['B', 'V', 'R', 'I']
    # Is our ref_mag good?
    if ref_mag in magnames:
        # Get ref_mag index
        mi = magnames.index(ref_mag) + 1
        # Open and read in mags table
        if os.path.exists(mag_file):
            data = []
            with open(mag_file, 'rb') as source:
                for line in source:
                    fields = line.split()
                    data.append(fields)
            # Now loop over data entries
            s2ns = []   # 400-900 nm
            s2n1 = []   # 400-500 nm
            s2n2 = []   # 500-600 nm
            s2n3 = []   # 600-700 nm
            s2n4 = []   # 700-800 nm
            s2n5 = []   # 800-900 nm
            mgs = []
            for line in data:
                if '#' in line[0]:  # Skip comment lines
                    continue
                if float(line[mi]) < 0:  # Skip missing mags
                    continue
                # Get relevant sp_*.npy file
                sp_file = 'sp_' + line[0] + '.npy'

                if os.path.exists(sp_file):
                    # print("Loading %s" % sp_file)
                    sp = np.load(sp_file)[0]
                    # Get normalized magnitude
                    if 'exptime' in sp:
                        expt = sp['exptime']
                        # A/B exposure needs to be doubled
                        if 'spectraB' in sp:
                            expt *= 2.
                    else:
                        print("No exptime found!")
                        continue
                    if only_exptime:
                        if expt != norm_exptime:
                            continue
                    magoff = 2.5 * math.log10(norm_exptime/expt)
                    mag = float(line[mi]) + magoff
                    mgs.append(mag)
                    # Get average signal-to-noise
                    s2n = sp['ph_10m_nm'] / np.sqrt(sp['var'])
                    ll = sp['nm']
                    roi = (ll > 400) & (ll < 900)
                    s2ns.append(np.nanmean(s2n[roi]))
                    roi = (ll > 400) & (ll < 500)
                    s2n1.append(np.nanmean(s2n[roi]))
                    roi = (ll > 500) & (ll < 800)
                    # roi = (ll > 500) & (ll < 600)
                    s2n2.append(np.nanmean(s2n[roi]))
                    roi = (ll > 800) & (ll < 900)
                    # roi = (ll > 600) & (ll < 700)
                    s2n3.append(np.nanmean(s2n[roi]))
                    roi = (ll > 700) & (ll < 800)
                    s2n4.append(np.nanmean(s2n[roi]))
                    roi = (ll > 800) & (ll < 900)
                    s2n5.append(np.nanmean(s2n[roi]))
                    print("%s: expt = %.1f, nmg = %.2f, s/n = %.2f" %
                          (sp_file, expt, mag, s2ns[-1]))
                else:
                    print("No file for object %s" % line[0])
                # Plot results
            if len(mgs) > 3:
                pl.rcParams.update({'font.size': 14})
                print("Plotting %d observations" % len(mgs))
                # pl.plot(mgs, s2ns, 'k.', label='400-900 nm')
                pl.plot(mgs, s2n1, 'm.', label='400-500 nm')
                pl.plot(mgs, s2n2, 'g.', label='500-800 nm')
                pl.plot(mgs, s2n3, 'r.', label='800-900 nm')
                # pl.plot(mgs, s2n1, 'm.', label='400-500 nm')
                # pl.plot(mgs, s2n2, 'b.', label='500-600 nm')
                # pl.plot(mgs, s2n3, 'g.', label='600-700 nm')
                # pl.plot(mgs, s2n4, 'r.', label='700-800 nm')
                # pl.plot(mgs, s2n5, 'y.', label='800-900 nm')
                pl.grid(True)
                pl.xlabel('Equivalent %s Magnitude at %.1f s' %
                          (ref_mag, norm_exptime))
                pl.ylabel('<SNR>')
                pl.legend()
                if not no_stamp:
                    plot_drp_ver()
                pl.show()
                mgs = np.array(mgs)
                s2n1 = np.array(s2n1)
                s2n2 = np.array(s2n2)
                s2n3 = np.array(s2n3)
                inds = mgs.argsort()
                mgs = mgs[inds]
                s2n1 = s2n1[inds]
                s2n2 = s2n2[inds]
                s2n3 = s2n3[inds]
                roi = (mgs > 16.) & (mgs <= 17.)
                print ("R 16-17: %f, %f, %f" % (np.nanmean(s2n1[roi]),
                                                np.nanmean(s2n2[roi]),
                                                np.nanmean(s2n3[roi])))
                roi = (mgs > 17.) & (mgs <= 18.)
                print ("R 17-18: %f, %f, %f" % (np.nanmean(s2n1[roi]),
                                                np.nanmean(s2n2[roi]),
                                                np.nanmean(s2n3[roi])))
                roi = (mgs > 18.) & (mgs <= 19.)
                print ("R 18-19: %f, %f, %f" % (np.nanmean(s2n1[roi]),
                                                np.nanmean(s2n2[roi]),
                                                np.nanmean(s2n3[roi])))
                roi = (mgs > 19.) & (mgs <= 20.)
                print ("R 19-20: %f, %f, %f" % (np.nanmean(s2n1[roi]),
                                                np.nanmean(s2n2[roi]),
                                                np.nanmean(s2n3[roi])))
            else:
                print("Not enough data points: %d" % len(mgs))
        else:
            print("Cannot find %s" % mag_file)
    else:
        print("Can't use %s mags" % ref_mag)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Analyze objects in mags.dat for signal-to-noise""",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--mag_file', type=str, default='mags.dat',
                        help='Input magnitudes file')
    parser.add_argument('--norm_exptime', type=float, default=3600.,
                        help='Exposure time to normalize magnitudes to (s)')
    parser.add_argument('--ref_mag', type=str, default='V',
                        help='Reference mag (B, V, R, or I)')
    parser.add_argument('--only_exptime', action="store_true", default=False,
                        help='Set to use norm_exptime to filter objects')
    parser.add_argument('--no_stamp', action="store_true", default=False,
                        help='Set to prevent plotting DRP version stamp')

    args = parser.parse_args()

    sig2noise(norm_exptime=args.norm_exptime, mag_file=args.mag_file,
              ref_mag=args.ref_mag, only_exptime=args.only_exptime,
              no_stamp=args.no_stamp)
