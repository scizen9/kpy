
import time
import shutil
import glob
import sys
import os
import pyfits as pf


def docp(src, dest):

    # Read FITS header
    f = pf.open(src)
    hdr = f[0].header
    f.close()

    # Record copies and standard star observations
    ncp = 0
    nstd = 0

    # Skip test and Focus images
    if 'test' not in hdr['OBJECT'] and 'Focus:' not in hdr['OBJECT']:

        # Copy with preserving metadata (date, etc.)
        shutil.copy2(src, dest)
        print 'copied %s to %s' % (src, dest)
        ncp = 1

        # Check for standard star observations
        if 'STD-' in hdr['OBJECT']:
            nstd = 1

    # Report skipping and type
    else:
        if 'test' in hdr['OBJECT']:
            print 'test file %s not copied' % src
        if 'Focus:' in hdr['OBJECT']:
            print 'Focus file %s not copied' % src

    return (ncp, nstd)


def proc_bias_crrs(reddir,ncp):

    # Default return value
    ret = False

    # Get new listing
    retcode = os.system("~/spy what ifu*.fits > what.list")
    if retcode == 0:

        # Generate new Makefile
        retcode2 = os.system("~/spy plan ifu*.fits")
        if retcode2 == 0:

            # Make bias + bias subtraction
            if ncp < 4:
                retcode3 = os.system("make bias")
            else:
                cmd = "make -j%d bias" % min([ncp,16])
                retcode3 = os.system(cmd)
            if retcode3 == 0:

                # Make CR rejection
                if ncp < 2:
                    retcode4 = os.system("make crrs")
                else:
                    cmd = "make -j%d crrs" % min([ncp,8])
                    retcode4 = os.system(cmd)

                # Success on all fronts!
                if retcode4 == 0:
                    print "bias, crrs processed for %d new images" % ncp
                    ret = True
                # Report failures
                else:
                    print "could not make crrs"
            else:
                print "could not make bias"
        else:
            print "could not make plan"
    else:
        print "could not make what.list"

    return ret

def proc_stds(reddir,ncp):

    # Default return value
    ret = False

    # Make new stds
    startTime = time.time()
    retcode = os.system("make newstds")
    procTime = int(time.time - startTime)

    if retcode == 0:
        print "%d new standard star observations processed in %d s" % (ncp, procTime)
        ret = True

    return ret

def cpnew(srcdir, destdir='./'):

    # Get most recent local ifu image
    lf = sorted(glob.glob(destdir+'/ifu*.fits'))[-1]

    # Get list of source files
    srcfiles = glob.glob(srcdir+'/ifu*.fits')

    # Record copies and standard star observations
    ncp = 0
    nstd = 0

    # Loop over source files
    for f in srcfiles:

        # Do we have a newer file in the source dir?
        if os.stat(f).st_mtime > os.stat(lf).st_mtime:

            # Get ifu image name
            fn = f.split('/')[-1]
            # Call copy
            nc,ns = docp(f,destdir+'/'+fn)
            # Record copies, stds
            ncp += nc
            nstd += ns

    # If we copied any files
    if ncp > 0:
        if not proc_bias_crrs(destdir,ncp):
            print "Error processing bias/crrs"
        elif nstd > 0:
            if not proc_stds(destdir,nstd):
                print "Error processing standard stars"
    else:
        print "no files copied"

    return ncp

def find_recent(destdir,fname):

    # Default return value
    ret = False

    # Get all but the most recent reduced data directories
    redlist = sorted([d for d in glob.glob('/scr2/sedm/redux/20??????') if os.path.isdir(d)])[0:-1]

    # Go back until we find our file
    for d in reversed(redlist):
        src = glob.glob(d+'/'+fname)
        if len(src) == 1:
            ncp = docp(src,distdir+'/'+fname)
            if ncp == 1:
                ret = True
                break

    return ret

def cpcal(dirlist, destdir='./'):

    # Have we got a calibration already?
    if os.path.isfile(destdir+'/cube.npy'):
        return (True, 0)

    # Default return value
    ret = False

    # Get current and previous dates
    sdate = dirlist[-1].split('/')[-1]
    pdate = dirlist[-2].split('/')[-1]

    # Most recent source directory
    srcdir = dirlist[-1]

    # Get list of current raw calibration files
    # (within 10 hours of day changeover)
    fspec = srcdir+"/ifu%s_0*.fits" % sdate
    flist = glob.glob(fspec)

    # Count calibration types
    bias = 0
    bias2 = 0
    dome = 0
    Xe = 0
    Hg = 0
    Cd = 0

    # Loop over file list
    for src in flist:

        # Read FITS header
        f = pf.open(src)
        hdr = f[0].header
        f.close()

        # Get OBJECT and ADCSPEED keywords
        obj = hdr['OBJECT']
        speed = hdr['ADCSPEED']

        # Filter Calibs and avoid test images
        if 'Calib' in obj and not 'test' in obj:
            if 'bias' in obj:
                if speed == 2.:
                    bias2 += 1
                elif speed == 0.1:
                    bias += 1
            if 'dome' in obj: dome += 1
            if 'Xe' in obj:     Xe += 1
            if 'Hg' in obj:     Hg += 1
            if 'Cd' in obj:     Cd += 1

            # Copy cal images
            imf = src.split('/')[-1]
            docp(src,destdir+'/'+imf)

    # Are we missing anything?
    if bias < 5 or bias2 < 5 or dome < 5 or Xe < 5 or Hg < 3 or Cd < 3:

        # If there is a previous night, get those files
        if (int(sdate)-int(pdate)) == 1:

            # Set the previous night as the source directory
            srcdir = dirlist[-2]

            # Get list of previous night's raw calibration files
            # (within four hours of day changeover)
            fspec = srcdir+"/ifu%s_2*.fits" % pdate
            flist = glob.glob(fspec)

            # Loop over file list
            for src in flist:

                # Read FITS header
                f = pf.open(src)
                hdr = f[0].header
                f.close()

                # Get OBJECT and ADCSPEED keywords
                obj = hdr['OBJECT']
                speed = hdr['ADCSPEED']

                # Filter Calibs and avoid test images
                if 'Calib' in obj and not 'test' in obj:
                    if 'bias' in obj:
                        if speed == 2.:
                            bias2 += 1
                        elif speed == 0.1:
                            bias += 1
                    if 'dome' in obj: dome += 1
                    if 'Xe' in obj:     Xe += 1
                    if 'Hg' in obj:     Hg += 1
                    if 'Cd' in obj:     Cd += 1

                    # Copy cal images
                    imf = src.split('/')[-1]
                    docp(src,destdir+'/'+imf)

        # Are we missing anything?
        if bias < 5:
            if find_recent(destdir,'bias0.1.fits'): bias = 10
        if bias2 < 5:
            if find_recent(destdir,'bias2.0.fits'): bias2 = 10
        if dome < 5:
            if find_recent(destdir,'dome.fits'): dome = 5
        if Xe < 5:
            if find_recent(destdir,'Xe.fits'): Xe = 5
        if Hg < 3:
            if find_recent(destdir,'Hg.fits'): Hg = 3
        if Cd < 3:
            if find_recent(destdir,'Xe.fits'): Cd = 3

        # We got them all!
        if bias >= 5 and bias2 >= 5 and dome >= 5 and Xe >= 5 and Hg >= 3 and Cd >= 3:
            ret = True

    # We got them all!
    else:
        ret = True

    return (ret, sum([bias,bias2,dome,Xe,Hg,Cd]))


def go():

    # Get all raw directories
    rawlist = sorted([d for d in glob.glob('/scr2/sedm/raw/20??????') if os.path.isdir(d)])

    # Source directory is most recent raw dir
    srcdir = rawlist[-1]

    # Outpur directory is based on source dir
    outdir = '/scr2/sedm/redux/' + srcdir.split('/')[-1]

    # Make sure reduced directory exists
    if os.path.exists(outdir) == False:
        os.mkdir(outdir)

    # Go there
    os.chdir(outdir)

    # Copy calibration files
    stat, ncp = cpcal(rawlist,outdir)
    if not stat:
        sys.exit("Could not copy calibrations, stopping")

    # Process calibrations
    if ncp > 0:
        startTime = time.time()
        if not proc_bias_crrs(outdir,20):
            sys.exit("Could not do bias,crrs processing, stopping")

        procbTime = int(time.time() - startTime)

        # Process cube
        startTime = time.time()
        retcode = os.system("make cube.npy")
        if retcode != 0:
            sys.exit("Could not generate cube.npy, stopping")

        proccTime = int(time.time() - startTime)

        # Report times
        print "Calibration processing took %d s (bias,crrs) and %d s (cube)" % (procbTime, proccTime)

    # loop and copy new files
    doit = True
    try:
        while doit:
            print "waiting...",
            sys.stdout.flush()
            time.sleep(60)
            print "checking for new ifu images..."
            sys.stdout.flush()
            startTime = time.time()
            ncp = cpnew(srcdir,outdir)
            if ncp > 0:
                procTime = int(time.time() - startTime)
                print "%d new ifu images process in %d s" % (ncp,procTime)
                sys.stdout.flush()

    except KeyboardInterrupt:
        sys.exit("Exiting")


if __name__ == '__main__':
    go()
