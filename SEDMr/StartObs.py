"""Conduct automatic reduction of SEDM data in sedm@pharos.

Usage:
    ~/spy ~/kpy/SEDMr/StartObs.py
    no arguments are required, this will run the 'go' function.
Functions:
    go      - outer loop waits for new data directory in /scr2/sedm/raw
    ObsLoop - one night observing loop: processes calibrations and science data
    cpcal   - copies calibration images into redux directory
    cpprecal    - copies calibration images from previous day directory
    find_recent - finds the most recent processed calibration file
    cpnew       - copies new images files into redux directory
    proc_stds   - processes standard star observations
    proc_bias_crrs  - processes biases and CR rejection
    docp        - low level copy routine that avoids test and focus images
    cal_proc_ready  - check if all required raw cal images are present
    cal_ready   - check if all required processed cal files are present
    cal_reset   - reset the status of calibration files
Globals:
    nbias,nbias2,nXe,ndome,nHg,nCd  - current count of raw cal files
    CalProcReady    - are all required raw cal images present?
    CalReady        - are all required processed cal files present?
"""
import time
import shutil
import glob
import sys
import os
import pyfits as pf
import argparse

nbias2 = 0
nbias = 0
nXe = 0
ndome = 0
nHg = 0
nCd = 0
CalProcReady = False
CalReady = False
CalPrevious = False
BiasReady = False


def cal_reset():
    """Reset counts of raw calibration files, and calibration status."""

    # Global variables
    global nbias, nbias2, ndome, nXe, nHg, nCd, \
            CalProcReady, CalReady, CalPrevious, BiasReady
    # Reset calibration file counters
    nbias2 = 0
    nbias = 0
    nXe = 0
    ndome = 0
    nHg = 0
    nCd = 0
    CalProcReady = False
    CalReady = False
    BiasReady = False
    CalPrevious = False


def cal_ready(reddir='./'):
    """Check for all required calibration files in input reduced directory."""

    # Global variables
    global CalReady, BiasReady
    # Do we have all the calibration files?
    if (os.path.exists(os.path.join(reddir,'bias0.1.fits')) and
        os.path.exists(os.path.join(reddir,'bias2.0.fits'))):
        BiasReady = True
        if (os.path.exists(os.path.join(reddir,'cube.npy')) and 
            os.path.exists(os.path.join(reddir,'flat-dome-700to900.npy'))):
            CalReady = True
        else:
            CalReady = False
    else:
        CalReady = False
        BiasReady = False


def cal_proc_ready():
    """Check counts for all required raw cal file types in current directory."""

    # Global variables
    global nbias, nbias2, ndome, nXe, nHg, nCd, CalProcReady
    # Do we have all the calibration files?
    if (nbias2 > 5 and nbias > 5 and nXe >= 5 and ndome >= 5 and 
            nHg >= 3 and nCd >= 3):
        CalProcReady = True


def docp(src, dest):
    """Low level copy from raw directory to redux directory

    Call:
        docp(src, dest)
    Inputs:
        src     - source file
        dest    - destination directory
    Procedure:
        Checks for raw cal files, standard star observations and
        increments their respective counters, while avoiding any
        test and focus images. Uses shutil.copy2 to do the copying.
    Returns:
        (ncp, nstd) - number of images copied, number of standard
        star images copied
    """

    # Global variables
    global nbias, nbias2, ndome, nXe, nHg, nCd
    # Read FITS header
    f = pf.open(src)
    hdr = f[0].header
    f.close()
    # Get OBJECT and ADCSPEED keywords
    obj = hdr['OBJECT']
    speed = hdr['ADCSPEED']
    # Record copies and standard star observations
    ncp = 0
    nstd = 0
    # Skip test and Focus images
    if 'test' not in obj and 'Focus:' not in obj:
        # Copy with preserving metadata (date, etc.)
        shutil.copy2(src, dest)
        print 'copied %s to %s' % (src, dest)
        ncp = 1
        # Check for standard star observations
        if 'STD-' in obj:
            nstd = 1
        # Check for calibration files
        elif 'Calib' in obj:
            if 'bias' in obj:
                if speed == 2.0: nbias2 += 1
                if speed == 0.1: nbias += 1
            if 'Xe' in obj: nXe += 1
            if 'dome' in obj: ndome += 1
            if 'Hg' in obj: nHg += 1
            if 'Cd' in obj: nCd += 1
    # Report skipping and type
    else:
        if 'test' in hdr['OBJECT']:
            print 'test file %s not copied' % src
        if 'Focus:' in hdr['OBJECT']:
            print 'Focus file %s not copied' % src

    return (ncp, nstd)


def proc_bias_crrs(reddir='./',ncp=1):
    """Process biases and CR rejection

    Call:
        proc_bias_crrs(reddir, ncp)
    Inputs:
        reddir  - redux directory in which to process
        ncp     - number of files to debias and CR clean
    Returns:
        True if processing was successful, otherwise False
    """

    # Are we using old calib files?
    global CalPrevious
    # Default return value
    ret = False
    # Get new listing
    retcode = os.system("~/spy what ifu*.fits > what.list")
    if retcode == 0:
        # Generate new Makefile
        # Are we using a previous calibration set?
        if CalPrevious:
            retcode2 = os.system("~/spy plan2 ifu*.fits")
        # This calibration set has been generated here
        else:
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
    """Process standard star observations

    Call:
        proc_stds(reddir, ncp)
    Inputs:
        reddir  - redux directory in which to process
        ncp     - number of files to debias and CR clean
    Returns:
        True if processing was successful, otherwise False
    """

    # Default return value
    ret = False
    # Make new stds
    startTime = time.time()
    retcode = os.system("make newstds")
    procTime = int(time.time() - startTime)
    # Did it work?
    if retcode == 0:
        print("%d new standard star observations processed in %d s" % 
                (ncp, procTime))
        ret = True

    return ret


def cpnew(srcdir, destdir='./'):
    """Copies new files from srcdir to destdir.

    Call:
        cpnew(srcdir, destdir) # destdir defaults to current dir ('./')
    Inputs:
        srcdir  - source directory (typically in /scr2/sedm/raw)
        destdir - destination directory (typically in /scr2/sedm/redux)
    Procedure:
        Searches for most recent ifu image in destdir and looks for and
        copies any ifu images in srcdir that are newer and complete.
        Then bias subtracts and CR rejects the images.  If any are standard
        star observations, process them as well.
    Returns:
        ncp     - Number of ifu images actually copied
    """

    # Global variables
    global CalReady, BiasReady
    # Get files in destination directory
    dflist = sorted(glob.glob(os.path.join(destdir,'ifu*.fits')))
    # Are there any files yet?
    if len(dflist) > 0:
        # Get most recent local ifu image
        lf = dflist[-1]
    else:
        lf = None
    # Record copies and standard star observations
    ncp = 0
    nstd = 0
    # Get list of source files
    srcfiles = glob.glob(os.path.join(srcdir,'ifu*.fits'))
    # Loop over source files
    for f in srcfiles:
        # Do we copy?
        doCopy = False
        # Is our source file complete?
        if os.stat(f).st_size >= 8398080:
            if lf is not None:
                # Do we have a newer file in the source dir?
                if os.stat(f).st_mtime > os.stat(lf).st_mtime:
                    doCopy = True
            # No files yet in dest, so all in source are needed
            else:
                doCopy = True
        if doCopy:
            # Get ifu image name
            fn = f.split('/')[-1]
            # Call copy
            nc,ns = docp(f,destdir+'/'+fn)
            # Record copies, stds
            ncp += nc
            nstd += ns
    # We copied files
    print "Copied %d files" % ncp
    # Are we ready to process biases?
    if BiasReady and ncp > 0:
        # Do bias subtraction, CR rejection
        if not proc_bias_crrs(destdir,ncp):
            print "Error processing bias/crrs"
        # Process any standard stars
        elif CalReady and nstd > 0:
            if not proc_stds(destdir,nstd):
                print "Error processing standard stars"
        elif not CalReady and nstd > 0:
            print("Copied %d std star obs, but need cube and flat "
                  "to process further" % nstd)
    # Bias not ready yet
    elif not BiasReady and ncp > 0:
        print "Need biases to process further"

    return ncp


def find_recent(redd,fname,destdir):
    """Find the most recent version of fname and copy it to destdir

    Call:
        find_recent(redd,fname,destdir)
    Inputs:
        redd    - reduced directory (something like /scr2/sedm/redux)
        fname   - what file to look for
        destdir - where the file should go
    Procedure:
        Look through sorted list of redux directories to find most recent
        version of the input file.  Copy it to the destination directory
    Returns:
        True if file found and copied, False otherwise.
    """

    # Default return value
    ret = False
    # Make sure the file doesn't already exist in destdir
    local_file = glob.glob(os.path.join(destdir,fname))
    if len(local_file) == 1:
        print "%s already exists in %s" % (fname, destdir)
        ret = True
    # Search in redd for file
    else:
        # Get all but the most recent reduced data directories
        fspec = os.path.join(redd,'20??????')
        redlist = sorted([d for d in glob.glob(fspec) 
                            if os.path.isdir(d)])[0:-1]
        # Go back in reduced dir list until we find our file
        for d in reversed(redlist):
            src = glob.glob(os.path.join(d,fname))
            if len(src) == 1:
                shutil.copy(src[0], destdir)
                ret = True
                print("Found %s in directory %s, copying to %s" % 
                        (fname, d, destdir))
                break

    return ret


def cpprecal(dirlist, destdir='./'):
    """Copy raw cal files from previous date's directory

    Call:
        cpprecal(dirlist, destdir)  # destdir defaults to current dir ('./')
    Inputs:
        dirlist - a list of raw directories (typically in /scr2/sedm/raw)
        destdir - where to put the files
    Procedure:
        Make sure we only look in previous day directory for files created
        within four hours of the day changeover for possible raw calibration
        files required for the current night's calibration.  Copy any such
        files into the destination directory.
    Returns:
        ncp     - number of images actually copied
    """

    # Reset calibration file counters
    cal_reset()
    # Get current and previous dates
    sdate = dirlist[-1].split('/')[-1]
    pdate = dirlist[-2].split('/')[-1]
    # Record how many images copied
    ncp = 0
    # If there is a previous night, get those files
    if (int(sdate)-int(pdate)) == 1:
        # Set the previous night as the source directory
        srcdir = dirlist[-2]
        # Get list of previous night's raw calibration files
        # (within four hours of day changeover)
        fspec = os.path.join(srcdir, "ifu%s_2*.fits" % pdate)
        flist = glob.glob(fspec)
        # Loop over file list
        for src in flist:
            if os.stat(src).st_size >= 8400000:
                # Read FITS header
                f = pf.open(src)
                hdr = f[0].header
                f.close()
                # Get OBJECT and ADCSPEED keywords
                obj = hdr['OBJECT']
                # Filter Calibs and avoid test images
                if 'Calib' in obj and not 'test' in obj:
                    # Copy cal images
                    imf = src.split('/')[-1]
                    nc, ns = docp(src,os.path.join(destdir,imf))
                    ncp += nc
            else:
                print "Truncated file: %s" % src
    # Check if we got all the calibration files
    cal_proc_ready()

    return ncp


def cpcal(srcdir, destdir='./'):
    """Copy raw cal files from srcdir into destdir.

    Call:
        cpcal(srcdir, destdir)  # destdir defaults to current directory ('./')
    Inputs:
        srcdir  - source for raw cal images, typically in /scr2/sedm/raw
        destdir - place to put the cal images typically in /scr2/sedm/redux
    Procedure:
        Find calibration files taken within 10 hours of the day changeover
        and copy them to the destination directory.
    Returns:
        ncp     - number of images actually copied
    """

    # Get current date
    sdate = srcdir.split('/')[-1]
    # Get list of current raw calibration files
    # (within 10 hours of day changeover)
    fspec = os.path.join(srcdir,"ifu%s_0*.fits" % sdate)
    flist = glob.glob(fspec)
    # Record number copied
    ncp = 0
    # Loop over file list
    for src in flist:
        if os.stat(src).st_size >= 8400000:
            # Read FITS header
            f = pf.open(src)
            hdr = f[0].header
            f.close()
            # Get OBJECT and ADCSPEED keywords
            try:
                obj = hdr['OBJECT']
            except:
                obj = ''
            # Filter Calibs and avoid test images
            if 'Calib' in obj and not 'test' in obj:
                # Copy cal images
                imf = src.split('/')[-1]
                nc, ns = docp(src,os.path.join(destdir,imf))
                ncp += nc
        else:
            print "Truncated file: %s" % src
    # Check if calibrations are read
    cal_proc_ready()

    return ncp


def ObsLoop(rawlist=None, redd=None):
    """One night observing loop: processes calibrations and science data

    Call:
        ObsLoop(rawlist,redd)
    Inputs:
        rawlist - list of raw data directories (usually in /scr2/sedm/raw)
        redd    - reduced directory (something like /scr2/sedm/redux)
    Procdure:
        Copy raw cal files until we are ready to process the night's
        calibrations.  When ready, process them.  If we get to UT = 3:00,
        the sun has set and we then retrieve the processed cal files from
        the most recent night.  Next, enter a loop that waits for new ifu
        images, copies them and performs the basic bias removal and CR
        rejection processing.  If there are standard star observations,
        re-generate the standard star calibration.  If we get no new ifu
        files and we get to UT = 15:00, then the sun is up and we exit
        the loop.
    Returns:
        True if night completed normall, False otherwise
    Note:
        KeyboardInterrupt handler exits gracefully with a ctrl-C.
    """

    # Global variables
    global CalProcReady, CalReady, CalPrevious, BiasReady
    # Default return value
    ret = False
    # Source directory is most recent raw dir
    srcdir = rawlist[-1]
    # Outpur directory is based on source dir
    outdir = os.path.join(redd,srcdir.split('/')[-1])
    # Do we have a new directory?  This tells us we are observing tonight
    if os.path.exists(outdir) == False:
        # Make it
        os.mkdir(outdir)
    # Go there
    os.chdir(outdir)
    # Are there good calibrations there?
    cal_ready(outdir)
    # If not, process them
    if not CalReady:
        # Copy calibration files from previous date directory
        npre = cpprecal(rawlist, outdir)
        if npre > 0:
            print("bias2.0: %d, bias0.1: %d, dome: %d, Xe: %d, Hg: %d, Cd: %d" %
                    (nbias2, nbias, ndome, nXe, nHg, nCd))
        # Now loop until we have calibrations we need
        while not CalProcReady:
            ncp = cpcal(srcdir, outdir)
            print("bias2.0: %d, bias0.1: %d, dome: %d, Xe: %d, Hg: %d, Cd: %d" %
                    (nbias2, nbias, ndome, nXe, nHg, nCd))
            # Not ready yet
            if not CalProcReady:
                # Wait a minute
                print "checking %s for new calibration files..." % srcdir,
                sys.stdout.flush()
                time.sleep(60)
                # Check to see if we are definitely after sunset
                gm = time.gmtime()
                if gm.tm_hour >= 3:
                    # Get earlier calibrations so we can proceed
                    print("It's getting late! UT = %02d:%02d >= 03:00" %
                            (gm.tm_hour, gm.tm_min))
                    print "Let's get our calibrations from a previous night"
                    ncf = find_recent(redd,'fine.npy',outdir)
                    ncc = find_recent(redd,'cube.npy',outdir)
                    ncd = find_recent(redd,'flat-dome-700to900.npy',outdir)
                    ncb = find_recent(redd,'bias0.1.fits',outdir)
                    nc2 = find_recent(redd,'bias2.0.fits',outdir)
                    # Check for failure
                    if not ncf or not ncc or not ncd or not ncb or not nc2:
                        msg = "Calibration stage failed: fine = %s, " \
                                "cube = %s, flat = %s, bias0.1 = %s, " \
                                "bias2.0 = %s, " \
                                "stopping" % (ncf, ncc, ncd, ncb, nc2)
                        sys.exit(msg)
                    # If we get here, we are done
                    CalReady = True
                    BiasReady = True
                    CalPrevious = True
                    break
                else:
                    print("UT = %02d:%02d, still less than 03:00, "
                            "so keep waiting" % (gm.tm_hour, gm.tm_min))
        # Process calibrations if we are using them
        if CalProcReady and not CalReady:
            # bias subtract and CR reject
            startTime = time.time()
            if not proc_bias_crrs(outdir,20):
                sys.exit("Could not do bias,crrs processing, stopping")
            procbTime = int(time.time() - startTime)
            # check status
            cal_ready(outdir)
            # Process cube
            startTime = time.time()
            retcode = os.system("make cube.npy")
            proccTime = int(time.time() - startTime)
            if (os.path.exists(os.path.join(outdir,'cube.npy'))):
                # Process flat
                startTime = time.time()
                retcode = os.system("make flat-dome-700to900.npy")
                if not (os.path.exists(
                        os.path.join(outdir,'flat-dome-700to900.npy'))):
                    print "Making of flat-dome-700to900.npy failed!"
            else:
                print "Making of fine.npy and cube.npy failed!"
            procfTime = int(time.time() - startTime)
            # check status
            cal_ready(outdir)
            if not CalReady:
                print "These calibrations failed!"
                print "Let's get our calibrations from a previous night"
                ncf = find_recent(redd,'fine.npy',outdir)
                ncc = find_recent(redd,'cube.npy',outdir)
                ncd = find_recent(redd,'flat-dome-700to900.npy',outdir)
                ncb = find_recent(redd,'bias0.1.fits',outdir)
                nc2 = find_recent(redd,'bias2.0.fits',outdir)
                # Check for failure
                if not ncf or not ncc or not ncd or not ncb or not nc2:
                    msg = "Calibration stage failed: fine = %s, cube = %s, " \
                            "flat = %s, bias0.1 = %s, bias2.0 = %s, " \
                            "stopping" % (ncf, ncc, ncd, ncb, nc2)
                    sys.exit(msg)
                # If we get here, we are done
                CalReady = True
                BiasReady = True
                CalPrevious = True
            else:
                # Report times
                print("Calibration processing took "
                    "%d s (bias,crrs), %d s (cube), and %d s (flat)" % 
                    (procbTime, proccTime, procfTime))
        else:
            print("Using previous calibration files bias0.1.fits, bias2.0.fits,"
                    " cube.npy, flat-dome-700to900.npy")
    else:
        print "Calibrations already present in %s" % outdir
    # Keep track of no copy
    nnc = 0
    # loop and copy new files
    doit = True
    try:
        while doit:
            # Wait a minute
            print "waiting for new ifu images...",
            sys.stdout.flush()
            time.sleep(60)
            # Check for new ifu images
            print "checking %s for new ifu images..." % srcdir
            sys.stdout.flush()
            # Record starting time for new file processing
            startTime = time.time()
            ncp = cpnew(srcdir,outdir)
            # We copied some new ones so report processing time
            if ncp > 0:
                procTime = int(time.time() - startTime)
                print "%d new ifu images processed in %d s" % (ncp,procTime)
                sys.stdout.flush()
                nnc = 0
            else:
                nnc += 1
            # Have we been waiting for a while?
            if nnc > 3:
                # Check time
                gm = time.gmtime()
                if gm.tm_hour > 15:
                    # No new observations and sun is probably up!
                    print("No new images for %d minutes and UT = %02d:%02d > "
                            "15:59 so sun is probably up!" % 
                            (nnc, gm.tm_hour, gm.tm_min))
                    print "Time to wait until we have a new raw directory"
                    doit = False
                    # Normal termination
                    ret = True
                else:
                    print("No new image for %d minutes but UT = %02d:%02d < "
                            "16:00, so keep waiting" % 
                            (nnc, gm.tm_hour, gm.tm_min))
    # Handle a ctrl-C
    except KeyboardInterrupt:
       sys.exit("Exiting")

    return ret


def go(rawd='/scr2/sedm/raw', redd='/scr2/sedm/redux'):
    """Outermost infinite loop that watches for a new raw directory.

    Call:
        go(rawd='/scr2/sedm/raw', redd='/scr2/sedm/redux')
    Inputs:
        rawd    - raw directory (string), should be /scr2/sedm/raw
        redd    - reduced directory (string), should be like /scr2/sedm/redux
    Procedure:
        Keep a list of raw directories in redd and fire off
        the ObsLoop procedure when a new directory appears.  Check for
        a new raw directory every 10 minutes.
    Returns:
        None
    Note:
        KeyboardInterrupt handler exits gracefully with a ctrl-C.
    """

    # Infinite loop
    dobs = True
    # Keep track of iterations
    its = 0
    # Get all raw directories
    fspec = os.path.join(rawd,'20??????')
    rawlist = sorted([d for d in glob.glob(fspec) if os.path.isdir(d)])
    nraw = len(rawlist)
    print("Found %d raw directories in %s: putting reduced data in %s" %
            (nraw, rawd, redd))
    try:
        while dobs:
            stat = ObsLoop(rawlist, redd)
            if stat:
                its += 1
                print("Finished SEDM observing iteration %d in raw dir %s" %
                        (its, rawlist[-1]))
                print "Now we wait until we get a new raw directory"
                waiting = True
                while waiting:
                    print "waiting for new raw directory..."
                    sys.stdout.flush()
                    time.sleep(600)
                    # Get all raw directories
                    new_rawlist = sorted([d for d in glob.glob(fspec)
                                            if os.path.isdir(d)])
                    new_nraw = len(new_rawlist)
                    if new_nraw > nraw:
                        waiting = False
                        sys.stdout.flush()
                        rawlist = new_rawlist
                        nraw = new_nraw
                        print("Starting next SEDM observing iteration with "
                                "raw dir %s" % rawlist[-1])
                    else:
                        gm = time.gmtime()
                        print("UT = %02d:%02d No new directories yet, "
                                "so keep waiting" % (gm.tm_hour, gm.tm_min))
                        sys.stdout.flush()
    # Handle a ctrl-C
    except KeyboardInterrupt:
       sys.exit("Exiting")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
            """StartObs.py

            """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--rawdir', type=str, default='/scr2/sedm/raw',
            help='Input raw directory (/scr2/sedm/raw)')
    parser.add_argument('--reduxdir', type=str, default='/scr2/sedm/redux',
            help='Output reduced directory (/scr2/sedm/redux)')

    args = parser.parse_args()

    go(rawd=args.rawdir, redd=args.reduxdir)

