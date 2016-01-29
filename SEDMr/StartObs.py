
import time
import shutil
import glob
import sys
import os
import pyfits as pf

nbias2 = 0
nbias = 0
nXe = 0
ndome = 0
nHg = 0
nCd = 0
CalProcReady = False
CalReady = False


def cal_reset():

    # Global variables
    global nbias, nbias2, ndome, nXe, nHg, nCd, CalProcReady, CalReady
    # Reset calibration file counters
    nbias2 = 0
    nbias = 0
    nXe = 0
    ndome = 0
    nHg = 0
    nCd = 0
    CalProcReady = False
    CalReady = False


def cal_ready(reddir='./'):

    # Global variables
    global CalReady
    # Do we have all the calibration files?
    if (os.path.exists(os.path.join(reddir,'cube.npy')) and 
            os.path.exists(os.path.join(reddir,'flat-dome-700to900.npy'))):
        CalReady = True
    else:
        CalReady = False


def cal_proc_ready():

    # Global variables
    global nbias, nbias2, ndome, nXe, nHg, nCd, CalProcReady
    # Do we have all the calibration files?
    if (nbias2 > 5 and nbias > 5 and nXe >= 5 and ndome >= 5 and 
            nHg >= 3 and nCd >= 3):
        CalProcReady = True


def docp(src, dest):

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
    procTime = int(time.time() - startTime)
    # Did it work?
    if retcode == 0:
        print("%d new standard star observations processed in %d s" % 
                (ncp, procTime))
        ret = True

    return ret


def cpnew(srcdir, destdir='./'):

    # Global variables
    global nbias, nbias2, ndome, nXe, nHg, nCd, CalReady
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
    if CalReady and ncp > 0:
        # Do bias subtraction, CR rejection
        if not proc_bias_crrs(destdir,ncp):
            print "Error processing bias/crrs"
        # Process any standard stars
        elif nstd > 0:
            if not proc_stds(destdir,nstd):
                print "Error processing standard stars"
    # We didn't copy files
    else:
        print "no files copied"

    return ncp


def find_recent(destdir,fname):

    # Default return value
    ret = False
    # Get all but the most recent reduced data directories
    redlist = sorted([d for d in glob.glob('/scr2/sedm/redux/20??????') 
                        if os.path.isdir(d)])[0:-1]
    # Go back in reduced dir list until we find our file
    for d in reversed(redlist):
        src = glob.glob(os.path.join(d,fname))
        if len(src) == 1:
            shutil.copy2(src[0], destdir)
            ret = True
            print "Found %s in directory %s" % (fname, d)
            break

    return ret


def cpprecal(dirlist, destdir='./'):

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
    # Check if we got all the calibration files
    cal_proc_ready()

    return ncp


def cpcal(srcdir, destdir='./'):

    # Get current date
    sdate = srcdir.split('/')[-1]
    # Get list of current raw calibration files
    # (within 10 hours of day changeover)
    fspec = os.path.join(srcdir,"ifu%s_0*.fits" % sdate)
    flist = glob.glob(fspec)
    print "n files: %d, filespec: %s" % (len(flist), fspec)
    # Record number copied
    ncp = 0
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
            # Copy cal images
            imf = src.split('/')[-1]
            nc, ns = docp(src,os.path.join(destdir,imf))
            ncp += nc
    # Check if calibrations are read
    cal_proc_ready()

    return ncp


def ObsLoop(rawlist=None):

    # Global variables
    global CalProcReady, CalReady
    # Default return value
    ret = False
    # Source directory is most recent raw dir
    srcdir = rawlist[-1]
    # Outpur directory is based on source dir
    outdir = os.path.join('/scr2/sedm/redux/',srcdir.split('/')[-1])
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
            print("bias2.0:%d, bias0.1:%d, dome:%d, Xe:%d, Hg:%d, Cd:%d" %
                    (nbias2, nbias, ndome, nXe, nHg, nCd))
        # Now loop until we have calibrations we need
        while not CalProcReady:
            ncp = cpcal(srcdir, outdir)
            print("bias2.0:%d, bias0.1:%d, dome:%d, Xe:%d, Hg:%d, Cd:%d" %
                    (nbias2, nbias, ndome, nXe, nHg, nCd))
            # Not ready yet
            if not CalProcReady:
                # Wait a minute
                print "waiting for new calibration files...",
                sys.stdout.flush()
                time.sleep(60)
                # Check to see if we are definitely after sunset
                gm = time.gmtime()
                if gm.tm_hour >= 3:
                    # Get earlier calibrations so we can proceed
                    print("It's getting late! UT = %02d:%02d >= 03:00" %
                            (gm.tm_hour, gm.tm_min))
                    print "Let's get our calibrations from a previous night"
                    ncc = find_recent(outdir,'cube.npy')
                    ncf = find_recent(outdir,'flat-dome-700to900.npy')
                    # Check for failure
                    if not ncc or not ncf:
                        msg = "Calibration stage failed: cube = %s, " \
                                "flat = %s, stopping" % (ncc, ncf)
                        sys.exit(msg)
                    # If we get here, we are done
                    CalReady = True
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
            # Process cube
            startTime = time.time()
            retcode = os.system("make cube.npy")
            if retcode != 0:
                sys.exit("Could not generate cube.npy, stopping")
            proccTime = int(time.time() - startTime)
            # Process flat
            startTime = time.time()
            retcode = os.system("make flat-dome-700to900.npy")
            if retcode != 0:
                sys.exit("Could not generate flat-dome-700to900.npy, stopping")
            procfTime = int(time.time() - startTime)
            # We are done!
            CalReady = True
            # Report times
            print("Calibration processing took "
                  "%d s (bias,crrs), %d s (cube), and %d s (flat)" % 
                  (procbTime, proccTime, procfTime))
        else:
            print("Using previous calibration files "
                  "cube.npy, flat-dome-700to900.npy")
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
            print "checking for new ifu images..."
            sys.stdout.flush()
            # Record starting time for new file processing
            startTime = time.time()
            ncp = cpnew(srcdir,outdir)
            # We copied some new ones so report processing time
            if ncp > 0:
                procTime = int(time.time() - startTime)
                print "%d new ifu images process in %d s" % (ncp,procTime)
                sys.stdout.flush()
                nnc = 0
            else:
                nnc += 1
            # Have we been waiting for a while?
            if nnc > 3:
                # Check number of no copies and time
                gm = time.gmtime()
                if gm.tm_hour > 15:
                    # No new observations and sun is probably up!
                    print("No new images for %d minutes and UT = %02d:%02d > "
                            "15:00 so sun is probably up!" % 
                            (nnc, gm.tm_hour, gm.tm_min))
                    print "Time to hibernate until we have a new raw directory"
                    doit = False
                    # Normal termination
                    ret = True
                else:
                    print("No new image for %d minutes but UT = %02d:%02d < "
                            "16:00, so keep waiting" % (nnc, gm.tm_hour, gm.tm_min))
    # Handle a ctrl-C
    except KeyboardInterrupt:
       sys.exit("Exiting")

    return ret


def go():

    # Infinite loop
    dobs = True
    # Keep track of iterations
    its = 0
    # Get all raw directories
    rawlist = sorted([d for d in glob.glob('/scr2/sedm/raw/20??????') 
                        if os.path.isdir(d)])
    nraw = len(rawlist)
    try:
        while dobs:
            stat = ObsLoop(rawlist)
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
                    new_rawlist = sorted([d for d in 
                        glob.glob('/scr2/sedm/raw/20??????') 
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
    go()

