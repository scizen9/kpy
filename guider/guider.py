# TELESCOPE GUIDER
# AUTHOR: Donal O Sullivan (dosulliv@caltech.edu)
# DATE: 24/04/2014
# Modified Mar 2015 by Nick Konidaris to not use SFTP

# DESCRIPTION
#   This code uses Telnet and local disk to communicate with the SED Machine on the Palomar 60'' Telescope.
#   Its purpose is to constantly download and analyze images for tracking purposes, sending pointing
#   adjustments if the telescope drifts.
#   All connection/locally-dependent variables are extracted outside the method 'auto()' so as to
#   avoid any hard-coded parameters within the method. As such, you should only ever have to change
#   these input parameters to use it on a different platform or even for a different instrument.

#   Most packages used are default in most installations of Python, except perhaps the following:
#   pyfits - fits image manipulator
#   pysftp - secure SSH file transfer protocol interface
#   telnetlib - telnet interface

#   The package 'sedmtools' is written by myself and should be included in the same directory. It
#   contains a few functions I wrote such as star-finding algorithms, offset measurements and
#   focus-finding. 

#   The directory defined by the parameter 'handled_dir' will be where files are downloaded to and 
#   where local copies of logs will be kept, although logs will also be uploaded to the server.
#   Make sure user has read/write access to this directory.

#   As far as I know this code handles itself well enough in terms of errors, but if it does bug
#   out, no information should be lost. Logs are updated and saved to the server and locally after 
#   each file download/check.

#1. SETUP

#1.1 IMPORTS
import os
import sys
import telnetlib
import sedmtools
import shutil
import numpy as np
import pyfits as pf
from time import sleep
from time import localtime

# 1.2 CODE PARAMETERS
params = {}
params["use_telnet"] = True # Boolean: True - send guider commands to Telescope; False - run code but do not send commands
params["telnet_ip"] = "198.202.125.194" # String: IP Address for Telnet commands
params["telnet_port"] = 49302 # Integer: Port for Telnet commands
params["handled_dir"] = "/home/sedm/guider/auto2/" # Directory to store handled files and logs
params["debug"] = False # Set to true to activate debugging mode output


def touch(path):
    """touch a file name at path"""

    f = open(path, "w")
    f.close()

# END 1. SETUP


def auto(params): # Run this method to begin the guider

    # 2 GLOBAL VARS
    global _obName,_obRA,_obDec,_dRA_0,_dDec_0,_ra_off,_dec_off,_star,_starBox,_sBs,\
    _log, _stored,_months,_host,_user,_pass,_handled_dir,_disconnect,_t,_use_telnet

    _debug = params["debug"] # Activates/de-activates certain print statements for debugging

    # 2.1 Telnet for sending commands to guider
    _use_telnet = params["use_telnet"]
    if _use_telnet:
        _t = telnetlib.Telnet(params["telnet_ip"], params["telnet_port"])
        if _debug:
            print "Connecting to telnet"

    if _debug:
        print "Creating global variables"

    # 2.2 Setup directories
    _handled_dir = params["handled_dir"]
    if not os.path.exists(_handled_dir):
        raise Exception("Path %s does not exist" % _handled_dir)
    _disconnect = False 

    # 2.3 Misc Variables
    _months = {1:"01",2:"02",3:"03",4:"04",5:"05",6:"06",7:"07",
                8:"08",9:"09",10:"10",11:"11",12:"12"}
    _stored = {}
    _obName = ""  # Object name
    _obRA = ""  # Object RA coordinate
    _obDec = ""  # Object DEC coordinate
    _dRA_0,_dDec_0 = 0, 0  # Calibration variables for offset, see below
    _ra_off, _dec_off = 0, 0  # Offset of actual image centre from target coordinates
    _star = (-1, -1)  # Star position in image data. Default (-1,-1) serves as error flag
    _starBox = []  # Will store 2D data around star
    _sBs = 20  # Size (pixels) of box to take around guide star
    _log = ""

    # END 2. GLOBAL VARS

    # 3. METHODS
    if _debug:
        print "Method defintiions"

    # 3.1 Log data into a file
    def log(data, filePath):
        log_file = open(filePath, 'a')
        log_file.write(data)
        log_file.close()

    # 3.2 Dual output: log & stdout
    def output(string, filePath):
        print(string),
        log(string, filePath)

    # 3.3 Get list of new files (that contain a given substring, e.g. '.fits' or 'rc') from server
    def getNew(handleddir,remotedir,substring):
        have = {}
        new_files = []
        handled_files = os.listdir(handleddir)
        try:
            remote_files = os.listdir(remotedir)
        except:
            print remotedir
            output("Error reading file list.",_log)
            return []

        for x in handled_files:
            have[x] = True

        for x in remote_files:
            if have.has_key(x):
                ignore = True
            else:
                ignore = False
            if not ignore and substring in x:
                new_files.append(x)

        return new_files

    # 3.4 Get whether the file is a daytime image, a new object or a new offset
    def getStatus(path):

        global _obName,_obRA, _obDec, _ra_off, _dec_off
        # print "......Status: %s" % path
        try:
            fitsImage = pf.open(path)
            status = {"daytime":False, "new obj":False, "moved":False}
        except:
            # wait 5s
            sleep(5)
            try:
                fitsImage = pf.open(path)
                status = {"daytime":False, "new obj":False, "moved":False}
            except:
                status = {"daytime":False, "new obj":False, "moved":False}
                return status
        try: 
            time = fileName[fileName.index('_')+1:fileName.index('.')]  # Extract time from filename
            nTime = float(time[0:2])*100 + float(time[3:5])  # Convert to a number: '20:01' -> 2001
        except: output("%20s Error parsing time from filename." % time, _log)
        try:
            # Extract Header Info
            obRA = fitsImage[0].header['OBJRA']
            obDec = fitsImage[0].header['OBJDEC']
            RAoff = fitsImage[0].header['RA_OFF']
            Decoff = fitsImage[0].header['DEC_OFF']
            obName = fitsImage[0].header['OBJECT']

        except:
            print "ERROR"
            output("%20s header error.\n" % time, _log)
            status["moved"] = False
            return status

        # Is image from before sunset/after sunrise? (e.g. flats)
        # if 530 < nTime < 2030: status["daytime"] = True

        # Are target coordinates different to current global values?
        if obRA!=_obRA or obDec!=_obDec or obName!=_obName:

            status["new obj"] = True
            _obName = obName
            _obRA = obRA
            _obDec = obDec
            _ra_off = RAoff
            _dec_off = Decoff

        # If not a new object, has the offset changed by more than 4.5'' in either direction?
        elif max( abs(_ra_off - RAoff), abs(_dec_off - Decoff) ) > 3.5:

            status["moved"] = True
            _obName = obName
            _ra_off = RAoff
            _dec_off = Decoff
        elif "(NS)" in obName:
            status["moved"] = True
            # Object is an asteroid or other small body

        return status

    # 3.5 Formatted output for coordinates, offsets and seeing conditions
    def outputCoords(fileName,path):
        global _obName,_obRA,_obDec,_star,_ra_off,_dec_off, _log
        output("\n%s" % (_obName) + 
            "\nCoordinates: %s %s %s" % (_obRA,_obDec,_star) +
            "\nOffsets: %.2f %.2f\n\n" % (_ra_off,_dec_off) +
            "%20s %10s %10s %10s %10s %10s %10s %10s\n" % 
                (fileName[0:fileName.index('_')+1],"dRA","dDec","sent",
                "fwhm_x","fwhm_y","fwhm_x:y","fwhm_avg"),
                path)

    # END 3. METHOD DEFINITIONS

    # 4. MAIN

    if _debug: 
        print "Beginning main loop"
        print "Disconnect = %s" % str(_disconnect)

    # 4.2 START MAIN LOOP
    while not _disconnect:

        if _debug: print "..Creating local variables"

        # 4.3 GET LOCAL VARIABLES: FILES, TIME, DATE, BOOLEANS

        # 4.3.1 COMMENT/UNCOMMENT TO SWITCH TO MANUAL DATE SELECTION
        # date = raw_input("Date: ")
        # date = date.split(' ')
        # year,month,day = int(date[0]),int(date[1]),int(date[2])

        # 4.3.2 COMMENT/UNCOMMENT TO SWITCH TO AUTOMATIC DATE SELECTION
        now = localtime()
        year,month,day = now.tm_year, now.tm_mon, now.tm_mday

        # 4.3.3 PARSE DATE INFORMATION
        day_s = str(day) if day >=10 else '0' + str(day)
        month_s = str(month) if month >=10 else '0'+str(month)
        monthName = _months[month]

        # 4.3.4 CREATE USEFUL STRINGS
        remote_dir = os.path.join("/home/sedm/sedm_data", str(year)+monthName+day_s) #Directory containing today's files on server
        _log = "%s_%s_%s_log.txt" % (str(year),month_s, day_s) #Log file for general program output
        seeing_log = "%s_%s_%s_seeing.txt" % (str(year),month_s, day_s) #Log file for seeing-condition data

        if _debug: print "..Getting new file list"

        # 4.4 GET NEW FILE LIST
        new_files = getNew(_handled_dir, remote_dir,"rc")

        if _debug: print ".... %s" % new_files
        # 4.5 DECIDE WHETHER TO GUIDE TELESCOPE OR NOT
        guide_telescope = True if (len(new_files)<=1 and _use_telnet) else False #Don't guide until up-to-date (i.e. <=1 new file)

        if _debug: print "..Beginning file loop" 
        sleep(5)
        # 4.6 RUN THROUGH NEW FILES
        for fileName in new_files:

            if _debug: print "....Creating local vars(2)"

            # 4.6.1 GET LOCAL VARIABLES
            handled_path = os.path.join(_handled_dir,fileName) #local path for file
            remote_path = os.path.join(remote_dir,fileName) #remote path for file
            time = fileName[fileName.index('_')+1:fileName.index('.')] #time of image
            dRA,dDec = "-","-" #Change in RA/Dec
            shift_sent = False #Boolean: move command sent to telescope?

            # 4.6.2 COPY FILES TO PROCESSED DIRECTORY
            touch(handled_path)

            # 4.6.3 GET 'STATUS' OF FILE (see method definition for details)
            if _debug: print "....Getting Status"
            status = getStatus(remote_path)
            image = pf.open(remote_path)

            # CAN'T REMEMBER WHY THIS IS HERE; DID WE ONLY WANT 60s EXPOSURES?
            if image[0].header['EXPTIME']!=30.0: continue

            if _debug: print "....Checking Status"

            # 4.6.4 Check status of object.

            # Ignore all day-time images
            if status["daytime"]: 

                if _debug: print "......Status = Day-time" 
                _stored[fileName] = True
                continue

            # If new object, moved, or no star found last time, then update star position
            elif status["new obj"] or status["moved"] or _star==(-1,-1):

                if _debug: 
                    print "......Status = Object/Telescope moved"
                    print "......Finding star"

                _star = sedmtools.findstar(image,100) #Get new position for star, returns (-1,-1) if not found

                if _star == (-1,-1):

                    _stored[fileName] = True
                    output("%20s: No star found in image.\n" % (time), _handled_dir+_log)
                    continue

                if _debug: print "......Updating zero-shift" 

                # Get (sBs x sBs) sized image of star at current position (deep copy from fits data)
                _starBox = np.empty_like(image[0].data[_star[0]-_sBs:_star[0]+_sBs,
                                        _star[1]-_sBs:_star[1]+_sBs])
                _starBox[:] = image[0].data[_star[0]-_sBs:_star[0]+_sBs,
                                    _star[1]-_sBs:_star[1]+_sBs]

                # Calibrate offset calculations (i.e. get what the method returns for two perfectly aligned images)
                _dRA_0,_dDec_0 = sedmtools.getshift(_starBox,_starBox)

                # Output coordinates
                if "(NS)" not in _obName: outputCoords(fileName,_handled_dir+_log)

            else: # If not moved and not a new object

                if _debug: print "......Status = None, calculating shift."

                # Get x,y position of star in previous image
                x,y = _star 

                # Get shift in the star's position by checking new image against old _starBox
                dRA, dDec = sedmtools.getshift(_starBox, 
                            image[0].data[x-_sBs: x+_sBs, y-_sBs:y+_sBs])

                # Apply offset calibrations
                dRA -= _dRA_0
                dDec -= _dDec_0

                # Try to telnet the command to adjust if we are currently guiding
                if guide_telescope:
                    if (dRA < 4) or (dDec < 4):
                        if _debug: print "......Sending Telnet shift to telescope" 

                        try:
                            _t.write("GM %s %s 10 10\n" % (dRA, dDec))
                            r = int(_t.expect(["-?\d"], 60)[2])
                            if r==0: shift_sent = True
                            else: 
                                if r != -3:
                                    output("Shift not sent. Error = %i\n" % r , _handled_dir+_log)
                        except:
                            shift_sent = False
                            output("Error sending shift via telnet.", _handled_dir+_log)

            # 4.6.5 GET SEEING CONDITIONS
            if _debug: print "....Getting FWHM"

            fwhm = sedmtools.getfwhm(image,_star) #Get FWHM of star

            # If FWHM measurement didn't work
            if fwhm == -1:

                if _debug: print "....FWHM Failed, getting new star"

                # Try to get new star-coordinates
                _star = sedmtools.findstar(image,100)

                # If starfind failed, output error, store and move on
                if _star == (-1,-1):
                    if _debug: print "....new star failed" 
                    output("%20s no usable star found." % time, _handled_dir+_log)
                    _stored[fileName] = True
                    continue

                # If starfind succeeded, output and try fwhm again
                else:
                    if _debug: print "....Getting FWHM (2)" 
                    outputCoords(fileName,_handled_dir+_log)  
                    fwhm = sedmtools.getfwhm(image,_star)

                    # If still failed, output error, store and move on
                    if fwhm == -1:
                        output("%20s FWHM measurement failed." % time, _handled_dir+_log)
                        _stored[fileName] = True
                        continue

            # 4.6.6 OUTPUT
            if _debug: print "....Output stage"

            # Normal results output stage
            f_x = "%.2f" % fwhm[0]
            f_y = "%.2f" % fwhm[1]
            f_r = "%.2f" % fwhm[2]
            f_a = "%.2f" % fwhm[3]
            if type(dRA)!=str: dRA = "%.2f" % dRA
            if type(dDec)!=str: dDec = "%.2f" % dDec
            output("%20s %10s %10s %10s %10s %10s %10s %10s\n" 
            % (time,dRA,dDec,shift_sent,f_x,f_y,f_r,f_a), _handled_dir+_log)
            sys.stdout.flush()
            _stored[fileName] = True
            f = open(_handled_dir+seeing_log,'a')
            f.write("%f-%s\n" % (fwhm[3],time))
            f.close()

            # Upload logs
            shutil.copyfile(_handled_dir+_log,remote_dir+_log)
            shutil.copyfile(_handled_dir+seeing_log,remote_dir+seeing_log) 

        if _debug: print "....Sleeping"

        sleep(3) # Wait 3 seconds before rerunning

auto(params)



