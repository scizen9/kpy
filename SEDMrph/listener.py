# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 22:22:15 2016

@author: nadiablago
"""

import matplotlib
matplotlib.use("Agg")
import os, time
import recenter_ifu
import socket
import sextractor
import sys
import sao
import subprocess
import datetime
import logging

def start_listening_loop():
    '''
    Start accepting connections from pylos.
    '''
        
    #Set address
    ip = "pharos.caltech.edu"
    port = 5006
    
    #Log into a file
    FORMAT = '%(asctime)-15s %(levelname)s [%(name)s] %(message)s'
    root_dir = "/tmp/"
    timestamp=datetime.datetime.isoformat(datetime.datetime.utcnow())
    timestamp=timestamp.split("T")[0]
    logging.basicConfig(format=FORMAT, filename=os.path.join(root_dir, "listener_{0}.log".format(timestamp)), level=logging.INFO)
    logger = logging.getLogger('listener')
    
    #bind socket
    s = socket.socket(socket.AF_INET,socket.SOCK_STREAM)
    try:
        s.bind((ip,port))
    except socket.error:
        return
    s.listen(10)
    
    #create continous while loop to listen for request
    while True:
        connection,caddress = s.accept()
        time.sleep(1)
        cmd = "touch /tmp/sedm_listener_alive"
        subprocess.call(cmd, shell=True)

        while True:
            cmd = "touch /tmp/sedm_listener_alive"
            subprocess.call(cmd, shell=True)
            data = connection.recv(2048)
            logger.info( "Incoming command: %s " % data)
            
            command = data.split(",")[0]
            
            if "FOCUS" in command:
                logger.info( "Finding the best focus.")
                lfiles = data.split(",")[1:]
                for i, lf in enumerate(lfiles):
                    lf = lf.replace("raw", "phot")
                    lfiles[i] = lf
                #Wait until the image is available.                
                while (not os.path.isfile(lfiles[-1])):
                        time.sleep(0.5)
                focus, sigma = sextractor.get_focus(lfiles, plot=False)
                connection.sendall("%.2f,%.2f\n"%(focus, sigma))
                logger.info("Selected focus: %.2f,%.2f\n"%(focus, sigma))
            elif "SAO" in command:
                try:
                    logger.info( "Looking for a nice SAO star.")
                    name, ra, dec = sao.get_sao()
                    connection.sendall("%d,%s,%s,%s\n"%(0,name, ra, dec))
                    logger.info("Found star. Returning: %d,%s,%s,%s\n"%(0,name, ra, dec))
                except Exception as e:
                    with open("/tmp/e", "w") as f:
                        f.write(str(sys.exc_info()[0]))
                    logger.error(str(sys.exc_info()[0]))
                    logger.error(e)
                    connection.sendall("%d,%s,%s,%s\n"%(-1,"null", (datetime.datetime.utcnow()[4]+3)*15, 40))

            elif "STATS" in command:
                logger.info( "Finding the statistics for the image.")
                lfile = data.split(",")[1]
                lfile = lfile.strip()
                lfile = lfile.replace("raw", "phot")
                while (not os.path.isfile(lfile)):
                        time.sleep(0.5)
                try:
                    nsources, fwhm, ellipticity = sextractor.get_image_pars(lfile)
                    connection.sendall("%d,%d,%d,%.3f\n"%(0, nsources, fwhm, ellipticity))
                except Exception as e:
                    with open("/tmp/e", "w") as f:
                        f.write(str(sys.exc_info()[0]))
                    logger.error(str(sys.exc_info()[0]))
                    logger.error(e)
                    connection.sendall("%d,%d,%d,%.3f\n"%(-1, 0, 0, 0))

            elif "OFFSET" in command:
                try:
                    isABstr = data.split(",")[1]
                    isAB = (isABstr =="AB")
                    image = data.split(",")[2].rstrip()
                    logger.info("Get Offsets AB=%s for image %s"%(isAB,image))
                    
                    astrofile = os.path.basename(image)
                    date = astrofile.split("_")[0].replace("rc","")
                    astrofile = astrofile.replace("rc", "a_rc").replace(".new", ".fits")
                    endpath = "/scr2/sedm/phot/%s/%s"%(date, astrofile)
                    if (not os.path.isdir(os.path.dirname(endpath))):
                        os.makedirs(os.path.dirname(endpath))
                    os.system('scp -i /home/sedm/.ssh/guider_rsa developer@p200-guider.palomar.caltech.edu:%s %s'%(image, endpath))
                    astro = False
                    #Only run astrometry if image is unavailable.                
                    if (not os.path.isfile(endpath)):
                        astro=True
                        endpath = endpath.replace(".new", ".fits")
                        logger.error("Astrometry resolved image %s could not be copied into %s. Setting astrometry to True."%(image, astrofile))
                    res = recenter_ifu.main(endpath, isAB, astro=astro, plot=True)
                    retcode = res[0]
                    offsets = res
                    logger.info( "OFFSETS %s  Return code %d"%(offsets,retcode))
                    if(isAB):
                        connection.sendall("%d,%s,%s,%s,%s\n"%(retcode,offsets[1], offsets[2], offsets[3], offsets[4]))
                    else:
                        connection.sendall("%d,%s,%s\n"%(retcode,offsets[1], offsets[2]))
                except Exception as e:
                    logger.error( "Error occurred when processing command  " + str(data))
                    with open("/tmp/offsets_error", "w") as f:
                        f.write(str(sys.exc_info()[0]))
                        f.write(str(e))
                    logger.error(str(sys.exc_info()[0]))
                    logger.error(e)
                    connection.sendall("%d,%s,%s\n"%(-1,0,0))
            else:
                break
            
if __name__ == '__main__':
    '''
    If the program has not been running since 1 minute, then the main funciton relaunches the main loop.
    Otherwise it does nothing.
    '''
    
    #If the file was modified last time more than 60s ago, relaunch the listener.   
    modified = datetime.datetime.strptime(time.ctime(os.path.getmtime("/tmp/sedm_listener_alive")), "%a %b  %d %H:%M:%S %Y") 
    now = datetime.datetime.now()
    
    if ( (now - modified).seconds > 60):
        start_listening_loop()

                
