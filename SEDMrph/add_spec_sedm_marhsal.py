#!/usr/local/EPD/epd-7.2.2/bin/python


import os
import pg
import ptfmarshal as ptf
import numpy as np
import argparse
import shutil
import re
import sys

def dbconnect():
    """
    Try to connect to the database
    return the connection if successful or none otherwise
    """

    try:
        db = pg.DB(host='localhost',dbname='ptftransient',user='ptftransient',passwd='discovery')
    except pg.InternalError:
        db = None

    return db



def get_sourceid(db, sourcename):
    query = "SELECT id FROM sources \
            WHERE name='%s';" % sourcename

    cursor = db.query(query)
    if cursor.ntuples():
        return cursor.dictresult()[0]['id']
    else:
        return None
        
def get_sedm_paramters(sedmfile):
    '''
        # TELESCOPE: P60
        # INSTRUMENT: SED-Machine
        # USER: sedmdrp
        # OBJECT: PTF16eqp
        # OUTFILE: PTF16eqp_SEDM.txt
        # OBSUTC: 2016 08 08 06:11:18.2
        # EXPTIME: 2700
        # SKYSUB: On
        # QUALITY: 2
        # AIRMASS: 1.13

    '''
    
    sedmp = {}
    
    if os.path.isfile(sedmfile):
        with open(sedmfile, "r") as f:
    
            for l in f.readlines():
                if (l.startswith("#")):
                    l = l.replace("# ","").replace(" ","")
                    vals = l.strip().split(":")
                    key = vals[0]
                    value = "".join(vals[1:])
                    sedmp[key] = value
        return  sedmp
    else:
        print "File %s foes not exist!"%sedmfile
        return None
    

def upload_sedm(sedmfile, db):
    
    ''' id | sourceid | programid | instrumentid |         obsdate         | exptime | snr | snrwav | minwav  | maxwav  | 
    dataformat |                   datapath                   |      observer       | reducedby | classification |     lastmodified        | spectype | uploadedby '''
    
    
    
    #Retrieve information from the SEDM file header.
    sedmp = get_sedm_paramters(sedmfile)
       
    #Read the spectrum
    spec = np.genfromtxt(sedmfile, comments="#")

    #Fill the parameters for the DB with the right values.
    params = {}
    params['id'] = -1
    params['sourceid'] =  get_sourceid(db, sedmp["OBJECT"].replace("PTF", ""))
    params['programid'] = 4
    params['instrumentid'] = 65
    params['obsdate'] = sedmp["OBSUTC"][0:8]
    params['exptime'] = sedmp["EXPTIME"]
    params['snr'] = ''
    params['snrwav'] = ''
    params['minwav'] = np.nanmin(spec[:,0])
    params['maxwav'] = np.nanmax(spec[:,0])
    params['dataformat'] = 'ascii'
    params['datapath'] = os.path.join('ptf/spectra/data/', os.path.basename(sedmfile))
    params['observer'] = 'SEDM robot'
    params['reducedby'] = 'SEDM reducer'
    params['classification'] = ''
    params['spectype'] = 'object'
    params['uploadedby'] = 'bookkeeper'

    
    # commit the entry to the DB
    newrow = ptf.altertable(db, 'spec', -1, params)
    
    print newrow
            





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''

        Runs astrometry.net on the image specified as a parameter and returns 
        the offset needed to be applied in order to center the object coordinates 
        in the reference pixel.
            
        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('-f', '--sedmfile', type=str, dest="specfile", help='SEDMfile in ascii format.', default=None)
    parser.add_argument('-c', '--copytodest', action='store_true', dest="copy", default=True, help='Whether to copy this file to its marshal destination.')
    parser.add_argument('-o', '--overwrite', action='store_true', dest="overwrite", default=False, help='Whether to overwrite the current file with a newer version.')


    args = parser.parse_args()
    
    sedmfile = args.specfile
    copy = args.copy
    overwrite = args.overwrite
    
    print "Uploading SEDM file :", sedmfile
    
    
    #We copy the file to its final marshal destination.
    #True by default.
    if (copy):
        sedmfiledest = os.path.join("/scr/apache/htdocs/marshals/transient/ptf/spectra/data", os.path.basename(sedmfile))
        
        #If the file does not exist in destination, we copy it
        if (not os.path.isfile(sedmfiledest)):
            print "File %s does not exist in destination as %s. Moving it there and updating the marshal."%(sedmfile, sedmfiledest)
            shutil.move(sedmfile, sedmfiledest)
        #If the file does exist
        else:
            #If it is said to overwrite, we create a new version.
            if (overwrite):
                #Retrieve the current version and sum one
                match = re.search(r'_v\(?([0-9]+)\)?', os.path.basename(sedmfile))
                version = int(match.group(1))
                sedmfiledest = os.path.join(os.path.dirname(sedmfiledest), os.path.basename(sedmfile).replace("_v%d."%version, "_v%d."%(version+1)))
                print "Moving file %s to %s"%(sedmfile, sedmfiledest)
                shutil.move(sedmfile, sedmfiledest)
            #If the file already exists and we don't overwrite, we just delete it from the reception folter.
            #Generally this means that the spectrum is already in the marshal.
            else:
                print "File %s already exists in destination as %s. No overwrite selected. Skipping."%(sedmfile, sedmfiledest)
                os.remove(sedmfile)
                sys.exit(0)                
    else:
        sedmfiledest = sedmfile
        
    try:
        db = dbconnect()
        upload_sedm(sedmfiledest, db)
    except IOError:
        print "Exception when uploading the file to the DB"
	pass