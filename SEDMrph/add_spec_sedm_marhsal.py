#!/usr/local/EPD/epd-7.2.2/bin/python


import os
import pg
import ptfmarshal as ptf
import numpy as np
import argparse


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

    args = parser.parse_args()
    
    sedmfile = args.specfile
    print "Uploading SEDM file :", sedmfile
    
    try:
        db = dbconnect()
        upload_sedm(sedmfile, db)
    except IOError:
        print "Exception when uploading the file to the DB"
	pass