# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 11:23:00 2016

@author: nblago
"""
import datetime
import glob
import os
import argparse
import subprocess
import re

def parse_and_fill(spec, snidoutput):
    '''
    Receive the SEDM file and the output from snid. It parses the output from snid, fills a dictionary with the results
    and fills some of the comments in the original file, so that these could be used in the marshal.
    '''
    
    pars = {"zmed":-1, "zmederr":-1, "agem": -1, "agemerr": -1, "Ia":0, "Ib":0, "Ic":0, "II":0, "NotSN":0}
    
    with open(snidoutput, "r") as snid:
        
        lines = snid.readlines()
        pars["zmed"] = float(lines[35].split(" ")[1])
        pars["zmed"] = float(lines[35].split(" ")[2])
        pars["agem"] = float(lines[36].split(" ")[1])
        pars["agemerr"] = float(lines[36].split(" ")[2])
    
def run_snid(specdir):

    for f in glob.glob(os.path.join(specdir, "*.txt")):

        #retrieve the quality of the spectra.
        with open(f, "r") as sf:
            a = sf.readlines()
            
            clas = [ai for ai in a if "TYPE" in ai]
            
            #If the file has been classified, move to the next
            if len(clas) > 0:
                continue
            
        #Else, we run the classification with snid
        cmd = "snid wmin=4500 wmax=9500 skyclip=1 medlen=20 aband=1 rlapmin=4 plot=0 %s"%f                
        print cmd
        try:
            subprocess.call(cmd, shell=True)
        except:
            print "Error running snid"
            continue
        
        snidoutput = f.replace(".txt", "_snid.output")
        parse_and_fill(f, snidoutput)
            
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
        '''

        Copies the output ascii spectra to the iPTF marshal

        Checks to be sure spectra have quality > 3 and only
        copies spectra for objects with "PTF" prefix.

        Specify input directory with -d, or --specdir parameters.
        If none specified, use current date directory in /scr2/sedmdrp/redux/
            
        ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-d', '--specdir', type=str, dest="specdir",
                        help='Directory with output PTF*_SEDM.txt files.',
                        default=None)

    args = parser.parse_args()
    
    specdir = args.specdir
    print "Directory where reduced spectra reside: ", specdir
    
    if specdir is None:
        timestamp = datetime.datetime.isoformat(datetime.datetime.utcnow())
        timestamp = timestamp.split("T")[0].replace("-", "")
        specdir = os.path.join("/scr2/sedmdrp/redux/", timestamp)
    else:
	specdir = os.path.abspath(specdir)
	timestamp = os.path.basename(specdir)
    os.chdir(specdir)
    print os.getcwd()
    sedmfiles = glob.glob("PTF*.txt")
    print "Copying", sedmfiles
    
    for f in sedmfiles:
        
        qual = -1
        # retrieve the quality of the spectra.
        with open(f, "r") as sf:
            a = sf.readlines()
            
            qual = [ai for ai in a if "QUALITY" in ai]
            
            if len(qual) > 0:
                match = re.search(r'\(?([0-9]+)\)?', qual[0])
                qual = int(match.group(1))

        # Only write the spectra that have good qualities.
        if qual < 3:
            newname = os.path.basename(f).split("_")[0].replace("PTF", "")
            newname += "_%s_P60_v1.ascii" % timestamp
            cmd = "rcp %s nblago@yupana.caltech.edu:/scr/apache/htdocs/" \
                  "marshals/transient/ptf/spectra/sedm_to_upload/%s" % (f,
                                                                        newname)
            print cmd
            try:
                subprocess.call(cmd, shell=True)
            except:
                print "Error copying the file"
                pass
