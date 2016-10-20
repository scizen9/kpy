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


def find_line_match(lines, matchDict):
    """
    Returns the line number of first line matching the position-word pairs
    provided within matchDict, if present within lines.
    If no match is found, -1 is returned.
    """
    lineNumber = -1

    for index in range(0, len(lines)):
        for positionKey in matchDict:
            if (len(lines[index].split()) == 0 or 
                    len(lines[index].split()) <= positionKey or 
                    lines[index].split()[positionKey] != 
                    matchDict[positionKey]):
                pass
            else:
                lineNumber = index
                break

    return lineNumber


def parse_and_fill(spec, snidoutput):
    """
    Receive the SEDM file and the output from snid. It parses the output
    from snid, fills a dictionary with the results and fills some of the
    comments in the original file, so that these could be used in the marshal.
    """
    
    pars = {"zmed": -1, "zmederr": -1, "agem": -1, "agemerr": -1,
            "Ia": 0, "Ib": 0, "Ic": 0, "II": 0, "NotSN": 0, "rlap": 0,
            "bestMatchType": "", "bestMatchSubtype": "", "bestMatchRedshift": 0}

    if os.exists(snidoutput):
        with open(snidoutput, "r") as snid:

            lines = snid.readlines()

            zmedLine = find_line_match(lines, {0: 'zmed'})
            agemLine = find_line_match(lines, {0: 'agem'})
            typeIaLine = find_line_match(lines, {0: 'Ia'})
            typeIbLine = find_line_match(lines, {0: 'Ib'})
            typeIcLine = find_line_match(lines, {0: 'Ic'})
            typeIILine = find_line_match(lines, {0: 'II'})
            typeNotSNLine = find_line_match(lines, {0: 'NotSN'})
            bestMatchLine = find_line_match(lines, {0: '1'})

            pars["zmed"] = float(lines[zmedLine].split()[1])
            pars["zmederr"] = float(lines[zmedLine].split()[2])
            pars["agem"] = float(lines[agemLine].split()[1])
            pars["agemerr"] = float(lines[agemLine].split()[2])
            pars["Ia"] = float(lines[typeIaLine].split()[2])
            pars["Ib"] = float(lines[typeIbLine].split()[2])
            pars["Ic"] = float(lines[typeIcLine].split()[2])
            pars["II"] = float(lines[typeIILine].split()[2])
            pars["NotSN"] = float(lines[typeNotSNLine].split()[2])
            pars["rlap"] = float(lines[bestMatchLine].split()[4])
            pars["bestMatchType"] = \
                lines[bestMatchLine].split()[2].split("-")[0]
            pars["bestMatchSubtype"] = \
                lines[bestMatchLine].split()[2].split("-")[1]
            pars["bestMatchRedshift"] = float(lines[bestMatchLine].split()[5])

    print pars

    with open(spec, "r") as specIn:
        specLines = specIn.readlines()
        
        commentLinesCount = 0
        for line in specLines:
            if line.split()[0] == '#':
                commentLinesCount += 1
            else:
                break

        if os.exists(snidoutput):
            specLines.insert(commentLinesCount,
                             "# SNIDFRAC_NOTSN: " + str(pars["NotSN"]) + "\n")
            specLines.insert(commentLinesCount,
                             "# SNIDFRAC_II: " + str(pars["II"]) + "\n")
            specLines.insert(commentLinesCount,
                             "# SNIDFRAC_IC: " + str(pars["Ic"]) + "\n")
            specLines.insert(commentLinesCount,
                             "# SNIDFRAC_IB: " + str(pars["Ib"]) + "\n")
            specLines.insert(commentLinesCount,
                             "# SNIDFRAC_IA: " + str(pars["Ia"]) + "\n")
            specLines.insert(commentLinesCount,
                             "# SNIDAGEMERR: " + str(pars["agemerr"]) + "\n")
            specLines.insert(commentLinesCount,
                             "# SNIDAGEM: " + str(pars["agem"]) + "\n")
            specLines.insert(commentLinesCount,
                             "# SNIDZMEDERR: " + str(pars["zmederr"]) + "\n")
            specLines.insert(commentLinesCount,
                             "# SNIDZMED: " + str(pars["zmed"]) + "\n")
            specLines.insert(commentLinesCount,
                             "# SNIDMATCHREDSHIFT: " +
                             str(pars["bestMatchRedshift"]) + "\n")
            specLines.insert(commentLinesCount,
                             "# SNIDMATCHRLAP: " + str(pars["rlap"]) + "\n")
            specLines.insert(commentLinesCount,
                             "# SNIDMATCHSUBTYPE: " +
                             pars["bestMatchSubtype"] + "\n")
            specLines.insert(commentLinesCount,
                             "# SNIDMATCHTYPE: " + pars["bestMatchType"] + "\n")
        else:
            specLines.insert(commentLinesCount,
                             "# SNIDMATCHTYPE: NONE\n")

    with open(spec, "w") as specOut:
        specOut.write("".join(specLines))


def run_snid(specdir, overwrite=False):
    """
    Runs snid in batch mode on all the PTF*_SEDM.txt files found in the given
    directory.  If a given file was already classified, it skips it, unless
    overwrite is requested.
    """
    
    for f in glob.glob(os.path.join(specdir, "PTF*_SEDM.txt")):

        # retrieve the quality of the spectra.
        with open(f, "r") as sf:
            a = sf.readlines()

            qual = [ai for ai in a if "QUALITY" in ai]

            if len(qual) > 0:
                match = re.search(r'\(?([0-9]+)\)?', qual[0])
                qual = int(match.group(1))
            else:
                qual = 5

            if qual < 3:
                clas = [ai for ai in a if "TYPE" in ai]

                # If the file has been classified, move to the next
                if len(clas) > 0 and not overwrite:
                    continue
            else:
                continue
            
        # Else, we run the classification with snid
        cmd = "snid wmin=4000 wmax=9500 skyclip=1 medlen=20 aband=1 rlapmin=4 inter=0 plot=0 %s"%f
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
        """

        Copies the output ascii spectra to the iPTF marshal

        Checks to be sure spectra have quality > 3 and only
        copies spectra for objects with "PTF" prefix.

        Specify input directory with -d, or --specdir parameters.
        If none specified, use current date directory in /scr2/sedmdrp/redux/

        """, formatter_class=argparse.RawTextHelpFormatter)

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

    # Run snid on extracted spectra
    print "Running snid on PTF*_SEDM.txt files in %s" % specdir
    run_snid(specdir=specdir)

    sedmfiles = glob.glob("PTF*.txt")
    print "Copying", sedmfiles
    
    versions = {}
    
    for f in sedmfiles:
        
        qual = 5
        # retrieve the quality of the spectra.
        with open(f, "r") as sf:
            a = sf.readlines()
            
            qual = [ai for ai in a if "QUALITY" in ai]
            
            if len(qual) > 0:
                match = re.search(r'\(?([0-9]+)\)?', qual[0])
                qual = int(match.group(1))
            else:
                qual = 5

        # Only write the spectra that have good qualities.
        if qual < 3:
            newname = os.path.basename(f).split("_")[0].replace("PTF", "")
            
            version = 1
            if not versions.has_key(newname):
                versions[newname] = 1
            else:
                versions[newname] += 1
                version = versions[newname]
                
            newname += "_%s_P60_v%d.ascii" % (timestamp, version)
            cmd = "rcp %s nblago@yupana.caltech.edu:/scr/apache/htdocs/" \
                  "marshals/transient/ptf/spectra/sedm_to_upload/%s" % (f,
                                                                        newname)
            print cmd
            cmd = "ls %s" % f
            try:
                subprocess.call(cmd, shell=True)
            except:
                print "Error copying the file"
                pass
