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


def find_line_match(lines, match_dict):
    """
    Returns the line number of first line matching the position-word pairs
    provided within matchDict, if present within lines.
    If no match is found, -1 is returned.
    """
    line_number = -1

    for index in range(0, len(lines)):
        for position_key in match_dict:
            if (len(lines[index].split()) == 0 or 
                    len(lines[index].split()) <= position_key or
                    lines[index].split()[position_key] !=
                    match_dict[position_key]):
                pass
            else:
                line_number = index
                break

    return line_number


def parse_and_fill(spec, snidoutput):
    """
    Receive the SEDM file and the output from snid. It parses the output
    from snid, fills a dictionary with the results and fills some of the
    comments in the original file, so that these could be used in the marshal.
    """
    
    pars = {"zmed": -1, "zmederr": -1, "agem": -1, "agemerr": -1,
            "Ia": 0, "Ib": 0, "Ic": 0, "II": 0, "NotSN": 0, "rlap": 0,
            "bestMatchType": "", "bestMatchSubtype": "", "bestMatchRedshift": 0}

    if os.path.exists(snidoutput):
        with open(snidoutput, "r") as snid:

            lines = snid.readlines()

            zmed_line = find_line_match(lines, {0: 'zmed'})
            agem_line = find_line_match(lines, {0: 'agem'})
            type_ia_line = find_line_match(lines, {0: 'Ia'})
            type_ib_line = find_line_match(lines, {0: 'Ib'})
            type_ic_line = find_line_match(lines, {0: 'Ic'})
            type_ii_line = find_line_match(lines, {0: 'II'})
            type_not_sn_line = find_line_match(lines, {0: 'NotSN'})
            best_match_line = find_line_match(lines, {0: '1'})

            pars["zmed"] = float(lines[zmed_line].split()[1])
            pars["zmederr"] = float(lines[zmed_line].split()[2])
            pars["agem"] = float(lines[agem_line].split()[1])
            pars["agemerr"] = float(lines[agem_line].split()[2])
            pars["Ia"] = float(lines[type_ia_line].split()[2])
            pars["Ib"] = float(lines[type_ib_line].split()[2])
            pars["Ic"] = float(lines[type_ic_line].split()[2])
            pars["II"] = float(lines[type_ii_line].split()[2])
            pars["NotSN"] = float(lines[type_not_sn_line].split()[2])
            pars["rlap"] = float(lines[best_match_line].split()[4])
            pars["bestMatchType"] = \
                lines[best_match_line].split()[2].split("-")[0]
            if len(lines[best_match_line].split()[2].split("-")) > 1:
                pars["bestMatchSubtype"] = \
                    lines[best_match_line].split()[2].split("-")[1]
            else:
                pars["bestMatchSubtype"] = '-'
            pars["bestMatchRedshift"] = float(lines[best_match_line].split()[5])

    print pars

    with open(spec, "r") as specIn:
        spec_lines = specIn.readlines()
        
        comment_lines_count = 0
        for line in spec_lines:
            if line.split()[0] == '#':
                comment_lines_count += 1
            else:
                break

        if os.path.exists(snidoutput):
            spec_lines.insert(comment_lines_count,
                              "# SNIDFRAC_NOTSN: " + str(pars["NotSN"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDFRAC_II: " + str(pars["II"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDFRAC_IC: " + str(pars["Ic"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDFRAC_IB: " + str(pars["Ib"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDFRAC_IA: " + str(pars["Ia"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDAGEMERR: " + str(pars["agemerr"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDAGEM: " + str(pars["agem"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDZMEDERR: " + str(pars["zmederr"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDZMED: " + str(pars["zmed"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDMATCHREDSHIFT: " +
                              str(pars["bestMatchRedshift"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDMATCHRLAP: " + str(pars["rlap"]) + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDMATCHSUBTYPE: " +
                              pars["bestMatchSubtype"] + "\n")
            spec_lines.insert(comment_lines_count,
                              "# SNIDMATCHTYPE: " +
                              pars["bestMatchType"] + "\n")
        else:
            spec_lines.insert(comment_lines_count,
                              "# SNIDMATCHTYPE: NONE\n")

    with open(spec, "w") as specOut:
        specOut.write("".join(spec_lines))


def run_snid(spec_dir='./', overwrite=False):
    """
    Runs snid in batch mode on all the PTF*_SEDM.txt files found in the given
    directory.  If a given file was already classified, it skips it, unless
    overwrite is requested.
    """
    
    for fl in glob.glob(os.path.join(spec_dir, "PTF*_SEDM.txt")):

        # retrieve the quality of the spectra.
        with open(fl, "r") as sfl:
            l = sfl.readlines()

            q = [li for li in l if "QUALITY" in li]

            if len(q) > 0:
                token = re.search(r'\(?([0-9]+)\)?', q[0])
                q = int(token.group(1))
            else:
                q = 5

            # If quality is good, check for previous classification
            if q < 3:
                clas = [li for li in l if "TYPE" in li]

                # If the file has been classified, move to the next
                if len(clas) > 0 and not overwrite:
                    continue
            else:
                continue
            
        # If we are here, we run the classification with snid
        cm = "snid wmin=4000 wmax=9500 skyclip=1 medlen=20 aband=1" \
             " rlapmin=4 inter=0 plot=0 %s" % fl
        print cm
        try:
            subprocess.call(cm, shell=True)
        except:
            print "Error running snid"
            continue
        
        snidoutput = fl.replace(".txt", "_snid.output")
        parse_and_fill(fl, snidoutput)
            
    
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
    run_snid(spec_dir=specdir)

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
            try:
                subprocess.call(cmd, shell=True)
            except:
                print "Error copying the file"
                pass
