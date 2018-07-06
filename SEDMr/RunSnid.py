# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 11:23:00 2016

@author: nblago
"""

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
            # Test type
            typstr = lines[best_match_line].split()[2]
            sntype = ''.join(c for c in typstr if ord(c) >= 32)
            if len(sntype) > 0:
                pars["bestMatchType"] = sntype.split("-")[0]
                if len(sntype.split("-")) > 1:
                    pars["bestMatchSubtype"] = sntype.split("-")[1]
                else:
                    pars["bestMatchSubtype"] = '-'
            else:
                pars["bestMatchType"] = "None"
                pars["bestMatchSubtype"] = '-'
            pars["bestMatchRedshift"] = float(lines[best_match_line].split()[5])

    print ("SNID RESULTS: Type=%(bestMatchType)s, Rlap=%(rlap).2f, "
           "Age=%(agem).2f+-%(agemerr)s day, "
           "z=%(zmed).4f+-%(zmederr).4f" % pars)

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

    return pars["bestMatchType"], pars


def run_snid(spec_file=None, overwrite=False):
    """
    Runs snid in batch mode on the input file.  If a given file was already
    classified, it skips it, unless overwrite is requested.
    """

    ran = False
    if spec_file is not None:
        fl = spec_file
        # retrieve the quality and classification of the spectra.
        with open(fl, "r") as sfl:
            l = sfl.readlines()

            q = [li for li in l if "QUALITY" in li]

            if len(q) > 0:
                token = re.search(r'\(?([0-9]+)\)?', q[0])
                q = int(token.group(1))
            else:
                q = 5
            print("quality = %d" % q)

            # check for previous classification
            clas = [li for li in l if "TYPE" in li]
            print("classification: ", clas)

        if q < 3 and (len(clas) <= 0 or overwrite):
            # If we are here, we run the classification with snid
            cm = "snid wmin=4500 wmax=9500 skyclip=1 medlen=20 aband=1" \
                 " rlapmin=4 inter=0 plot=2 %s" % fl
            print(cm)
            try:
                subprocess.call(cm, shell=True)
                ran = True
            except:
                print("Error running snid")
        else:
            if q >= 3:
                print("low quality spectrum")
            if len(clas) > 0:
                print("already classified")
    return ran
    # END run_snid


def record_snid(spec_file=None):
    """
    Records snid results in specfile. If a given file was already
    classified, it skips it, unless overwrite is requested.
    """

    if spec_file is not None:
        fl = spec_file
        try:
            snidoutput = fl.replace(".txt", "_snid.output")
            snid_type, pars = parse_and_fill(fl, snidoutput)
            psoutput = fl.replace(".txt", "_comp0001_snidflux.ps")
            if os.path.exists(psoutput):
                pngfile = fl.replace(".txt", "_" + snid_type + ".png")
                cm = "convert -flatten -rotate 90 " + psoutput + " " + pngfile
                subprocess.call(cm, shell=True)
            res = ("SNID RESULTS: Type=%(bestMatchType)s, Rlap=%(rlap).2f, "
                   "Age=%(agem).2f+-%(agemerr)s day, "
                   "z=%(zmed).4f+-%(zmederr).4f" % pars)
        except:
            print("Error recording snid")
            res = ""
        return res


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
        """

        Classifies an input SEDM ascii spectrum with SNID.

        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-d', '--specfile', type=str, dest="specfile",
                        help='*_SEDM.txt ascii spectrum file.',
                        default=None)
    parser.add_argument('--overwrite', action="store_true", default=False,
                        help='Overwrite existing classification')
    parser.add_argument('--update', action="store_true", default=False,
                        help='Update existing classification')

    args = parser.parse_args()

    # Run snid on extracted spectra and record results in specfile
    if args.update:
        print("Updating snid results in %s" % args.specfile)
        record_snid(spec_file=args.specfile)
    else:
        print("Running snid on and recording results in %s" % args.specfile)
        good = run_snid(spec_file=args.specfile, overwrite=args.overwrite)
        if good:
            record_snid(spec_file=args.specfile)
