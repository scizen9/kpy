"""Generate Makefile for reducing ifu data"""

import os
import sys

import astropy.io.fits as pf

import NPK.Bar as Bar
import NPK.Standards as Stds


def extract_info(infiles):

    headers = []

    print("-- Ingesting headers --")
    update_rate = int(len(infiles) / (Bar.setup() - 1))
    if update_rate <= 0:
        update_rate = 1
    for ix, ifile in enumerate(infiles):
        if ix % update_rate == 0:
            Bar.update()
        FF = pf.open(ifile)
        FF[0].header['filename'] = ifile
        if 'JD' not in FF[0].header:
            # print("Skipping %s" % ifile)
            continue
        headers.append(FF[0].header)
        FF.close()

    Bar.done()

    return sorted(headers, key=lambda x: x['JD'])


def identify_observations(headers):
    """Return a list of object name, observation number, and list of files.

    e.g. will return:

    {'STD-BD+25d4655': {1: ['...']}, {2: ['...']}, 
           'PTF14dvo': {1: ['...', '...']}}
    
    where STD-BD+25d4655 was observed at the beginning and end of night. SN
    14dov was observed once with A-B.
    """
    JD = 0.

    objcnt = {}
    objs = {}
    calibs = {}

    for header in headers:
        if header['JD'] < JD:
            raise Exception("Headers not sorted by JD")
        JD = header['JD']

        fname = header['filename']
        obj = header['OBJECT']
        name = header['NAME']
        exptime = header['exptime']
        adcspeed = header['ADCSPEED']
        if "test" in obj:
            continue
        if "Calib" in obj or "bias" in obj:

            def appendToCalibs(Str):

                if Str in obj:
                    if "bias" in Str and exptime == 0.:
                        Str = "%s%1.1f" % (Str, adcspeed)
                        prefix = ""
                        suffix = ""
                    elif "Xe" in Str or "Hg" in Str or "Cd" in Str or \
                                    "Ne" in Str or "dome" in Str:
                        prefix = "b_"
                        suffix = ""
                    else:
                        prefix = "crr_b_"
                        suffix = ""

                    if "bias" in Str and exptime != 0.:
                        print("Mis-labeled bias with exptime > 0: %9.1f" %
                              exptime)
                    else:
                        calibs[Str] = calibs.get(Str, [])
                        calibs[Str].append(prefix + fname + suffix)

            appendToCalibs("bias")
            appendToCalibs("dome")
            appendToCalibs("Xe")
            appendToCalibs("Hg")
            appendToCalibs("Cd")
            appendToCalibs("Ne")
            appendToCalibs("twilight")

        if "Focus:" in obj:
            continue
        if "dark" in obj:
            continue
        if "Calib" in obj:
            continue
        if "STOW" in name:
            continue
        if obj.rstrip() == "":
            continue
        name = name.replace(" ", "_")
        name = name.replace(")", "_")
        name = name.replace("(", "_")
        name = name.replace("[", "_")
        name = name.replace("]", "_")
        name = name.replace("/", "_")
        name = name.replace(":", "_")

        # The 'A' position defines the start of an object set
        if '[A]' in obj or name not in objcnt:
            cnt = objcnt.get(name, 0) + 1
            vals = objs.get(name, {})
            objcnt[name] = cnt
            vals[cnt] = [fname]
            objs[name] = vals
        else:
            cnt = objcnt[name]
            objs[name][cnt].append(fname)

    print("\n-- Calibrations --")
    for k, v in calibs.items():
        print("%15s : %2.0i" % (k, len(v)))

    print("\n-- Standard Star Sets --")
    for k, v in objs.items():
        if "STD-" in k:
            print("%20s : %2.0i" % (k, len(v)))

    print("\n-- Science Object Sets --")
    for k, v in objs.items():
        if "STD-" not in k:
            print("%20s : %2.0i" % (k, len(v)))

    return objs, calibs


make_preamble = """
PY = ~/spy
PYC = ~/kpy/SEDMr
PYP = ~/kpy/SEDMrph
EXTSINGLE =  $(PY) $(PYC)/Extractor.py
ATM =  $(PY) $(PYC)/AtmCorr.py
EXTPAIR =  $(PY) $(PYC)/Extractor.py
FLEXCMD = $(PY) $(PYC)/Flexure.py
IMCOMBINE = $(PY) $(PYC)/Imcombine.py
PLOT = $(PY) $(PYC)/Check.py
REPORT = $(PY) $(PYC)/DrpReport.py
SPCCPY = $(PY) $(PYP)/sedmspeccopy.py
PTFREPORT = $(PY) $(PYC)/PtfDrpReport.py

BSUB = $(PY) $(PYC)/Debias.py
BGDSUB =  $(PY) $(PYC)/SubtractBackground.py
CRRSUB =  $(PY) $(PYC)/CosmicX.py

SRCS = $(wildcard ifu*fits)
BIAS = $(addprefix b_,$(SRCS))
CRRS = $(addprefix crr_,$(BIAS))
BACK = $(addsuffix .gz,$(addprefix bs_,$(CRRS)))
EXTR = $(subst .fits.gz,.npy,$(BACK))
FLEX = $(subst .fits,.npy,$(addprefix flex_,$(BACK)))
FIGS = $(subst .npy,_SEDM.pdf,$(EXTR))

mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
current_dir := $(notdir $(patsubst %/,%,$(dir $(mkfile_path))))

crr_b_% : b_%
	$(CRRSUB) --niter 4 --sepmed --gain 1.0 --readnoise 5.0 --objlim 1.8 \\
		--sigclip 8.0 --fsmode convolve --psfmodel gaussy --psffwhm=2 \\
		$< $@ mask$@

bs_crr_b_%.gz : crr_b_%
	$(BGDSUB) fine.npy $< --gausswidth=100

flex_bs_crr_b_%.npy : bs_crr_b_%.fits.gz
	$(FLEXCMD) cube.npy $< --outfile $@

bs_crr_b_%.npy : bs_crr_b_%.fits.gz flex_bs_crr_b_%.npy
	$(EXTSINGLE) cube.npy --A $< --outname $@ --flat_correction flat-dome-700to900.npy --Aoffset flex_$@ --specExtract

%_SEDM.pdf : sp_%.npy
	$(PLOT) --spec $< --savefig --interact

.PHONY: cleanstds newstds report ptfreport finalreport

bias: bias0.1.fits bias2.0.fits $(BIAS)
crrs: $(CRRS) 
back: $(BACK)
extr: $(EXTR)
figs: $(FIGS)

$(BIAS): bias0.1.fits bias2.0.fits
	$(BSUB) $(subst b_,,$@)

$(CRRS): 
	$(CRRSUB) --niter 4 --sepmed --gain 1.0 --readnoise 5.0 --objlim 1.8 \\
		--sigclip 8.0 --fsmode convolve --psfmodel gaussy --psffwhm=2 \\
		$(subst crr_,,$@) $@ mask$@

$(BACK): 
	$(BGDSUB) fine.npy $(subst .gz,,$(subst bs_,,$@)) --gausswidth=100


seg_dome.fits: dome.fits
	$(PY) $(PYC)/SexLamps.py dome.fits

seg_Hg.fits: Hg.fits
	$(PY) $(PYC)/SexSpectra.py Hg.fits

dome.fits_segments.npy: seg_dome.fits
	$(PY) $(PYC)/FindSpectra.py seg_dome.fits dome.fits dome.fits_segments --order 1

rough.npy: dome.fits_segments.npy seg_Hg.fits
	$(PY) $(PYC)/Wavelength.py rough --hgfits Hg.fits --hgcat cat_Hg.fits.txt --dome dome.fits_segments.npy --outname rough 

fine.npy: rough.npy Cd.fits Xe.fits
	$(PY) $(PYC)/Wavelength.py fine --cdfits Cd.fits --xefits Xe.fits --hgfits Hg.fits --hgassoc assoc_Hg.npy --outname fine

cube.npy: fine.npy
	$(PY) $(PYC)/Cube.py fine.npy --step make --outname cube.npy
	$(PLOT) --cube cube.npy --savefig
	$(PLOT) --cube cube.npy --lambdarms --savefig

bs_twilight.fits.gz: twilight.fits fine.npy
	$(BGDSUB) fine.npy twilight.fits --gausswidth=100

bs_dome.fits.gz: dome.fits fine.npy
	$(BGDSUB) fine.npy dome.fits --gausswidth=100

dome.npy: cube.npy dome.fits
	$(PY) $(PYC)/Extractor.py cube.npy --A dome.fits --outname dome --extflat

flat-dome-700to900.npy: dome.npy
	$(PY) $(PYC)/Flat.py dome.npy
    
wave: fine.npy
cube: cube.npy

flex: back $(FLEX)

$(FLEX): cube.npy
	$(eval OUTNAME = $(subst .gz,,$@))
	$(FLEXCMD) cube.npy $(subst flex_,,$(subst npy,fits,$@)) --outfile $(OUTNAME)

stds: flat-dome-700to900.npy std-correction.npy

cleanstds:
	rm -f std-correction.npy Standard_Correction.pdf

newstds: cleanstds stds

report:
	$(REPORT) | tee report.txt

upload:
	$(SPCCPY) --specdir $(dir $(mkfile_path))

ptfreport: upload
	$(PTFREPORT) 

finalreport: ptfreport
	$(REPORT) | tee report.txt | mail -s "SEDM DRP Report for $(current_dir)" neill@srl.caltech.edu,rsw@astro.caltech.edu,nblago@caltech.edu

"""


def MF_imcombine(objname, files, dependencies=""):

    filelist = " ".join(["%s " % ifile for ifile in files])
    first = "%s.fits: %s %s\n" % (objname, filelist, dependencies)

    if len(files) > 7:
        reject = "--Nlo 1 --Nhi 1"
    else:
        reject = ""
    if "bias" in objname:
        second = "\t$(IMCOMBINE) --outname %s.fits --listfile %s.lst %s --files %s\n" % (
            objname, objname, reject, filelist)
    else:
        second = "\t$(IMCOMBINE) --outname %s.fits --listfile %s.lst %s --files %s\n" % (
            objname, objname, reject, filelist)

    if "bias" not in objname and "dome" not in objname:
        second += "\n%s.npy : cube.npy %s.fits\n\t$(EXTSINGLE) cube.npy --A %s.fits --outname %s.npy --flat_correction flat-dome-700to900.npy --nosky\n" % (
            objname, objname, objname, objname)

    return first + second + "\n"


def MF_single(objname, obsnum, ifile, standard=None):
    """Create the MF entry for a observation with a single file. """

    # print(objname, obsnum, ifile)

    tp = {'objname': objname, 'obsfile': "bs_crr_b_%s" % ifile}
    tp['num'] = '_obs%s' % obsnum
    tp['outname'] = "%(objname)s%(num)s.npy" % tp
    tp['name'] = "%(objname)s%(num)s" % tp
    tp['specnam'] = "sp_%(objname)s%(num)s.npy" % tp

    if standard is None:
        tp['STD'] = ''
    else:
        tp['STD'] = "--std %s" % standard

    if 'PTF' in objname:
        tp['interact'] = '--interact'
    else:
        tp['interact'] = ''

    tp['flexname'] = "flex_bs_crr_b_%s.npy" % os.path.splitext(ifile)[0]

    first = """# %(outname)s
%(outname)s: cube.npy %(flexname)s %(obsfile)s.gz
\t$(EXTSINGLE) cube.npy --A %(obsfile)s.gz --outname %(outname)s %(STD)s --flat_correction flat-dome-700to900.npy --Aoffset %(flexname)s

sp_%(outname)s: %(outname)s
\t$(EXTSINGLE) cube.npy --A %(obsfile)s.gz --outname %(outname)s %(STD)s --flat_correction flat-dome-700to900.npy --Aoffset %(flexname)s --specExtract --autoExtract
\t$(PLOT) --spec %(specnam)s --savespec --savefig %(interact)s

redo_%(name)s:
\t@echo re-make-ing sp_%(outname)s
\t$(EXTSINGLE) cube.npy --A %(obsfile)s.gz --outname %(outname)s %(STD)s --flat_correction flat-dome-700to900.npy --Aoffset %(flexname)s --specExtract
\t$(PLOT) --spec %(specnam)s --savespec --savefig --interact

cube_%(outname)s.fits: %(outname)s
\t$(PY) $(PYC)/Cube.py %(outname)s --step extract --outname cube_%(outname)s.fits
""" % tp
    second = """corr_%(outname)s: %(outname)s
\t$(ATM) CORR --A %(outname)s --std %(objname)s --outname corr_%(outname)s\n""" % tp
    fn = "%(outname)s" % tp

    if standard is None:
        return first + "\n", fn
    else:
        return first + second + "\n", fn


def MF_standard(objname, obsnum, ifile, standard=None):
    """Create the MF entry for a standard star observation. """

    # print(objname, obsnum, ifile)

    tp = {'objname': objname, 'obsfile': "bs_crr_b_%s" % ifile}
    tp['num'] = '_obs%s' % obsnum
    tp['outname'] = "%(objname)s%(num)s.npy" % tp
    tp['name'] = "%(objname)s%(num)s" % tp
    tp['specnam'] = "sp_%(objname)s%(num)s.npy" % tp

    if standard is None:
        tp['STD'] = ''
        tp['specplot'] = ''
    else:
        tp['STD'] = "--std %s" % standard
        tp['specplot'] = "\t$(PLOT) --spec %(specnam)s --savespec --savefig" % tp

    tp['flexname'] = "flex_bs_crr_b_%s.npy" % os.path.splitext(ifile)[0]

    first = """# %(outname)s
%(outname)s: cube.npy %(flexname)s %(obsfile)s.gz
\t$(EXTSINGLE) cube.npy --A %(obsfile)s.gz --outname %(outname)s %(STD)s --flat_correction flat-dome-700to900.npy --Aoffset %(flexname)s
%(specplot)s

sp_%(outname)s: %(outname)s
\t$(EXTSINGLE) cube.npy --A %(obsfile)s.gz --outname %(outname)s %(STD)s --flat_correction flat-dome-700to900.npy --Aoffset %(flexname)s --specExtract
%(specplot)s

cube_%(outname)s.fits: %(outname)s
\t$(PY) $(PYC)/Cube.py %(outname)s --step extract --outname cube_%(outname)s.fits
""" % tp
    second = """corr_%(outname)s: %(outname)s
\t$(ATM) CORR --A %(outname)s --std %(objname)s --outname corr_%(outname)s\n""" % tp
    fn = "%(outname)s" % tp

    if standard is None:
        return first + "\n", fn
    else:
        return first + second + "\n", fn


def MF_AB(objname, obsnum, A, B):
    """Create the MF entry for an A-B observation"""

    # print(objname, obsnum, A, B)
    tp = {'objname': objname, 'A': "bs_crr_b_" + A, 'B': "bs_crr_b_" + B}
    if obsnum == 1:
        tp['num'] = ''
    else:
        tp['num'] = '_obs%i' % obsnum
    tp['outname'] = "%(objname)s%(num)s.npy" % tp
    tp['name'] = "%(objname)s%(num)s" % tp
    tp['specnam'] = "sp_%(objname)s%(num)s.npy" % tp
    # we only use the flexure from the A image
    tp['flexname'] = "flex_bs_crr_b_%s.npy" % os.path.splitext(A)[0]

    tp['bgdnameA'] = "bgd_%s.npy" % os.path.splitext(A)[0]
    tp['bgdnameB'] = "bgd_%s.npy" % os.path.splitext(B)[0]

    return """# %(outname)s\n%(outname)s: cube.npy %(A)s.gz %(B)s.gz %(flexname)s
\t$(EXTPAIR) cube.npy --A %(A)s.gz --B %(B)s.gz --outname %(outname)s --flat_correction flat-dome-700to900.npy --Aoffset %(flexname)s

sp_%(outname)s: %(outname)s
\t$(EXTPAIR) cube.npy --A %(A)s.gz --B %(B)s.gz --outname %(outname)s --flat_correction flat-dome-700to900.npy --Aoffset %(flexname)s --specExtract
\t$(PLOT) --spec %(specnam)s --savespec --savefig --interact

redo_%(name)s:
\t@echo re-make-ing sp_%(outname)s
\t$(EXTPAIR) cube.npy --A %(A)s.gz --B %(B)s.gz --outname %(outname)s --flat_correction flat-dome-700to900.npy --Aoffset %(flexname)s --specExtract
\t$(PLOT) --spec %(specnam)s --savespec --savefig --interact
\n\n""" % tp, "%(outname)s" % tp


def to_makefile(objs, calibs):

    MF = ""

    all = ""
    stds = ""
    stds_dep = ""
    sci = ""
    oth = ""
    auto = ""

    flexures = ""

    for calibname, files in calibs.items():

        if "bias" not in calibname:
            pass
        MF += MF_imcombine(calibname, files)
        all += "%s.fits " % calibname

    flatfiles = []
    for objname, observations in objs.items():

        objname = objname.replace(" ", "_")
        objname = objname.replace(")", "_")
        objname = objname.replace("(", "_")
        objname = objname.replace("[", "_")
        objname = objname.replace("]", "_")

        for obsnum, obsfiles in observations.items():
            flatfiles.append(obsfiles)

            # Handle Standard Stars
            if objname.startswith("STD-"):
                pred = objname[4:].rstrip().lower().replace("+", "").replace("-", "_")
                if pred in Stds.Standards:
                    standard = pred

                    obs = obsnum

                    for obsfile in obsfiles:
                        m, a = MF_standard(objname, "%i" % obs, obsfile,
                                           standard=standard)
                        # in case we have a STD A/B pair, we break them
                        # up into two separate observations
                        obs += 1
                        MF += m
                        # don't need these in all: dependants of target "stds"
                        # all += a + " "
                        stds_dep += a + " "

                else:
                    standard = None

                    for obsfile in obsfiles:
                        m, a = MF_single(objname, "%i" % obsnum, obsfile,
                                         standard=standard)
                        MF += m
                        oth += "sp_" + a + " "
                        auto += "sp_" + a + " "
                continue

            # Handle science targets
            if len(obsfiles) == 2:
                m, a = MF_AB(objname, obsnum, obsfiles[0], obsfiles[1])

                MF += m
                all += a + " "

                if not objname.startswith("STD-"):
                    if objname.startswith("PTF"):
                        sci += "sp_" + a + " "
                    else:
                        oth += "sp_" + a + " "
            else:
                for obsfile in obsfiles:
                    standard = None

                    m, a = MF_single(objname, "%i" % obsnum, obsfile)

                    if standard is not None:
                        stds += "corr_%s " % a

                    MF += m
                    all += a + " "

                    if not objname.startswith("STD-"):
                        if objname.startswith("PTF"):
                            sci += "sp_" + a + " "
                        else:
                            oth += "sp_" + a + " "
                            auto += "sp_" + a + " "
    stds += " "

    preamble = make_preamble

    f = open("Makefile", "w")
    clean = "\n\nclean:\n\trm %s %s" % (all, stds)
    science = "\n\nscience: %s report\n" % sci
    other = "\n\nother: %s report\n" % oth
    automatic = "\n\nauto: %s report\n" % auto
    corr = "std-correction.npy: %s \n\t$(ATM) CREATE --outname std-correction.npy --files sp_STD*npy \n" % stds_dep

    f.write(preamble + corr + "\nall: stds %s%s%s%s%s" % (all, clean, science,
                                                          other, automatic) +
            "\n" + MF + "\n" + flexures)
    f.close()


def make_plan(headers):
    """Convert headers to a makefile, assuming headers sorted by JD."""

    objs, calibs = identify_observations(headers)
    to_makefile(objs, calibs)


if __name__ == '__main__':

    files = sys.argv[1:]
    to_process = extract_info(files)

    make_plan(to_process)
