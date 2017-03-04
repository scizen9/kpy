import numpy as np
import matplotlib.pyplot as pl
from matplotlib import gridspec
from matplotlib.widgets import Cursor

import scipy
import scipy.spatial

from numpy.polynomial.chebyshev import chebval


def get_ellipse_xys(ell):
    a = ell[0]
    b = ell[1]

    pts = np.zeros((361, 2))
    beta = -ell[4] * np.pi / 180.
    sin_beta = np.sin(beta)
    cos_beta = np.cos(beta)
    alpha = np.radians(np.r_[0.:360.:1j * 361])

    sin_alpha = np.sin(alpha)
    cos_alpha = np.cos(alpha)

    pts[:, 0] = ell[2] + (a * cos_alpha * cos_beta - b * sin_alpha * sin_beta)
    pts[:, 1] = ell[3] + (a * cos_alpha * sin_beta + b * sin_alpha * cos_beta)

    return pts


class MouseCross(object):
    """ Draw a cursor with the mouse cursor """

    def __init__(self, ax, ellipse=None, nosky=False, **kwargs):
        self.ax = ax
        self.radius_as = ellipse[0]
        self.nosky = nosky
        self.ellipse = ellipse

        print("semimajor axis is %s arcsec" % self.radius_as)
        print("x - expand ap, z - shrink ap, y - toggle sky/host sub")

        marker = get_ellipse_xys(self.ellipse)
        self.line, = self.ax.plot(marker[:, 0], marker[:, 1], '-',
                                  visible=True, color='red', linewidth=2.,
                                  **kwargs)

    def show_cross(self, event):
        if event.inaxes == self.ax:
            self.line.set_visible(False)
            self.ellipse = (self.ellipse[0], self.ellipse[1],
                            event.xdata, event.ydata, self.ellipse[4])
            marker = get_ellipse_xys(self.ellipse)
            self.line, = self.ax.plot(marker[:, 0], marker[:, 1], '-',
                                      visible=True, color='red', linewidth=2.)
        else:
            self.line.set_visible(False)

        pl.draw()

    def size_cross(self, event):
        self.line.set_visible(False)
        if event.key == "z":
            self.radius_as -= 0.2
        elif event.key == "x":
            self.radius_as += 0.2
        elif event.key == "y":
            if self.nosky:
                self.nosky = False
            else:
                self.nosky = True

        print("semimajor axis is %s arcsec" % self.radius_as)

        self.ellipse = (self.radius_as,
                        self.radius_as * (self.ellipse[1]/self.ellipse[0]),
                        event.xdata, event.ydata, self.ellipse[4])
        marker = get_ellipse_xys(self.ellipse)
        self.line, = self.ax.plot(marker[:, 0], marker[:, 1], '-',
                                  visible=True, color='red', linewidth=2.)
        self.line.set_visible(True)

        if self.nosky:
            print("Sky subtraction off")
        else:
            print("Sky subtraction on")

        pl.draw()


class PositionPicker(object):
    """ This class is used to select an extraction point in a data cube """

    spectra = None
    Xs = None
    Ys = None
    Vs = None
    pointsize = None
    picked = None
    radius_as = None
    bgd_sub = False
    xc = None
    yc = None
    nosky = None
    scaled = None
    ellipse = None

    def __init__(self, spectra=None, pointsize=35, bgd_sub=False, ellipse=None,
                 objname=None, scaled=False,
                 lmin=600, lmax=650, cmin=-300, cmax=300, nosky=False):
        """ Create spectum picking gui.

        Args:
            spectra: SEDMr.Spectra object
        """

        self.spectra = spectra
        self.pointsize = pointsize
        self.scaled = scaled
        self.lmin = lmin
        self.lmax = lmax
        self.cmin = cmin
        self.cmax = cmax
        self.objname = objname
        self.bgd_sub = bgd_sub
        self.nosky = nosky
        self.radius_as = ellipse[0]
        self.ellipse = ellipse

        self.Xs, self.Ys, self.Vs = spectra.to_xyv(lmin=lmin, lmax=lmax)

        if bgd_sub:
            self.Vs -= np.nanmedian(self.Vs)

        pl.ioff()
        pl.title("%s Image from %s to %s nm" % (self.objname,
                                                self.lmin,
                                                self.lmax))
        self.figure = pl.figure(1)

        self.figure.canvas.mpl_connect("button_press_event", self)

        self.draw_cube()

    def draw_cube(self):

        if self.scaled:
            dv_min = self.cmin
            dv_max = self.cmax
        else:
            # get middle value
            if self.bgd_sub:
                v_mid = 0.
            else:
                v_mid = np.nanmedian(self.Vs)

            # get standard deviation
            v_std = np.nanstd(self.Vs)
            if 0 < v_std < 100:
                dv_min = v_mid - 3.0 * v_std
                dv_max = v_mid + 3.0 * v_std
            else:
                dv_min = -300
                dv_max = 300

        # plot (may want to use cmap=pl.cm.Spectral)
        print("scaling image display between %d and %d" % (dv_min, dv_max))
        pl.scatter(self.Xs, self.Ys, c=self.Vs, s=self.pointsize, linewidth=0,
                   vmin=dv_min, vmax=dv_max, cmap=pl.get_cmap('jet'))

        pl.ylim(-14, 14)
        pl.xlim(14, -14)
        pl.xlabel("-RA offset [asec]")
        pl.ylabel("Dec offset [asec]")
        pl.colorbar()

        # c = Cursor(self.figure.gca(), useblit=True)

        cross = MouseCross(self.figure.gca(), ellipse=self.ellipse,
                           nosky=self.nosky)
        self.figure.canvas.mpl_connect('motion_notify_event', cross.show_cross)
        self.figure.canvas.mpl_connect("key_press_event", cross.size_cross)
        pl.show()
        self.radius_as = cross.radius_as
        self.nosky = cross.nosky
        self.ellipse = cross.ellipse

    def __call__(self, event):
        """Event call handler for Picker gui."""

        if event.name == 'button_press_event':
            print("X = %+10.5f, Y = %+10.5f" % (event.xdata, event.ydata))
            self.picked = (event.xdata, event.ydata)
            self.xc = event.xdata
            self.yc = event.ydata
            pl.close(self.figure)


class ScaleCube(object):
    """ This class is used to scale a data cube """

    spectra = None
    Xs = None
    Ys = None
    Vs = None
    pointsize = None
    bgd_sub = False
    lmin = None
    lmax = None
    cmin = None
    cmax = None

    def __init__(self, spectra=None, pointsize=35, bgd_sub=False,
                 objname=None, lmin=600, lmax=650):
        """ Create scaling gui.

        Args:
            spectra: SEDMr.Spectra object
        """

        print("First scale cube using keys to change limits:")
        if bgd_sub:
            print("> - increase upper/lower spread by 200")
            print("< - decrease upper/lower spread by 200")
        else:
            print("> - to increase upper limit by 100")
            print("< - to decrease upper limit by 100")
            print(". - to increase lower limit by 100")
            print(", - to decrease lower limit by 100")
        print("x - exit")
        print("q - to abandon scaling")

        self.spectra = spectra
        self.pointsize = pointsize
        self.lmin = lmin
        self.lmax = lmax
        self.objname = objname
        self.bgd_sub = bgd_sub
        self.scaled = False
        self.scat = None
        self.cb = None

        self.Xs, self.Ys, self.Vs = spectra.to_xyv(lmin=lmin, lmax=lmax)

        if bgd_sub:
            self.Vs -= np.nanmedian(self.Vs)

        # get standard deviation
        v_std = np.nanstd(self.Vs)

        # get middle value
        if self.bgd_sub:
            v_mid = 0.
            if 0 < v_std < 100:
                self.cmin = v_mid - 3.0 * v_std
                self.cmax = v_mid + 3.0 * v_std
            else:
                self.cmin = -300
                self.cmax = 300
        else:
            v_mid = np.nanmedian(self.Vs)
            self.cmin = v_mid - 3.0 * v_std
            self.cmax = v_mid + 3.0 * v_std

        print("mid, std: %f, %f" % (v_mid, v_std))

        pl.ioff()

        self.figure = pl.figure(1)

        self.figure.canvas.mpl_connect("key_press_event", self)

        self.draw_cube()

    def draw_cube(self):

        pl.title("Scaling %s Image from %s to %s nm\nfrom %.1f to %.1f int" %
                 (self.objname, self.lmin, self.lmax, self.cmin, self.cmax))

        self.scat = pl.scatter(self.Xs, self.Ys, c=self.Vs, s=self.pointsize,
                               linewidth=0, vmin=self.cmin, vmax=self.cmax,
                               cmap=pl.get_cmap('jet'))

        pl.ylim(-14, 14)
        pl.xlim(14, -14)
        pl.xlabel("RA offset [asec]")
        pl.ylabel("Dec offset [asec]")
        self.cb = self.figure.colorbar(self.scat)
        pl.show()

    def update_cube(self):

        ax = self.figure.gca()
        ax.set_title("Scaling %s Image from %s to %s nm\nfrom %.1f to %.1f Irr"
                     % (self.objname, self.lmin, self.lmax,
                        self.cmin, self.cmax))

        self.scat.remove()

        self.scat = ax.scatter(self.Xs, self.Ys, c=self.Vs,
                               s=self.pointsize, linewidth=0,
                               vmin=self.cmin, vmax=self.cmax,
                               cmap=pl.get_cmap('jet'))
        self.cb.set_clim(self.cmin, self.cmax)
        self.cb.update_normal(self.scat)
        self.figure.canvas.draw()

    def __call__(self, event):
        """Event call handler for scaling gui."""

        if event.key == 'x':
            self.scaled = True
            print("Scaling between %f and %f" % (self.cmin, self.cmax))
            pl.close(self.figure)
        if event.key == 'q':
            self.scaled = False
            print("Using default scaling")
            pl.close(self.figure)
        elif event.key == '>':
            if self.bgd_sub:
                if self.cmax > 100. and self.cmin < -100.:
                    self.cmax += 100.
                    self.cmin -= 100.
                else:
                    self.cmax += 10.
                    self.cmin -= 10.
            else:
                self.cmax += 100.
            self.update_cube()
        elif event.key == '<':
            if self.bgd_sub:
                if self.cmax > 100. and self.cmin < -100.:
                    self.cmax -= 100.
                    self.cmin += 100.
                elif self.cmax > 10. and self.cmin < -10.:
                    self.cmax -= 10.
                    self.cmin += 10.
            else:
                if self.cmax > 100.:
                    self.cmax -= 100.
            self.update_cube()
        elif event.key == '.':
            if not self.bgd_sub:
                self.cmin += 100.
                self.update_cube()
        elif event.key == ',':
            if not self.bgd_sub:
                self.cmin -= 100.
                self.update_cube()


class WaveFixer(object):
    """ This class is used to fix bad wavelength solutions """

    cube = None  # Raw data cube spectra
    KT = None  # KDTree object
    X1 = []
    X2 = []
    Y1 = []

    pointsize = None
    picked = None

    state = "Display"

    fig = None
    ax_cube = None
    ax_spec = None

    def __init__(self, cube=None, pointsize=35):
        """ Create spectum picking gui.

        Args:
            cube: Data cube list
                
        """

        self.actions = {"m": self.mode_switch}

        self.cube = cube
        self.pointsize = pointsize

        for ix, s in enumerate(self.cube):
            if s.xrange is None:
                continue
            xs = np.arange(s.xrange[0], s.xrange[1], .1)

            if s.lamcoeff is not None:
                lls = chebval(xs, s.lamcoeff)
                ha1 = np.argmin(np.abs(lls - 656.3)) / 10.0 + s.xrange[0]

            if s.mdn_coeff is not None:
                lls = chebval(xs, s.mdn_coeff)
                ha2 = np.argmin(np.abs(lls - 656.3)) / 10.0 + s.xrange[0]

            self.X1.append(ha1)
            self.X2.append(ha2)
            self.Y1.append(s.yrange[0])

        self.X1, self.X2, self.Y1 = map(np.array, [self.X1, self.X2,
                                                   self.Y1])

        ok = (np.abs(self.X1 - self.X2) < 2) & np.isfinite(self.X2) & \
             np.isfinite(self.Y1)

        self.good_cube = self.cube[ok]
        locs = np.array([self.X2[ok], self.Y1[ok]]).T
        self.KT = scipy.spatial.KDTree(locs)

        assert (len(locs) == len(self.good_cube))

        # Setup drawing
        gs = gridspec.GridSpec(1, 2, width_ratios=[1, 2.5])
        self.fig = pl.figure(1, figsize=(22, 5.5))
        self.ax_cube = pl.subplot(gs[0])
        self.ax_spec = pl.subplot(gs[1])

        self.ax_cube.set_xlim(-100, 2200)
        self.ax_cube.set_ylim(-100, 2200)

        pl.ion()
        self.draw_cube()

        print("Registering events")
        self.fig.canvas.mpl_connect("button_press_event", self)
        self.fig.canvas.mpl_connect("key_press_event", self)

        pl.show()

    def mode_switch(self):
        """ Toggle operating mode between Display and Select """

        if self.state == "Display":
            self.state = "Select"
        else:
            self.state = "Display"

        print(self.state)
        self.draw_cube()

    def draw_spectra(self):
        """ Draw nearest spectra """
        if self.picked is None:
            return

        print("Drawing spectra")
        # xl = self.ax_spec.get_xlim()
        # yl = self.ax_spec.get_ylim()
        self.ax_spec.cla()
        # self.ax_spec.set_xlim(xl)
        # self.ax_spec.set_ylim(yl)

        x, y = self.X2[self.picked], self.Y1[self.picked]
        objs = self.KT.query_ball_point((x, y), 70)

        print("Query around %s found: %s" % ((x, y), objs))

        spec = self.cube[self.picked]
        ix = np.arange(*spec.xrange)
        fiducial_ll = chebval(ix, spec.lamcoeff)
        # self.ax_spec.plot(ix-spec.xrange[0], fiducial_ll, linewidth=3)
        self.ax_spec.step(fiducial_ll, spec.spec, linewidth=3)

        for spec in self.good_cube[objs]:
            try:
                ix = np.arange(*spec.xrange)
                ll = chebval(ix, spec.lamcoeff)
                # self.ax_spec.plot(ix-spec.xrange[0], ll-fiducial_ll)
                self.ax_spec.step(ll, spec.spec)
            except:
                pass

        # self.ax_spec.set_ylim(-30,30)
        self.ax_spec.set_xlim(370, 700)
        self.fig.show()

    def draw_cube(self):
        """ Draw the data cube """
        print("drawing cube")

        # Draw cube

        xl = self.ax_cube.get_xlim()
        yl = self.ax_cube.get_ylim()
        self.ax_cube.cla()
        self.ax_cube.set_xlim(xl)
        self.ax_cube.set_ylim(yl)

        self.ax_cube.plot(self.X2, self.Y1, 'o',
                          markersize=8, marker='h', linewidth=0)

        bad = np.abs(self.X1 - self.X2) > 4
        self.ax_cube.plot(self.X2[bad], self.Y1[bad], 'ro',
                          markersize=7, marker='h', linewidth=0)

        if self.picked is not None:
            print(self.X2[self.picked], self.Y1[self.picked])

            self.ax_cube.plot([self.X2[self.picked]],
                              [self.Y1[self.picked]],
                              'o', ms=12, color='yellow',
                              alpha=0.4, visible=True)

        tit = "State: %s | Press ? for help" % self.state
        self.ax_cube.set_title(tit)

        self.fig.show()

    def handle_button_press(self, event):

        if event.inaxes == self.ax_cube:
            """Clicked In Data Cube Display"""

            if self.state == 'Display':
                """ Display state (not pick state, ignore) """
                return

            dists = np.abs(self.X2 - event.xdata) + \
                    np.abs(self.Y1 - event.ydata)

            ix = np.nanargmin(dists)
            print(dists[ix])

            if dists[ix] > 20:
                self.picked = None
            else:
                self.picked = ix

            self.draw_cube()
            self.draw_spectra()

    def __call__(self, event):
        """Event call handler for Picker gui."""

        print(event.name)

        if event.name == 'pick_event':
            import pdb
            pdb.set_trace()

        elif event.name == 'button_press_event':
            """ Note order of if statement to skip button over pick event"""

            self.handle_button_press(event)

        elif event.name == 'key_press_event':
            key = event.key

            if key in self.actions:
                to_call = self.actions[key]
                to_call()

            if key == "?":
                for k, v in self.actions.items():
                    print("%s: %s" % (k, v.__doc__))
