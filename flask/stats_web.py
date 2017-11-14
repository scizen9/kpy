''' Modified from online verion of:

 _README: https://github.com/bokeh/bokeh/blob/master/examples/app/stocks/README.md
 
 
.. note::
    Running this example requires having the "stats.log" file.

Use the ``bokeh serve`` command to run the example by executing:

    bokeh serve stocks

at your command prompt. Then navigate to the URL

    http://localhost:5006/stocks
..

'''
try:
    from functools import lru_cache
except ImportError:
    # Python 2 does stdlib does not have lru_cache so let's just
    # create a dummy decorator to avoid crashing
    print ("WARNING: Cache for this example is available on Python 3 only.")
    def lru_cache():
        def dec(f):
            def _(*args, **kws):
                return f(*args, **kws)
            return _
        return dec

from os.path import dirname, join

import pandas as pd
import datetime
import numpy as np
import os

from bokeh.io import curdoc
from bokeh.layouts import row, column
from bokeh.models import ColumnDataSource, CDSView, GroupFilter, Label
from bokeh.models.widgets import PreText, Select
from bokeh.plotting import figure
from astropy.time import Time




@lru_cache()
def load_stats(statsfile='stats.log'):

    data = pd.read_csv(statsfile, header=None, 
                       names=['path', 'obj', 'jd', 'ns', 'fwhm', 'ellipticity', 'bkg', 'airmass', 'in_temp', 'imtype', 'out_temp', 'in_hum'])
    jds = data['jd']
    t = Time(jds, format='jd', scale='utc')
    date = t.utc.datetime
    day_frac_diff = datetime.timedelta(np.ceil((datetime.datetime.now() - datetime.datetime.utcnow() ).total_seconds())/3600/24)
    local_date = date + day_frac_diff
     
    data2 = data.assign(localdate=local_date)
    data2.set_index('localdate')
    return pd.DataFrame({'date':data2['localdate'], 'ns':data2['ns'], 'fwhm':data2['fwhm'], 'ellipticity':data2['ellipticity'], \
    'bkg':data2['bkg'], 'airmass':data2['airmass'], 'in_temp':data2['in_temp'], 'imtype':data2['imtype'],\
    'out_temp':data2['out_temp'], 'in_hum':data2['in_hum']})
    
    

@lru_cache()
def plot_stats(statsfile):


    source = ColumnDataSource(data=dict(date=[], ns=[], fwhm=[], ellipticity=[], bkg=[], airmass=[], in_temp=[], imtype=[], out_temp=[], in_hum=[]))
    source_static = ColumnDataSource(data=dict(date=[], ns=[], fwhm=[], ellipticity=[], bkg=[], airmass=[], in_temp=[], imtype=[], out_temp=[], in_hum=[]))
    viewScience = CDSView(source=source, filters=[GroupFilter(column_name='imtype', group='SCIENCE')])
    viewAcquisition = CDSView(source=source, filters=[GroupFilter(column_name='imtype', group='ACQUISITION')])
    viewGuider = CDSView(source=source, filters=[GroupFilter(column_name='imtype', group='GUIDER')])
    viewFocus = CDSView(source=source, filters=[GroupFilter(column_name='imtype', group='FOCUS')])
    
    tools = 'pan,wheel_zoom,xbox_select,reset'
    
    def update(selected=None):
    
        data = load_stats(statsfile)
        source.data = source.from_df(data[['date', 'ns', 'fwhm', 'ellipticity', 'bkg', 'airmass', 'in_temp', 'imtype', 'out_temp', 'in_hum']])
        source_static.data = source.data
    
        #update_stats(data, t1, t2)


    ns = figure(plot_width=450, plot_height=250, tools=tools, x_axis_type='datetime', active_drag="xbox_select")
    ns.line('date', 'ns', source=source_static)
    ns.circle('date', 'ns', size=1, source=source, color=None, selection_color="orange")
    ns.title.text =  "Number of bright sources extracted"
    
    bkg = figure(plot_width=450, plot_height=250, tools=tools, x_axis_type='datetime', active_drag="xbox_select")
    bkg.x_range = ns.x_range
    bkg.line('date', 'bkg', source=source_static)
    bkg.circle('date', 'bkg', size=1, source=source, color=None, selection_color="orange")
    bkg.title.text =  "Background (counts)"
    
    
    temp = figure(plot_width=450, plot_height=250, tools=tools, x_axis_type='datetime', active_drag="xbox_select")
    temp.x_range = ns.x_range
    temp.line('date', 'in_temp', source=source_static, color='blue', legend="Inside")
    temp.line('date', 'out_temp', source=source_static, color='green', legend="Outside")
    temp.circle('date', 'in_temp', size=1, source=source, color=None, selection_color="orange")
    temp.title.text =  "Temperature"
    
    
    fwhm = figure(plot_width=450, plot_height=250, tools=tools, x_axis_type='datetime', active_drag="xbox_select")
    fwhm.x_range = ns.x_range
    fwhm.circle('date', 'fwhm', source=source_static, color="green", legend="Focus", view=viewFocus)
    fwhm.circle('date', 'fwhm', source=source_static, color="red", legend="Science", view=viewScience)
    fwhm.circle('date', 'fwhm', source=source_static, color="blue", legend="Acquisition", view=viewAcquisition)
    fwhm.circle('date', 'fwhm', source=source_static, color="black", legend="Guider", view=viewGuider)
    
    fwhm.circle('date', 'fwhm', size=1, source=source, color=None, selection_color="orange")
    fwhm.title.text = "FWHM [arcsec]"
    
    airmass = figure(plot_width=450, plot_height=250, tools=tools, x_axis_type='datetime', active_drag="xbox_select")
    airmass.x_range = ns.x_range
    airmass.line('date', 'airmass', source=source_static)
    airmass.circle('date', 'airmass', size=1, source=source, color=None, selection_color="orange")
    airmass.title.text = "Airmass"
    
    ellipticity = figure(plot_width=450, plot_height=250, tools=tools, x_axis_type='datetime', active_drag="xbox_select")
    ellipticity.x_range = ns.x_range
    ellipticity.line('date', 'ellipticity', source=source_static)
    ellipticity.circle('date', 'ellipticity', size=1, source=source, color=None, selection_color="orange")
    ellipticity.title.text = "Ellipticity"
    
    
    humidity = figure(plot_width=450, plot_height=250, tools=tools, x_axis_type='datetime', active_drag="xbox_select")
    humidity.x_range = ns.x_range
    humidity.line('date', 'in_hum', source=source_static)
    humidity.circle('date', 'in_hum', size=1, source=source, color=None, selection_color="orange")
    humidity.title.text = "Inside Humidity"
    # set up callbacks
    
    
    left = column(ns, bkg,temp, humidity)
    right = column(fwhm, airmass, ellipticity)
    layout = row(left, right)
    
    # initialize
    update()
    
    curdoc().add_root(layout)
    curdoc().title = "Stats"
    return layout

    #output_file("templates/stats.html", title="Night Statistics")
    #show(layout)

    return curdoc()

@lru_cache()
def plot_not_found_message(day):
    not_found = figure(plot_width=900, plot_height=450, x_range=[0, 900], y_range=[0, 450])
    not_found.image(image=[np.zeros([900, 450])+0.1], x=0, y=0, dw=900, dh=450)
    citation = Label(x=50, y=225, x_units='screen', y_units='screen', text='No statistics found for today \n (likely we were weathered out...)')
    not_found.add_layout(citation)
    not_found.title.text = "Statistics not found for day %s"%(day)

    layout = column(not_found)
    curdoc().add_root(layout)
    curdoc().title = "Stats not found"
    
    
#if __name__ == '__main__':
day = ("%s"%(datetime.datetime.utcnow())).split()[0]
    
#This value is for testing, use the value below in production
#datadir = '/home/nblago/workspace/kpy/tests'
datadir = join(dirname('/scr2/sedm/phot/'), day, 'stats/stats.log')

daylog = join(datadir, 'stats.log')

if (os.path.isfile(daylog)):
    plot_stats(daylog)
else:
    plot_not_found_message(day)
