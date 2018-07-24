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
from bokeh.models import ColumnDataSource, Label, CDSView, GroupFilter, Range1d, LinearAxis
from bokeh.models import HoverTool
from bokeh.models.annotations import BoxAnnotation
from bokeh.models.widgets import PreText, Select
from bokeh.plotting import figure
from bokeh.core.properties import value
from bokeh.palettes import Paired


from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun, get_moon
import astropy.units as u

#from datetime import datetime, timedelta
import datetime
import model


@lru_cache()
def load_p48seeing(obsdate):

    time, seeing = model.get_p18obsdata(obsdate)
    day_frac_diff = datetime.timedelta(np.ceil((datetime.datetime.now() - datetime.datetime.utcnow() ).total_seconds())/3600/24)
    local_date = np.array(time) + day_frac_diff
    d = pd.DataFrame({'date':local_date, 'seeing':seeing})

    return d


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
def plot_stats(statsfile, mydate):


    source = ColumnDataSource(data=dict(date=[], ns=[], fwhm=[], ellipticity=[], bkg=[], airmass=[], in_temp=[], imtype=[], out_temp=[], in_hum=[]))
    source_static = ColumnDataSource(data=dict(date=[], ns=[], fwhm=[], ellipticity=[], bkg=[], airmass=[], in_temp=[], imtype=[], out_temp=[], in_hum=[]))

    viewScience = CDSView(source=source, filters=[GroupFilter(column_name='imtype', group='SCIENCE')])
    viewAcquisition = CDSView(source=source, filters=[GroupFilter(column_name='imtype', group='ACQUISITION')])
    viewGuider = CDSView(source=source, filters=[GroupFilter(column_name='imtype', group='GUIDER')])
    viewFocus = CDSView(source=source, filters=[GroupFilter(column_name='imtype', group='FOCUS')])
    source_p48 = ColumnDataSource(data=dict(date=[], seeing=[]))

    def update(selected=None):
    
        if statsfile:
            data = load_stats(statsfile)
            source.data = source.from_df(data[['date', 'ns', 'fwhm', 'ellipticity', 'bkg', 'airmass', 'in_temp', 'imtype', 'out_temp', 'in_hum']])
            source_static.data = source.data
    
        p48 = load_p48seeing(mydate)
        source_p48.data = source_p48.from_df(p48[['date', 'seeing']])
        source_static_p48.data = source_p48.data

    source_static_p48 = ColumnDataSource(data=dict(date=[], seeing=[]))        
    tools = 'pan,box_zoom,reset'
    
    
    p48seeing = figure(plot_width=425, plot_height=250, tools=tools, x_axis_type='datetime', active_drag="box_zoom")
    p48seeing.circle('date', 'seeing', source=source_static_p48, color="black")
    p48seeing.title.text = "P18 seeing [arcsec]"

    if statsfile:
        ns = figure(plot_width=425, plot_height=250, tools=tools, x_axis_type='datetime', active_drag="box_zoom")
        ns.line('date', 'ns', source=source_static)
        ns.circle('date', 'ns', size=1, source=source, color=None, selection_color="orange")
        ns.title.text =  "Number of bright sources extracted"
        

        bkg = figure(plot_width=425, plot_height=250, tools=tools, x_axis_type='datetime', active_drag="box_zoom")
        bkg.x_range = ns.x_range
        bkg.line('date', 'bkg', source=source_static)
        bkg.circle('date', 'bkg', size=1, source=source, color=None, selection_color="orange")
        bkg.title.text =  "Background (counts)"
        
        
        temp = figure(plot_width=425, plot_height=250, tools=tools, x_axis_type='datetime', active_drag="box_zoom")
        temp.x_range = ns.x_range
        temp.line('date', 'in_temp', source=source_static, color='blue', legend="Inside")
        temp.line('date', 'out_temp', source=source_static, color='green', legend="Outside")
        temp.circle('date', 'in_temp', size=1, source=source, color=None, selection_color="orange")
        temp.title.text =  "Temperature [C]"
        
        
        fwhm = figure(plot_width=425, plot_height=250, tools=tools, x_axis_type='datetime', active_drag="box_zoom")
        fwhm.x_range = ns.x_range
        fwhm.circle('date', 'fwhm', source=source_static, color="green", legend="Focus", view=viewFocus)
        fwhm.circle('date', 'fwhm', source=source_static, color="red", legend="Science", view=viewScience)
        fwhm.circle('date', 'fwhm', source=source_static, color="blue", legend="Acquisition", view=viewAcquisition)
        fwhm.circle('date', 'fwhm', source=source_static, color="black", legend="Guider", view=viewGuider)
        fwhm.circle('date', 'fwhm', size=1, source=source, color=None, selection_color="orange")
        fwhm.title.text = "P60 FWHM [arcsec]"
        
        airmass = figure(plot_width=425, plot_height=250, tools=tools, x_axis_type='datetime', active_drag="box_zoom")
        airmass.x_range = ns.x_range
        airmass.line('date', 'airmass', source=source_static)
        airmass.circle('date', 'airmass', size=1, source=source, color=None, selection_color="orange")
        airmass.title.text = "Airmass"
        
        ellipticity = figure(plot_width=425, plot_height=250, tools=tools, x_axis_type='datetime', active_drag="box_zoom")
        ellipticity.x_range = ns.x_range
        ellipticity.line('date', 'ellipticity', source=source_static)
        ellipticity.circle('date', 'ellipticity', size=1, source=source, color=None, selection_color="orange")
        ellipticity.title.text = "Ellipticity"
        
        
        humidity = figure(plot_width=425, plot_height=250, tools=tools, x_axis_type='datetime', active_drag="box_zoom")
        humidity.x_range = ns.x_range
        humidity.line('date', 'in_hum', source=source_static)
        humidity.circle('date', 'in_hum', size=1, source=source, color=None, selection_color="orange")
        humidity.title.text = "Inside Humidity [%]"
    
        p48seeing.x_range = ns.x_range

        left = column(fwhm, p48seeing, airmass)
        center = column(ellipticity, ns, bkg, )
        right = column(temp, humidity)
        layout = row(left, center, right)

    else:

        layout = row(column(p48seeing))
    
    # initialize
    update()
    
    curdoc().add_root(layout)
    curdoc().title = "Stats"


    return layout

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
    

@lru_cache()
def plot_stats_allocation(data):
    """
    Plots in the shape of bars the time available and spent for each active allocation.
    """


    #Create the first plot with the allocation hours
    alloc_names = data['allocations']
    categories = ["spent_hours", "free_hours"]
    colors = [ "#e84d60", "darkgreen"] #"#c9d9d3"

    N = len(alloc_names)

    source = ColumnDataSource(data=data)
    p = figure(x_range=alloc_names, plot_height=420, plot_width=80*8, title="Time spent/available for SEDM allocations this term",
               toolbar_location=None, tools="")

    p.vbar_stack(categories, x='allocations', width=0.9, color=colors, source=source, legend=["Spent", "Available"])
    p.y_range.start = 0
    p.x_range.range_padding = 0.1
    p.xgrid.grid_line_color = None
    p.axis.minor_tick_line_color = None
    p.outline_line_color = None
    p.legend.location = "top_right"
    p.legend.orientation = "horizontal"
    p.yaxis.axis_label = 'Hours'
    p.xaxis.major_label_orientation = 0.3

    #Create the second plot with the % spent
    alloc_names = data['allocations']
    percentage = (data["spent_hours"] / data["alloc_hours"]) * 100

    colors=N*['#084594']
    '''for i, p in enumerate(percentage):
        if p<50: colors[i] = '#22A784'
        elif p>50 and p<75: colors[i] = '#FD9F6C'
        else: colors[i] = '#DD4968'''

    source = ColumnDataSource(data=dict(alloc_names=alloc_names, percentage=percentage, color=colors))

    p2 = figure(x_range=alloc_names, y_range=(0,100), plot_height=420, plot_width=80*8, title="Percentage of time spent",
               toolbar_location=None, tools="")

    p2.vbar(x='alloc_names', top='percentage', width=0.9, color='color', source=source)

    p2.xgrid.grid_line_color = None
    p2.legend.orientation = "horizontal"
    p2.legend.location = "top_center"
    p2.yaxis.axis_label = '% time spent'
    p2.xaxis.major_label_orientation = 0.3

    
    #Create the pie charts
    pieColors = 10*["red", "green", "blue", "orange", "yellow", 'lime', 'brown', 'cyan', \
        'magenta', 'olive', 'black', 'teal', 'gold', 'crimson', 'moccasin', 'greenyellow', 'navy', 'ivory', 'lightpink']

    #First one with the time spent

    # define starts/ends for wedges from percentages of a circle
    percents_only = np.round( np.array(list(data["spent_hours"] / np.sum(data["spent_hours"])))*100, 1)
    percents = np.cumsum( [0] + list(data["spent_hours"] / np.sum(data["spent_hours"])))
    starts = [per*2*np.pi for per in percents[:-1]]
    ends = [per*2*np.pi for per in percents[1:]]

    p3 = figure(x_range=(-1, 2.5), y_range=(-1.1, 1.1), plot_height=420, plot_width=600, title="% spent")

    #Add individual wedges:
    for i in range(N):
        p3.wedge(x=0, y=0, radius=.9, start_angle=starts[i], end_angle=ends[i], color=pieColors[i], legend="[{0}%] {1}".format(percents_only[i], alloc_names[i]) )

    p3.xgrid.grid_line_color = None
    p3.ygrid.grid_line_color = None
    p3.legend.orientation = "vertical"
    p3.legend.location = "top_right"
    p3.legend.border_line_alpha = 0
    p3.legend.background_fill_color = None
    p3.xaxis.visible = False
    p3.yaxis.visible = False

    #Second one with the time allocated

    # define starts/ends for wedges from percentages of a circle
    percents_only = np.round( np.array(list(data["alloc_hours"] / np.sum(data["alloc_hours"])))*100, 1)
    percents = np.cumsum( [0] + list(data["alloc_hours"] / np.sum(data["alloc_hours"])))
    starts = [per*2*np.pi for per in percents[:-1]]
    ends = [per*2*np.pi for per in percents[1:]]

    p4 = figure(x_range=(-1, 2.5), y_range=(-1.1, 1.1), plot_height=420, plot_width=600, title="% time allocated to each program")
    #Add individual wedges:
    for i in range(N):
        p4.wedge(x=0, y=0, radius=.9, start_angle=starts[i], end_angle=ends[i], color=pieColors[i], legend="[{0}%] {1}".format(percents_only[i], alloc_names[i]) ) 

    p4.xgrid.grid_line_color = None
    p4.ygrid.grid_line_color = None
    p4.legend.orientation = "vertical"
    p4.legend.location = "top_right"
    p4.legend.border_line_alpha = 0
    p4.legend.background_fill_color = None
    p4.xaxis.visible = False
    p4.yaxis.visible = False



    layout =  row(column(p, p2), column(p4, p3))
    

    curdoc().add_root(layout)
    curdoc().title = "Allocation stats"

    return layout


def plot_visibility(ras, decs, names, allocs=[None], priorities=[5], endobs=[None], exptime=2430, date=None):
    ''' makes a visibility plot for one or many objects, highlighting observed patches if relevant
    all arguments are arrays, even if they are of size 1
    priorities:   integers
    obsd:         list/array of observed objects, should match 'names'
    endobs:       'YYYY-MM-DDTHH:MM:SS.ssssssssss' (as outputed from sql query)
    exptime:      in seconds
    date:         YYYYMMDD, conveniently matching the folder name'''
    
    allocpalette = Paired[12][1::2] + Paired[12][::2]
    priorities = np.array(priorities, dtype=np.int)
    allocs = np.asarray(allocs)
    names = np.asarray(names)

    p = figure(plot_width=700, plot_height=500, toolbar_location='above',
               y_range=(0, 90), y_axis_location="right")

    ### setup with axes, sun/moon, frames, background 
    palomar_mountain = EarthLocation(lon=243.1361*u.deg, lat=33.3558*u.deg, height=1712*u.m)
    utcoffset = -7 * u.hour  # Pacific Daylight Time
    
    if date == None:
        time = (Time.now() - utcoffset).datetime # date is based on local time
        time = Time(datetime.datetime(time.year, time.month, time.day)) 
    else:
        time = Time(datetime.datetime(int(date[:4]), int(date[4:6]), int(date[6:8])))
    midnight = time - utcoffset # 7am local time of correct date, midnight UTC

    if endobs[0] != None:
        print endobs[0]
        endobs = Time(np.array(endobs, dtype='|S32'), format='isot')
        print endobs[0]
        endobs.format = u'datetime'
        print endobs[0], 'finished prep'

    delta_midnight = np.linspace(-8, 8, 200) * u.hour
    t = midnight + delta_midnight
    abstimes = [i.datetime.strftime('%I:%M %p') for i in t + utcoffset]
    frame = AltAz(obstime=t, location=palomar_mountain)
    sun_alt  =  get_sun(t).transform_to(frame).alt
    moon_alt = get_moon(t).transform_to(frame).alt
    
    # shading for nighttime and twilight
    dark_times    = delta_midnight[sun_alt < 0].value
    twilit_times  = delta_midnight[sun_alt < -18 * u.deg].value
    plotted_times = delta_midnight[sun_alt <   5 * u.deg].value
    
    twilight = BoxAnnotation(left=min(twilit_times), right=max(twilit_times), bottom=0, 
                             fill_alpha=0.15, fill_color='black', level='underlay')
    night    = BoxAnnotation(left=min(dark_times),    right=max(dark_times),    bottom=0, 
                             fill_alpha=0.25, fill_color='black', level='underlay')
    earth    = BoxAnnotation(top=0, fill_alpha=0.8, fill_color='sienna')
    
    p.add_layout(night)
    p.add_layout(twilight)
    p.add_layout(earth)
    
    # sun and moon
    sun  = p.line(delta_midnight, sun_alt,  line_color='red', name="Sun", legend='Sun', line_dash='dashed')
    moon = p.line(delta_midnight, moon_alt, line_color='yellow', line_dash='dashed', 
                                                   name="Moon", legend='Moon')
    # labels and axes
    p.title.text = "Visibility for %s UTC" %midnight
    p.xaxis.axis_label = "Hours from PDT Midnight"
    p.x_range.start = min(plotted_times)
    p.x_range.end   = max(plotted_times)
    p.yaxis.axis_label = "Airmass"
    
    # primary airmass label on right
    airmasses = (1.01, 1.1, 1.25, 1.5, 2., 3., 6.)
    ticker = [90 - np.arccos(1./i) * 180/np.pi for i in airmasses]
    p.yaxis.ticker = ticker
    p.yaxis.major_label_overrides = {tick: str(airmasses[i]) for i, tick in enumerate(ticker)}
    
    # add supplementary alt label on left
    p.extra_y_ranges = {"altitude": Range1d(0, 90)}
    p.add_layout(LinearAxis(y_range_name="altitude", axis_label='Altitude [deg]'), 'left')


    ### adding data from the actual objects
    objs = SkyCoord(np.array(ras,  dtype=np.float), 
                    np.array(decs, dtype=np.float), unit="deg")
    alloc_color = {}
    for i, val in enumerate(np.unique(allocs)):
        alloc_color[val] = allocpalette[i % len(allocpalette)]
        
    tooltipped = [] # things with tooltips
    tooltips = [('obj',        '@name'), # make it #name when we get to bokeh 0.13
                ('time',       '@abstime'), 
                ('altitude',   u"@alt\N{DEGREE SIGN}"), 
                ('airmass',    '@airmass'),
                ('priority',   '@priority'), 
                ('allocation', '@alloc')]
        
    for i in np.array(allocs).argsort(): # go in order by alloc for an alphabetized legend
        color = alloc_color[allocs[i]]
        obj = objs[i].transform_to(frame)
        source = ColumnDataSource(    dict(times=delta_midnight, 
                                             alt=obj.alt,
                                         airmass=obj.secz,
                                         abstime=abstimes,
                                        priority=np.full(len(abstimes), priorities[i]),
                                           alloc=np.full(len(abstimes), allocs[i]),
                                            name=np.full(len(abstimes), names[i]))) # delete the name when we get to bokeh 0.13
        if allocs[i] == None: # single object
            legend = names[i]
            tooltips = tooltips[:4]
        else:
            legend = '{}'.format(allocs[i])
            if endobs[0] != None: # plot that highlights observed part of the night
                # full path of the night
                dotted = p.line('times', 'alt', color=color, source=source, line_dash='2 2',
                                name=names[i], line_width=1, legend=legend)
                # manually crop the source so only thick observed part has tooltips
                endtime = endobs[i]
                initime = endtime - exptime * u.second
                if i > 0:
                    initime = max(initime, endobs[i - 1])
                mask = np.logical_and(delta_midnight + midnight + utcoffset > initime,
                                      delta_midnight + midnight + utcoffset < endtime)
                source = ColumnDataSource(pd.DataFrame(source.data)[mask])
                priorities[i] += 3 # all it changes is the line width                    
                
        line = p.line('times', 'alt', color=color, source=source, name=''.format(names[i]),
                      line_width=priorities[i], legend=legend)
            
        tooltipped.append(line)
    
    p.legend.click_policy = 'hide'
    p.legend.location = 'bottom_right'
    p.add_tools(HoverTool(renderers=tooltipped, tooltips=tooltips))
    
    curdoc().add_root(p)
    curdoc().title = 'Visibility plot'
    
    return p
