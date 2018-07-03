"""

"""
import sys
import os, re
sys.path.append(os.path.abspath('../'))

from datetime import datetime, timedelta
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import math
import pygal
from pygal.style import Style
import json
from flask import Flask, request, flash, redirect, render_template, url_for, make_response, get_flashed_messages
from SEDMDb.SedmDb import SedmDB
from SEDMDb.SedmDb_tools import DbTools
from werkzeug.security import check_password_hash
import forms
from forms import *
# from flask.ext.appbuilder.charts.views import DirectByChartView
import flask
import flask_login
from flask import Flask, make_response, jsonify

from wtforms import Form, HiddenField, fields, validators
from astropy.io import fits
import requests as url_requests
from flask_table import Table, Col, BoolCol, DatetimeCol, ButtonCol
import StringIO

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from IPython.display import HTML
import matplotlib.pyplot as plt
import pandas as pd
from bokeh.embed import components

import stats_web
import model


# config
SECRET_KEY = 'secret'
USERNAME = 'admin'
PASSWORD = 'default'

db = SedmDB(host='localhost', dbname='sedmdb')
tools = DbTools(db)

app = Flask(__name__)
app.config.from_object(__name__)
login_manager = flask_login.LoginManager()
login_manager.init_app(app)

config = {
    'path':{
    'path_archive':'/scr2/sedmdrp/redux/',
    'path_phot':'/scr2/sedm/phot/',
    'path_raw' : '/scr2/sedm/raw/'}
}

prev_users = ['SEDM_Admin', 'rsw', 'SEDmCzar']

g_p_a_dict = {}  # {g.id: (g.designator, {p.id: (p.designator, {a.id: a.designator})}  )}
g_p_a = db.execute_sql("SELECT g.id, g.designator, p.id, p.designator, a.id, a.designator FROM groups g, program p, allocation a "
                       "WHERE p.group_id = g.id AND a.program_id = p.id;")
for alloc in g_p_a:
    allo=list(alloc)
    allo[0] = int(allo[0])
    allo[2] = int(allo[2])
    allo[4] = int(allo[4])
    if allo[0] in g_p_a_dict.keys():
        if allo[2] in g_p_a_dict[allo[0]][1].keys():
            g_p_a_dict[allo[0]][1][allo[2]][1][allo[4]] = allo[5]
        else:
            g_p_a_dict[allo[0]][1][allo[2]] = (allo[3], {allo[4]: allo[5]})
    else:
        g_p_a_dict[allo[0]] = (allo[1], {allo[2]: (allo[3], {allo[4]: allo[5]})})


def deg2hour(ra, dec, sep=":"):
    '''
    Returns the ra, dec in hours.
    '''
    from astropy.coordinates import SkyCoord

    if ( type(ra) is str and type(dec) is str ):
        return ra, dec
                        
    c = SkyCoord(ra, dec, frame='icrs', unit='deg')                
    ra = c.ra.to_string(unit=u.hourangle, sep=sep, precision=2, pad=True)
    dec = c.dec.to_string(sep=sep, precision=2, alwayssign=True, pad=True)
                                    
    return str(ra), str(dec)
                                        

def radec_str2rad(_ra_str, _dec_str):
    """

    :param _ra_str: 'H:M:S'
    :param _dec_str: 'D:M:S'
    :return: ra, dec in rad
    """
    # convert to rad:
    _ra = map(float, _ra_str.split(':'))
    _ra = (_ra[0] + _ra[1] / 60.0 + _ra[2] / 3600.0) * np.pi / 12.
    _dec = map(float, _dec_str.split(':'))
    _dec = (_dec[0] + _dec[1] / 60.0 + _dec[2] / 3600.0) * np.pi / 180.

    return _ra, _dec


class User(flask_login.UserMixin):
    pass


@login_manager.user_loader
def load_user(user_id):
    users = db.get_from_users(['id', 'username'], {'id': user_id})
    if not users:
        return None
    user = User()
    user.id = users[0][0]
    user.name = users[0][1]
    flask_login.login_user(user) # , remember=True)
    get_user_permissions()
    return user


@app.route('/data/<path:filename>')
@flask_login.login_required
def data_static(filename):
    '''
     Get files from the archive
    :param filename:
    :return:
    '''
    _p, _f = os.path.split(filename)

    if _f.startswith('rc') or _f.startswith('finder'):
        return flask.send_from_directory(os.path.join(config['path']['path_phot'], _p), _f)
    else:
        return flask.send_from_directory(os.path.join(config['path']['path_archive'], _p), _f)


def file_exists(filename):
    '''
    Checks in the local directory to see if the file within the main reduction directory exists.
    '''

    _p, _f = os.path.split(filename)
    return os.path.isfile( os.path.join(config['path']['path_archive'], _p, _f))

@flask_login.login_required
def create_gallery(imagelist, mydate, ncols=6, width=80, filetype='spec'):
    '''
    Creates the right formatting for the image gallery.
    '''

    nrows = int(math.ceil(len(imagelist)*1. / ncols)) 


    s = ""  
    for i in range(nrows):
        s = s + '<div class="row"> \n'
        for j in range(ncols):
            pos =i*ncols + j
            if ( pos < len(imagelist)):

                if filetype=='spec':
                    path = os.path.join(mydate, imagelist[pos])
                #Removing the part of the path which is redundant with the photometric folder.
                elif filetype == 'phot':
                    path = imagelist[pos].replace(config["path"]["path_phot"], "")
                elif filetype == 'finder':
                    path = imagelist[pos].replace(config["path"]["path_phot"], "")
                else:
                    print ("File type not found.")
                impath = flask.url_for('data_static', filename=path)
                print (impath)

                imdir, imname = os.path.split(path)
                imname = imname.split(".")[0]

                if filetype=='spec' and file_exists(path.replace(".png", ".pdf")):
                    impathlink = flask.url_for('data_static', filename=path.replace(".png", ".pdf"))
                else:
                    impathlink = impath

                s = s + '''
                  <div class="col-md-{0}">
                    <div class="thumbnail">
                       <div class="desc">{1}</div>
                      <a href="{2}">
                        <img src="{3}" style="width:{4}% height:"{5}%">
                      </a>
                    </div>
                  </div>\n'''.format(12/ncols, imname, impathlink, impath, width, width)
            else:
                s = s + '</div> \n'
                break
        s = s + '</div> \n'


    return s

@app.route('/data_access/<path:instrument>', methods=['GET'])
@flask_login.login_required
def data_access(instrument):
    '''
    Displays the data access page. 
    It accepts the "date" parameter as the day we want to see the data from.
    With no parameters, it just searches for the last day that the telescope was open.

    :return:

    '''
    gallerydic = {}
    files_table = [("", "")]
    gallery = ""
    gallery_phot = ""
    message = ""

    #Parse the date
    if 'date' not in flask.request.args:
        curdate = datetime.datetime.utcnow()
        mydate = "%d%02d%02d"%(curdate.year, curdate.month, curdate.day)
    else:
        mydate_in = flask.request.args['date']

        #Just making sure that we have only allowed digits in the date
        mydate = re.findall(r"(2\d{3}[0-1]\d{1}[0-3]\d{1})", mydate_in)
        if len(mydate) ==0:
            message = "Incorrect format for the date! Your input is: %s. Shall be YYYYMMDD. \n"%(mydate_in) 
            script, div = "", ""
            mydate = ""
        else:
            mydate = mydate[0]
            message = ""



    if instrument.lower() =='ifu':

        # get the weather stats
        if ( not 'date' in flask.request.args):
            reduxfiles, mydate = model.search_redux_files(None)
        else:
            reduxfiles, mydate = model.search_redux_files(mydate)

        if ( not 'date' in flask.request.args and reduxfiles is None):
            message=message + " No data found up to 100 days prior to today... Weather has been terrible lately!"
        elif ( not 'date' in flask.request.args and not reduxfiles is None):
            message = " Data reduction for the last opened day %s. \r To see a different date, type in the navigation bar: ?date=YYYYDDMM"%mydate
        
        elif ( 'date' in flask.request.args and reduxfiles is None):
            message=message + "No data found for the date %s."%(mydate)
        else:
            message = message + "Reduced files found for %s"%( mydate)

        if not reduxfiles is None:

            spec_files = np.array([i[0].endswith(".txt") for i in reduxfiles.values ])
            files_table = [("/data/%s/%s"%(mydate,i[0]), i[0]) for i in reduxfiles[spec_files].values]
            files_table.sort()
            images = [i[0] for i in reduxfiles.values if i[0].endswith(".png")]
            images.sort()
            gallery = create_gallery(images, mydate, ncols=2, width=100, filetype='spec')


    if instrument.lower() =='rc':

        photfiles = model.search_phot_files_by_imtype(mydate = mydate)

        if ( not 'date' in flask.request.args and photfiles is None):
            message=message + " No data found up to 100 days prior to today... Weather has been terrible lately!"
        elif ( not 'date' in flask.request.args and not photfiles is None):
            message = " Data reduction for the last opened day %s. \r To see a different date, type in the navigation bar: ?date=YYYYDDMM"%mydate
        
        elif ( 'date' in flask.request.args and photfiles is None):
            message=message + "No data found for the date %s."%(mydate)
        else:
            message = message + "Reduced files found for %s"%( mydate)

        if not photfiles is None:
            #images = [i[0] for i in photfiles.values if i[0].endswith(".png")]
            #gallery = create_gallery(images, mydate, ncols=6, width=150, filetype='spec')


            for k in photfiles.keys():
                if k.lower() != 'guider':
                    fs = photfiles[k]
                    fs.sort()
                    gallerydic[k] = create_gallery(fs, mydate, ncols=6, width=150, filetype='phot')    

    if instrument.lower() =='finders':

        photfiles, mydate = model.search_finder_files(mydate=mydate)

        #If there are no files on that date, search for them 100 days back.
        if photfiles is None:
            photfiles, mydate = model.search_finder_files()

        if ( not 'date' in flask.request.args and photfiles is None):
            message=message + " No finder charts found up to 100 days prior to today... Weather has been terrible lately!"
        elif ( not 'date' in flask.request.args and not photfiles is None):
            message = " Data reduction for the last opened day %s. \r To see a different date, type in the navigation bar: ?date=YYYYDDMM"%mydate
        
        elif ( 'date' in flask.request.args and photfiles is None):
            message=message + "No finder charts found for the date %s."%(mydate)
        else:
            message = message + "Finder charts found for %s"%( mydate)

        if not photfiles is None:
            fs = list(photfiles['filename'])
            fs.sort()
            print (fs)
            gallery = create_gallery(fs, mydate, ncols=4, width=150, filetype='finder')   



    return render_template('header.html', current_user=flask_login.current_user) + \
            render_template('data4date.html', data=files_table, gallery=gallery, gallerydic=gallerydic, instrument=instrument.upper(), message=message) + \
            render_template('footer.html')
    

@app.route('/weather_stats', methods=['GET'])
def weather_stats():
    '''
    Displays the weather statistics page. It accepts the "date" parameter as the day we want to see the stats from.
    with no parameters, it jsut searches for the last day that the telescope was open and we have a stats log for it.

    :return:

    '''


    if not 'date' in flask.request.args:
        # get the weather stats
        statsfile, mydate = model.search_stats_file()
        stats_plot = stats_web.plot_stats(statsfile, mydate)
        if (stats_plot is None):
            message=message + " No statistics log found up to 100 days prior to today... Weather has been terrible lately!"
            script, div = None, None
        else:
            message = " Weather statistics for last opened day: %s \r To try a different date, type in the navigation bar: ?date=YYYYDDMM"%(os.path.basename(os.path.dirname(os.path.dirname(statsfile))))
            script, div = components(stats_plot)
    else:
        mydate_in = flask.request.args['date']
        #Just making sure that we have only allowed digits in the date
        mydate = re.findall(r"(2\d{3}[0-1]\d{1}[0-3]\d{1})", mydate_in)
        if len(mydate) ==0:
            message = "Incorrect format for the date! Your input is: %s. Shall be YYYYMMDD. \n"%(mydate_in) 
            script, div = "", ""
        else:
            mydate = mydate[0]
            message = ""

            statsfile, mydate_out = model.search_stats_file(mydate)
            stats_plot = stats_web.plot_stats(statsfile, mydate)
            if (not statsfile):
                message=message + "No statistics log found for the date %s. Showing P18 data."%(mydate)
                script, div = components(stats_plot)

            else:
                stats_plot = stats_web.plot_stats(statsfile, mydate)
                message = message + "Weather statistics for selected day: %s"%(mydate)
                script, div = components(stats_plot)


    return render_template('header.html', current_user=flask_login.current_user) + \
            render_template('weather_stats.html', script=script, div=div, message=message) + \
            render_template('footer.html')

@app.route('/')
def index():
    sys.stdout.flush()  # send any stdout to the logfile
    
    if flask_login.current_user.is_authenticated:
        # retrieve all of the requests submitted by the user
        enddate = datetime.datetime.utcnow() + datetime.timedelta(days=1)
        inidate = datetime.datetime.utcnow() - datetime.timedelta(days=7, hours=8)

        dfreq = model.get_requests_for_user(flask_login.current_user.id, inidate, enddate)


        # organize requests into dataframes by whether they are completed or not
        complete = dfreq[(dfreq['status']=='COMPLETED') | (dfreq['status']=='REDUCED')]
        active = dfreq[(dfreq['status']=='PENDING') | (dfreq['status']=='ACTIVE')]
        expired = dfreq[(dfreq['status']=='EXPIRED')]

        # retrieve information about the user's allocations
        ac = model.get_allocations_user(flask_login.current_user.id)


        # generate the html and table titles
        request_tables = [HTML(active.to_html(escape=False, classes='table', index=False)), HTML(complete.to_html(escape=False, classes='table', index=False)), HTML(expired.to_html(escape=False, classes='table', index=False))]
        request_titles = ['', 'Active Requests for the last 7 days', 'Completed Requests in the last 7 days', 'Expired in the last 7 days']
        alloc_table = [HTML(ac.to_html(escape=False, classes='table table-striped', index=False, col_space=10))]
        alloc_titles = ['', 'Your Active Allocations']
        greeting = 'Hello %s!'%flask_login.current_user.name
        myimage = flask.url_for('static', filename='img/smile.jpg')
    else:  # if there is not user, set the lists as empty
        request_tables=[]
        request_titles=[]
        alloc_table = []
        alloc_titles = []
        greeting = 'Hello stranger! Please, log in to see Wonderland.'
        myimage=flask.url_for('static', filename='img/SEDM_ifu.jpg')

            #render_template('weather_stats.html', script=script, div=div, message=message) + \
    return render_template('header.html', current_user=flask_login.current_user) + \
            render_template('index.html', req_tables = request_tables, req_titles=request_titles, all_table = alloc_table, all_titles = alloc_titles, greeting=greeting, myimage=myimage) + \
            render_template('footer.html')

@app.route('/visibility')
def visibility():
    sys.stdout.flush()  # send any stdout to the logfile
    
    if flask_login.current_user.is_authenticated:
        # retrieve all of the requests submitted by the user
        enddate = datetime.datetime.utcnow() + datetime.timedelta(days=1)
        inidate = datetime.datetime.utcnow() - datetime.timedelta(days=7, hours=8)

        dfreq = model.get_requests_for_user(flask_login.current_user.id, inidate, enddate)


        # organize requests into dataframes by whether they are completed or not
        complete = dfreq[(dfreq['status']=='COMPLETED') | (dfreq['status']=='REDUCED')]
        active = dfreq[(dfreq['status']=='PENDING') | (dfreq['status']=='ACTIVE')]

        # retrieve information about the user's allocations
        ac = model.get_allocations_user(flask_login.current_user.id)


        # generate the html and table titles
        request_tables = [HTML(active.to_html(escape=False, classes='table', index=False))]
        request_titles = ['Active Requests for the last 7 days']


        visibility = plot_visibilities(active['RA'], active['DEC'], active['object'], active['allocation'], active['priority'])
    else:  # if there is not user, set the lists as empty
        return redirect(flask.url_for('index'))

            #render_template('weather_stats.html', script=script, div=div, message=message) + \
    return render_template('header.html', current_user=flask_login.current_user) + \
            render_template('requests_visibility.html', req_tables = request_tables, req_titles=request_titles, visibility=visibility) + \
            render_template('footer.html')

@app.route('/login', methods=['GET', 'POST'])
def login():
    form = LoginForm()
    if form.validate_on_submit():
        # Login and validate the user.index
        # user should be an instance of your `User` class
        username = form.username.data
        password = form.password.data
        user_pass = db.get_from_users(['username', 'password', 'id'], {'username': username})
        if not user_pass:  # if there is no user with that username
            message = "Incorrect username or password!"
            return (render_template('header.html') +
                    render_template('login.html', form=form, message=message) +
                    render_template('footer.html'))
        elif user_pass[0] == -1:
            message = user_pass[1]
            return (render_template('header.html') +
                    render_template('login.html', form=form, message=message) +
                    render_template('footer.html'))
        elif check_password_hash(user_pass[0][1], flask.request.form['password']):  # check password
            # log user in and set .id = id and .name = username
            user = User()
            user.id = user_pass[0][2]
            user.name = username
            flask_login.login_user(user)
            get_user_permissions()  # generate .groups, .program, and .allocation
            flash("Logged in as %s" % (username,))
            return redirect(flask.url_for('index'))
        else:
            message = "Incorrect username or password!"
            return (render_template('header.html') +
                    render_template('login.html', form=form, message=message) +
                    render_template('footer.html'))
    return (render_template('header.html') +#, current_user=flask_login.current_user) +
            render_template('login.html', form=form, message=None) +
            render_template('footer.html'))


@app.route('/manage_user', methods=['GET', 'POST'])
@flask_login.login_required
def manage_user():

    if flask_login.current_user.name != 'SEDM_admin':
        return redirect(flask.url_for('index'))

    message = ""
    form1 = SearchUserForm()
    form2 = UsersForm()

    old_groups = []
    new_groups = []
    allocations =[]

    #Case with no user at all
    if len(flask.request.args) ==0:
        message = "Introduce the search criteria for your user. For exact search try \"."


        return (render_template('header.html') +
                render_template('manage_users.html', form1=form1, from2=form2, message=message) +
                render_template('footer.html'))

    #Case with when we want to search for a specific user
    if 'search_user' in flask.request.form:
        username = form1.search_string.data
    elif('user' in flask.request.args):    
        username = flask.request.args['user']
    elif (form2.username.data):
        username = form2.username.data
        #username.replace("'", "").replace('"', '')
        #username = '"{0}"'.format(username)
    else:
        username = ''

    u = model.get_info_user(username)
    message = u["message"]

    print message

    if username =="" or not "username" in u.keys():

        form2.old_groups.choices = []
        form2.new_groups.choices = []

        if 'add_user' in flask.request.form:

            name = form2.name.data
            email = form2.email.data
            new_password = form2.pass_new.data
            new_password_conf = form2.pass_conf.data


            if form2.pass_new.data and new_password ==new_password_conf:
                status, mes = db.add_user({"username":username, "name":name, "email":email, "password":new_password})
                if status ==0:
                    flash("User created")
                else:
                    message = mes
            else:
                message = "New user requires a password!"

            return (render_template('header.html') +
                        render_template('manage_users.html', form1=form1, form2=form2, allocations=[], message=message) +
                        render_template('footer.html'))

        else:
            return (render_template('header.html') +
                        render_template('manage_users.html', form1=form1, form2=form2, allocations=[], message=message) +
                        render_template('footer.html'))
    else:
        form2.old_groups.choices = [(g[0], g[0]) for g in u["old_groups"]]
        form2.new_groups.choices = [(g[0], g[0]) for g in u["new_groups"]]
        allocations = u["allocations"]


    if 'search_user' in flask.request.form and form1.validate_on_submit():
        form2.username.data = u["username"]
        form2.name.data = u["name"]
        form2.email.data = u["email"]

    elif 'add_group' in flask.request.form :
        username = form2.username.data
        #username.replace("'", "").replace('"', '')
        u = model.get_info_user(username)#'"{0}"'.format(username))
        message = u["message"]

        g = flask.request.form['new_groups']
        model.add_group(u["id"], g)
        message = "Added group for user %s"%(form2.username.data)

    elif 'remove_group' in flask.request.form:

        username = form2.username.data
        #username.replace("'", "").replace('"', '')
        u = model.get_info_user(username)#'"{0}"'.format(username))
        message = u["message"]
        g = flask.request.form['old_groups']
        model.remove_group(u["id"], g)
        message = "Deleted group for user %s"%form2.username.data

    elif 'modify_user' in flask.request.form and form2.validate_on_submit():
        username = form2.username.data
        u = model.get_info_user(username) #'"{0}"'.format(form2.username.data))
        message = u["message"]

        name = form2.name.data
        email = form2.email.data
        new_password = form2.pass_new.data
        new_password_conf = form2.pass_conf.data

        status, mes = db.update_user({'id': u["id"], 'name':name, 'email':email})  
        flash("User with username %s updated with name %s, email %s. %s"%(username, name, email, mes))   

        #If there is any infoirmation in the password field, we update
        if form2.pass_new.data:
            db.update_user({'id': u["id"], 'password': new_password})
            flash("Password changed ")
    elif 'delete_user' in flask.request.form and form2.name.data:

        username = form2.username.data
        u = model.get_info_user(username)

        if 'username' in u.keys():
            status, mes = db.remove_user({'id': u["id"]})  
            flash("Deleted user with username %s. %s"%(username, mes))   
            return (render_template('header.html') +
                    render_template('manage_users.html', form1=form1, form2=None, allocations=[], message=message) +
                    render_template('footer.html'))

    else:
        print "NOTHING TO BE DONE"
        pass

    u = model.get_info_user(username)
    form2.old_groups.choices = [(g[0], g[0]) for g in u["old_groups"]]
    form2.new_groups.choices = [(g[0], g[0]) for g in u["new_groups"]]
    allocations = u["allocations"]

    return (render_template('header.html') +
                    render_template('manage_users.html', form1=form1, form2=form2, allocations=allocations, message=message) +
                    render_template('footer.html'))


@app.route("/delete_allocation", methods=['GET', 'POST'])
@flask_login.login_required

def delete_allocation():
    '''
    This handles when a user needs to cancel a reservation. 
    '''

    if flask_login.current_user.name != 'SEDM_admin':
        return redirect(flask.url_for('index'))

    id = int(request.args.get('id'))
    status, message = model.delete_allocation(id)

    return redirect(url_for('manage_allocation'))

@app.route("/manage_allocation", methods=['GET', 'POST'])
@flask_login.login_required
def manage_allocation():
    '''
    This page allows to add and delete program allocations.
    '''
    if flask_login.current_user.name != 'SEDM_admin':
        return redirect(flask.url_for('index'))


    message = ""
    # Declare the table
    class AllocationTable(Table):
        id = 0
        classes = ['table']
        name = Col('Name')
        program = Col('Program')
        description = Col('Description')
        inidate = DatetimeCol('Initial Date')
        enddate = DatetimeCol('End Date')
        time_allocated = Col('Time allocated')
        time_spent = Col('Time spent')
        active = BoolCol('Active')
        delete = ButtonCol('Delete', 'delete_allocation', url_kwargs=dict(id='id'), anchor_attrs={'class': 'btn btn-danger'})

    # Get some objects
    class Item(object):
        def __init__(self, id, name, program, description, inidate, enddate, time_allocated, time_spent, active):
            self.id = id
            self.name = name
            self.program = program
            self.description = description
            self.inidate = inidate
            self.enddate = enddate
            self.time_allocated = time_allocated
            self.time_spent = time_spent
            self.active = active

    form = AllocationForm()

    if request.method == 'GET':
        form.program_id.choices = [(row["id"], row["name"]) for row in model.get_all_programs()] 

    postform= flask.request.form

    if request.method == 'POST':

        pardic= {  'program_id': int(postform['program_id']),
                   'designator': postform['designator'],
                   'inidate' :  datetime.datetime.strptime(postform['inidate'], "%Y-%m-%d %H:%M:%S"),
                    'enddate' : datetime.datetime.strptime(postform['enddate'], "%Y-%m-%d %H:%M:%S"),
                    'time_allocated': float(postform['time_allocated'])*3600*24,
                    'time_spent': float(postform['time_spent'])*3600*24}
        status, message = db.add_allocation(pardic)

    # Or, more likely, load items from your database with something like
    # Populate the table
    items = model.get_allocations()
    table = AllocationTable(items)


    return (render_template('header.html', current_user=flask_login.current_user) +
            render_template('allocation.html', table=table, form_a = form, message=message) +
            render_template('footer.html'))
    


@app.route("/delete_program", methods=['GET', 'POST'])
@flask_login.login_required

def delete_program():
    '''
    This handles when a user needs to cancel a reservation. 
    '''

    if flask_login.current_user.name != 'SEDM_admin':
        return redirect(flask.url_for('index'))

    id = int(request.args.get('id'))
    status, message = model.delete_program(id)

    return redirect(url_for('manage_program'))

@app.route("/manage_program", methods=['GET', 'POST'])
@flask_login.login_required
def manage_program():
    '''
    This page allows to add and delete programs.
    '''
    if flask_login.current_user.name != 'SEDM_admin':
        return redirect(flask.url_for('index'))


    message = ""
    # Declare the table
    class ProgramTable(Table):
        id = 0
        classes = ['table']
        designator = Col('Designator')
        name = Col('Name')
        group = Col('Group')
        pi = Col('PI')
        priority = Col('Priority')
        delete = ButtonCol('Delete', 'delete_program', url_kwargs=dict(id='id'), anchor_attrs={'class': 'btn btn-danger'})

    # Get some objects
    class Item(object):
        def __init__(self, id, designator, name, group, pi, priority):
            self.id = id
            self.designator = designator
            self.name = name
            self.group = group
            self.pi = pi
            self.priority = priority

    form = ProgramForm()


    form.group.choices = [(row["id"], row["designator"]) for row in model.get_all_groups()] 

    postform= flask.request.form

    if request.method == 'POST':

        pardic= {  'group_id': int(postform['group']),
                   'designator': postform['designator'],
                   'name': postform['name'],
                   'PI': postform['pi'],
                   'priority': float(postform['priority']),
                    }
        print pardic
        status, message = db.add_program(pardic)

    # Populate the table
    items = model.get_programs()
    table = ProgramTable(items)


    return (render_template('header.html', current_user=flask_login.current_user) +
            render_template('program.html', table=table, form_a = form, message=message) +
            render_template('footer.html'))
    

@app.route("/delete_group", methods=['GET', 'POST'])
@flask_login.login_required

def delete_group():
    '''
    This handles when a user needs to delete a group. 
    '''

    if flask_login.current_user.name != 'SEDM_admin':
        return redirect(flask.url_for('index'))

    id = int(request.args.get('id'))
    status, message = model.delete_group(id)

    return redirect(url_for('manage_group'))

@app.route("/manage_group", methods=['GET', 'POST'])
@flask_login.login_required
def manage_group():
    '''
    This page allows to add and delete groups.
    '''
    if flask_login.current_user.name != 'SEDM_admin':
        return redirect(flask.url_for('index'))


    message = ""
    # Declare the table
    class GroupTable(Table):
        id = 0
        classes = ['table']
        designator = Col('Designator')
        delete = ButtonCol('Delete', 'delete_group', url_kwargs=dict(id='id'), anchor_attrs={'class': 'btn btn-danger'})

    # Get some objects
    class Item(object):
        def __init__(self, id, designator):
            self.id = id
            self.designator = designator

    form = GroupForm()

    postform= flask.request.form

    if request.method == 'POST':

        pardic= {  
                   'designator': postform['designator']
                    }

        status, message = db.add_group(pardic)

    # Populate the table
    items = model.get_all_groups()
    table = GroupTable(items)


    return (render_template('header.html', current_user=flask_login.current_user) +
            render_template('group.html', table=table, form_a = form, message=message) +
            render_template('footer.html'))

@app.route("/passchange", methods=['GET', 'POST'])
@flask_login.login_required
def password_change():
    form = forms.PassChangeForm()

    if form.validate_on_submit():
        # check for correct password and change if true
        user_id = flask_login.current_user.id
        password = form.password.data
        new_password = form.pass_new.data
        new_password_conf = form.pass_conf.data
        user_pass = db.get_from_users(['username', 'password', 'id'], {'id': user_id})
        if not user_pass:
            form = forms.PassChangeForm()
            message = "Incorrect username or password!"
            return (render_template('header.html') +
                    render_template('passchange.html', form=form, message=message) +
                    render_template('footer.html'))
        elif user_pass[0] == -1:
            form = forms.PassChangeForm()
            message = user_pass[1]
            return (render_template('header.html') +
                    render_template('passchange.html', form=form, message=message) +
                    render_template('footer.html'))
        elif check_password_hash(user_pass[0][1], flask.request.form['password']):
            db.update_user({'id': user_pass[0][2], 'password': new_password})
            flash("Password changed")
            return redirect(flask.url_for('logout'))
        else:
            form = forms.PassChangeForm()
            message = "Incorrect username or password!"
            return (render_template('header.html', current_user=flask_login.current_user) +
                    render_template('passchange.html', form=form, message=message) +
                    render_template('footer.html'))
    return (render_template('header.html', current_user=flask_login.current_user) +
            render_template('passchange.html', form=form, message=None) +
            render_template('footer.html'))


@app.route("/logout")
@flask_login.login_required
def logout():
    flask_login.logout_user()
    return redirect(flask.url_for('login'))


@app.route('/user_info')
@flask_login.login_required
def user_info():
    return (render_template('header.html', current_user=flask_login.current_user) +
            render_template('footer.html'))

@app.route('/view_request', methods=['GET', 'POST'])
@flask_login.login_required
def view_request():

    message = ""
    if flask_login.current_user.name not in prev_users:
        return redirect(flask.url_for('index'))

    request_id = request.args.get('id')
    if not request_id:
        request_id = 20180416223941629
        message = "Assuming default id for the request."

    req, flt = model.get_request_by_id(request_id)
    message = "seeing the request"

    return render_template('view_request.html', req=req, flt=flt, message=message)

@app.route('/update_request', methods=['GET', 'POST'])
@flask_login.login_required
def update_request():
    if flask_login.current_user.name not in prev_users:
        return redirect(flask.url_for('index'))

    if request.is_json:
        print request.json
        ret = model.parse_update(request.json)
    return jsonify(ret)

@app.route('/delete_request', methods=['GET', 'POST'])
@flask_login.login_required
def delete_request():
    if flask_login.current_user.name not in prev_users:
        return redirect(flask.url_for('index'))

    request_id = request.args.get('id')
    if not request_id:
        request_id = 20180416223941629
    ret = model.delete_request_by_id(request_id)

    return redirect('view_request?id=%s' % request_id)

@app.route('/request', methods=['GET', 'POST'])
@flask_login.login_required
def requests():
    """
    generate forms for request creation
    """
    form1 = RequestForm()
    # if form1.validate_on_submit():
    #    print 'form1'
    alloc = db.execute_sql("SELECT id, designator FROM allocation WHERE id IN %s AND active = 't';" % ('(' + str(flask_login.current_user.allocation)[1:-1] + ')',))
    alloc = model.get_allocations_user(flask_login.current_user.id)
    if alloc is None or len(alloc)==0:
        choices = [(0, "You have none active!")]
    else:
        choices = [z for z in zip(alloc['id'], alloc['allocation'])]
    form1.allocation.choices = choices
    user_id = flask_login.current_user.id
    print Time.now()
    if not form1.inidate.data:
        form1.inidate.data = datetime.datetime.today()
    if not form1.enddate.data:    
        form1.enddate.data = datetime.datetime.today() + datetime.timedelta(2)
    if form1.submit_req.data:
        pass
    if form1.submit_req.data and form1.validate_on_submit():
        # TODO: modify following to accomodate new forms
        print "request submitted"
        if request.method == 'POST':            
            # search for target object, create if it doesn't exist
            if form1.obj_ra.data and form1.obj_dec.data and form1.typedesig.data in ['f', 'v']:
                try:
                    ra = float(form1.obj_ra.data)
                    coord=SkyCoord(ra, form1.obj_dec.data, unit='deg')
                except:
                    coord=SkyCoord(form1.obj_ra.data, form1.obj_dec.data, unit=('hourangle', 'deg'))
                # search database for matching name+ra+dec
                objects = db.get_objects_near(coord.ra.value, coord.dec.value, .5, ['name','id'])
                id_list = [obj[1] for obj in objects if obj[0] == form1.obj_name.data]
                # if it isn't in the database, add it
                if not id_list:
                    add_obj = db.add_object({'name': form1.obj_name.data, 'ra': coord.ra.value, 'dec': coord.dec.value, 'typedesig': form1.typedesig.data})
                    print add_obj
                    obj_id = add_obj[0]
                    if obj_id == -1:
                        message = add_obj[1]
                        return (render_template('header.html', current_user=flask_login.current_user) +
                                render_template('request.html', form1=form1, sso_form=None, message=message) +
                                render_template('footer.html'))
                else:
                    obj_id = id_list[0]
                # if it is periodic, make sure there is an entry
                if form1.typedesig.data == 'v':
                    var = db.get_from_periodic(['id'], {'object_id': obj_id})    
                    if not var:
                        period_form = forms.PeriodicForm()
                        if period_form.validate():
                            obj =  db.add_object({'name': form1.obj_name.data, 'typedesig': 'v'})
                            print obj
                            if obj[0] == -1:
                                message = obj[1]
                                return (render_template('header.html', current_user=flask_login.current_user) +
                                        render_template('request.html', form1=form1, sso_form=None, period_form=period_form) +
                                        render_template('footer.html'))
                            else: 
                                obj_id = obj[0]
                                sso= db.add_periodic({'object_id': obj_id, 'mjd0': sso_form.mjd0.data, 'phasedays': sso_form.phasedays.data})
                                print sso
                                if sso[0] == -1:
                                    message = sso[1]
                                    return (render_template('header.html', current_user=flask_login.current_user) +
                                            render_template('request.html', form1=form1, sso_form=None, period_form=period_form) +
                                            render_template('footer.html'))
                    elif var[0] == -1:
                        message = var[1]
                        return (render_template('header.html', current_user=flask_login.current_user) +
                                render_template('request.html', form1=form1, sso_form=None, message=message) +
                                render_template('footer.html'))
            
            else:  # form requires obj_name and typedesig to be provided
                obj = db.get_from_object(['id'], {'name': form1.obj_name.data, 'typedesig': form1.typedesig.data})
                if not obj:
                    # create object (and orbit for SSO)
                    if form1.typedesig.data == 'e':
                        sso_form = forms.HeliocentricEllipticalForm()
                        if sso_form.validate():
                            obj =  db.add_object({'name': form1.obj_name.data, 'typedesig': 'e'})
                            print obj
                            if obj[0] == -1:
                                message = obj[1]
                                return (render_template('header.html', current_user=flask_login.current_user) +
                                        render_template('request.html', form1=form1, sso_form=sso_form, message=message) +
                                        render_template('footer.html'))
                            else: 
                                obj_id = obj[0]
                                sso= db.add_elliptical_heliocentric({'object_id': obj_id, 'inclination': sso_form.i.data, 'longascnode_O': sso_form.O.data, 'perihelion_o': sso_form.o.data, 
                                                                    'a': sso_form.a.data, 'n': sso_form.n.data, 'e': sso_form.e.data, 'M': sso_form.M.data, 'mjdepoch': Time(str(sso_form.E.data)).mjd, 
                                                                    'D': sso_form.D.data, 'M1': sso_form.M1.data, 'M2': sso_form.M2.data})
                                print sso
                                if sso[0] == -1:
                                    message = sso[1]
                                    return (render_template('header.html', current_user=flask_login.current_user) +
                                            render_template('request.html', form1=form1, sso_form=sso_form, message=message) +
                                            render_template('footer.html'))
                        else:
                            return (render_template('header.html', current_user=flask_login.current_user) +
                                    render_template('request.html', form1=form1, sso_form=sso_form, message=message) +
                                    render_template('footer.html'))
                    elif form1.typedesig.data == 'h':
                        sso_form = forms.HeliocentricHyperbolicForm()
                        if sso_form.validate():
                            obj =  db.add_object({'name': form1.obj_name.data, 'typedesig': 'h'})
                            if obj[0] == -1:
                                message = obj[1]
                                return (render_template('header.html', current_user=flask_login.current_user) +
                                        render_template('request.html', form1=form1, sso_form=sso_form, message=message) +
                                        render_template('footer.html'))
                            else: 
                                obj_id = obj[0]
                                sso= db.add_hyperbolic_heliocentric({'object_id': obj_id, 'inclination': sso_form.i.data, 'longascnode_O': sso_form.O.data, 'perihelion_o': sso_form.o.data, 
                                                                    'q': sso_form.q.data, 'e': sso_form.e.data, 'D': sso_form.D.data,  
                                                                    'T': str(sso_form.T.data), 'M1': sso_form.M1.data, 'M2': sso_form.M2.data})
                                print sso
                                if sso[0] == -1:
                                    message = sso[1]
                                    return (render_template('header.html', current_user=flask_login.current_user) +
                                            render_template('request.html', form1=form1, sso_form=sso_form, message=message) +
                                            render_template('footer.html'))
                        else:
                            return (render_template('header.html', current_user=flask_login.current_user) +
                                    render_template('request.html', form1=form1, sso_form=sso_form, message=message) +
                                    render_template('footer.html'))
                    elif form1.typedesig.data == 'p':
                        sso_form = forms.HeliocentricParabolicForm()
                        if sso_form.validate():
                            obj =  db.add_object({'name': form1.obj_name.data, 'typedesig': 'p'})
                            if obj[0] == -1:
                                message = obj[1]
                                return (render_template('header.html', current_user=flask_login.current_user) +
                                        render_template('request.html', form1=form1, sso_form=sso_form, message=message) +
                                        render_template('footer.html'))
                            else: 
                                obj_id = obj[0]
                                sso= db.add_parabolic_heliocentric({'object_id': obj_id, 'inclination': sso_form.i.data, 'longascnode_O': sso_form.O.data, 'perihelion_o': sso_form.o.data, 
                                                                   'q': sso_form.q.data, 'D': sso_form.D.data,  
                                                                   'T': str(sso_form.T.data), 'M1': sso_form.M1.data, 'M2': sso_form.M2.data})
                                print sso
                                if sso[0] == -1:
                                    message = sso[1]
                                    return (render_template('header.html', current_user=flask_login.current_user) +
                                            render_template('request.html', form1=form1, sso_form=sso_form, message=message) +
                                            render_template('footer.html'))
                        else:
                            return (render_template('header.html', current_user=flask_login.current_user) +
                                    render_template('request.html', form1=form1, sso_form=sso_form, message=message) +
                                    render_template('footer.html'))
                    elif form1.typedesig.data == 'E':
                        sso_form = forms.EarthSatelliteForm()
                        if sso_form.validate():
                            obj =  db.add_object({'name': form1.obj_name.data, 'typedesig': 'E'})
                            if obj[0] == -1:
                                message = obj[1]
                                return (render_template('header.html', current_user=flask_login.current_user) +
                                        render_template('request.html', form1=form1, sso_form=sso_form, message=message) +
                                        render_template('footer.html'))
                            else: 
                                obj_id = obj[0]
                                sso= db.add_earth_satellite({'object_id': obj_id, 'inclination': sso_form.i.data, 'ra': sso_form.ra.data, 'pedigree': sso_form.P.data, 
                                                            'M': sso_form.M.data, 'e': sso_form.e.data, 'n': sso_form.n.data,  
                                                            'T': str(sso_form.T.data), 'decay': sso_form.d.data, 'reforbit': sso_form.r.data})
                                print sso
                                if sso[0] == -1:
                                    message = sso[1]
                                    return (render_template('header.html', current_user=flask_login.current_user) +
                                            render_template('request.html', form1=form1, sso_form=sso_form, message=message) +
                                            render_template('footer.html'))
                        else:
                            return (render_template('header.html', current_user=flask_login.current_user) +
                                    render_template('request.html', form1=form1, sso_form=sso_form, message=message) +
                                    render_template('footer.html'))
                    elif form1.typedesig.data == 'v':
                        sso_form = forms.PeriodicForm()
                        if sso_form.validate():
                            obj =  db.add_object({'name': form1.obj_name.data, 'typedesig': 'v'})
                            if obj[0] == -1:
                                message = obj[1]
                                return (render_template('header.html', current_user=flask_login.current_user) +
                                        render_template('request.html', form1=form1, sso_form=sso_form, message=message) +
                                        render_template('footer.html'))
                            else: 
                                obj_id = obj[0]
                                sso= db.add_periodic({'object_id': obj_id, 'mjd0': sso_form.mjd0.data, 'phasedays': sso_form.phasedays.data})
                                print sso
                                if sso[0] == -1:
                                    message = sso[1]
                                    return (render_template('header.html', current_user=flask_login.current_user) +
                                            render_template('request.html', form1=form1, sso_form=sso_form, message=message) +
                                            render_template('footer.html'))
                        else:
                            return (render_template('header.html', current_user=flask_login.current_user) +
                                    render_template('request.html', form1=form1, sso_form=sso_form, message=message) +
                                    render_template('footer.html'))
                    elif form1.typedesig.data == 'f':
                        message = "No fixed object with that name, RA and DEC required"
                        return (render_template('header.html', current_user=flask_login.current_user) +
                                render_template('request.html', form1=form1, sso_form=None, message=message) +
                                render_template('footer.html'))
                elif obj[0] == -1:
                    message = obj[1]                    
                    return (render_template('header.html', current_user=flask_login.current_user) +
                            render_template('request.html', form1=form1, sso_form=None, message=message) +
                            render_template('footer.html'))
                else:
                    obj_id = obj[0][0]
                   
 
            # generate 'nexposures' from ifu and filter inputs
            if 'ifu' in form1 and form1.ifu.data:
                exposures = '{1, '
                obs_seq = '{1ifu,'
                exptime = '{2400,'
            else:
                exposures = '{'
                obs_seq = '{'
                exptime = '{'
            
            myseq = form1.filters_op.data
            myseq = myseq[0:-1]
            filt_seq = np.array([ int(d) for d in myseq.split(",")])

            for b, f in zip(filt_seq, ['u', 'g', 'r', 'i']):
                if b>0:
                    obs_seq += str(b) +f+','
                    exptime += '60,'


            exposures = exposures + form1.filters_op.data
            obs_seq = obs_seq[0:-1] + '}'
            exptime = exptime[0:-1] + '}'

            # generate the input dictionary for add_request
            pardic = {'object_id': obj_id, 'exptime': exptime,
                      'priority': form1.priority.data, 'inidate': str(form1.inidate.data),
                      'enddate': str(form1.enddate.data), 'user_id': user_id, 'allocation_id': form1.allocation.data,
                      'nexposures': exposures, 'obs_seq': obs_seq}
            # add in the optional parameters
            if form1.cadence.data:
                pardic['cadence'] = form1.cadence.data
            if form1.phasesamples.data:
                pardic['phasesamples'] = form1.phasesamples.data
            if form1.min_moon_dist.data:
                pardic['min_moon_dist'] = form1.min_moon_dist.data
            if form1.max_moon_illum.data:
                pardic['max_moon_illum'] = form1.max_moon_illum.data
            if form1.max_cloud_cover.data:
                pardic['max_cloud_cover'] = form1.max_cloud_cover.data
            req = db.add_request(pardic)
            # TODO: set up user_id and program_id=
            message = req[1]
            return (render_template('header.html', current_user=flask_login.current_user) +
                    render_template('request.html', form1=form1, sso_form=None, message=message) +
                    render_template('footer.html'))
    return (render_template('header.html') +#, current_user=flask_login.current_user) +
            render_template('request.html', form1=form1, sso_form=None, message=False) +
            render_template('footer.html'))
# TODO: set up ability to view and cancel requests
# TODO: allow searching for objects, return link to page with information/observations of it
# TODO: End of Night report as html


@app.route('/schedule')
def schedule():
    # TODO: create a table
    schedule = None  # TODO: decide how schedule will be stored (file/database?)
    return (render_template('header.html') +#, current_user=flask_login.current_user) +
            render_template('schedule.html') +
            render_template('footer.html'))

@app.route('/add_target', methods=['GET','POST'])
def add_json_target():
    import model


    try:
        textfields = ['jsonfile']
        fields = textfields[:]
        
        jsonfile = request.get_json(silent=True)
        if jsonfile:
            jsonfile = form['jsonfile']
            try:
                entryUpdate = json.loads(jsonfile.file.read())
                message = model.updateFollowupConfig(entryUpdate)
            except Exception, e:
                err = open('/home/sedm/kpy/flask/static/error_output%s.txt' % datetime.datetime.utcnow().strftime("%H_%M_%S"),'w')
                err.write(str(e))
                err.close()
        else:
            message = "30 No json file submitted."
        return message
    except Exception, e:
        err = open('/home/sedm/kpy/flask/static/error.log','a')
        err.write('%s: %s\n' % (datetime.datetime.utcnow().strftime("%H_%M_%S"),str(e)))
        err.close()




@app.route('/objects', methods=['GET', 'POST'])
@flask_login.login_required
def objects():
    form1 = SubmitObjectForm()
    form2 = FindObjectForm()
    form2.typedesig.data = 'f'
    message = None
    urls = []
    # add objects
    if form1.add_obj.data and form1.validate_on_submit():
        # TODO: modify following to accomodate new forms
        if request.method == 'POST':            
            # if the object is non-fixed, show the proper forms
            if form1.typedesig.data == 'v':
                period_form = forms.PeriodicForm()
                if not period_form.validate():
                    return (render_template('header.html', current_user=flask_login.current_user) +
                                render_template('object.html', form1=form1, form2=form2, sso_form=None, period_form=period_form) +
                                render_template('footer.html'))
            elif form1.typedesig.data == 'e':
                sso_form = forms.HeliocentricEllipticalForm()
                if not sso_form.validate():
                    return (render_template('header.html', current_user=flask_login.current_user) +
                            render_template('object.html', form1=form1, form2=form2, sso_form=sso_form, message=message) +
                            render_template('footer.html'))
            elif form1.typedesig.data == 'h':
                sso_form = forms.HeliocentricHyperbolicForm()
                if not sso_form.validate():
                    return (render_template('header.html', current_user=flask_login.current_user) +
                            render_template('object.html', form1=form1, form2=form2, sso_form=sso_form, message=message) +
                            render_template('footer.html'))
            elif form1.typedesig.data == 'p':
                sso_form = forms.HeliocentricParabolicForm()
                if not sso_form.validate():
                    return (render_template('header.html', current_user=flask_login.current_user) +
                            render_template('object.html', form1=form1, form2=form2, sso_form=sso_form, message=message) +
                            render_template('footer.html'))
            elif form1.typedesig.data == 'E':
                sso_form = forms.EarthSatelliteForm()
                if not sso_form.validate():
                    return (render_template('header.html', current_user=flask_login.current_user) +
                            render_template('object.html', form1=form1, form2=form2, sso_form=sso_form, message=message) +
                            render_template('footer.html'))
            elif form1.typedesig.data == 'f':
                form1.obj_ra.validators = [validators.input_required()]
                form1.obj_dec.validators = [validators.input_required()]
                if not form1.validate():
                    return (render_template('header.html', current_user=flask_login.current_user) +
                            render_template('object.html', form1=form1, form2=form2, ssoform=None, message=message) +
                            render_template('footer.html'))
                    
            # keep track of whether the object has been found/created
            obj_id = None  
            # search for target object, create if it doesn't exist
            if form1.obj_ra.data and form1.obj_dec.data and form1.typedesig.data in ['f', 'v']:
                # add the base object
                try:
                    ra = float(form1.obj_ra.data.strip())
                    dec = float(form1.obj_dec.data.strip())
                    coord=SkyCoord(ra, dec, unit='deg')
                except:
                    coord=SkyCoord(form1.obj_ra.data, form1.obj_dec.data, unit=('hourangle', 'deg'))
                
                add_obj = db.add_object({'name': form1.obj_name.data, 'ra': coord.ra.value, 'dec': coord.dec.value, 'typedesig': form1.typedesig.data})
                obj_id = add_obj[0]
                if obj_id == -1:
                    message = add_obj[1]
                    return (render_template('header.html', current_user=flask_login.current_user) +
                            render_template('object.html', form1=form1, form2=form2, sso_form=None, message=message) +
                            render_template('footer.html'))
                #else:
                #    obj_id = id_list[0]
            # if it is periodic, make sure there is an entry
            if form1.typedesig.data == 'v':
                # form is already validated above
                if obj_id:  # if the obj_id is already known
                    per= db.add_periodic({'object_id': obj_id, 'mjd0': period_form.mjd0.data, 'phasedays': period_form.phasedays.data})
                    if per[0] == -1:
                        message = per[1]
                        return (render_template('header.html', current_user=flask_login.current_user) +
                                render_template('object.html', form1=form1, form2=form2, sso_form=None, period_form=period_form, message=message) +
                                render_template('footer.html'))
                    else:
                        message = "Periodic entry '%s' added with object_id: %s and periodic_id: %s" % (form1.object_name.data, obj_id, per[0])
                        # reset the forms as an additional indication of completion
                        form1 = SubmitObjectForm()
                        period_form = forms.PeriodicForm()
                        return (render_template('header.html', current_user=flask_login.current_user) +
                                render_template('object.html', form1=form1, form2=form2, sso_form=None, period_form=period_form, message=message) +
                                render_template('footer.html'))
                    
                obj =  db.add_object({'name': form1.obj_name.data, 'typedesig': 'v'})
                if obj[0] == -1:
                    message = obj[1]
                    return (render_template('header.html', current_user=flask_login.current_user) +
                            render_template('object.html', form1=form1, form2=form2, sso_form=None, period_form=period_form, message=message) +
                            render_template('footer.html'))
                else: 
                    obj_id = obj[0]
                    sso= db.add_periodic({'object_id': obj_id, 'mjd0': sso_form.mjd0.data, 'phasedays': sso_form.phasedays.data})
                    if sso[0] == -1:
                        message = sso[1]
                        return (render_template('header.html', current_user=flask_login.current_user) +
                                render_template('object.html', form1=form1, form2=form2, sso_form=None, period_form=period_form) +
                                render_template('footer.html'))
                    
                    message = "Periodic entry for '%s' added with object_id: %s and periodic_id: %s" % (form1.object_name.data, obj_id, per[0])
                    # reset the forms as an additional indication of completion
                    form1 = SubmitObjectForm()
                    period_form = forms.PeriodicForm()
                    return (render_template('header.html', current_user=flask_login.current_user) +
                            render_template('object.html', form1=form1, form2=form2, sso_form=None, period_form=period_form, message=message) +
                            render_template('footer.html'))
            # Handle Solar System Objects
            elif form1.typedesig.data == 'e':
                # form validation checked above
                # create base object
                obj =  db.add_object({'name': form1.obj_name.data, 'typedesig': 'e'})
                if obj[0] == -1:
                    message = obj[1]
                    return (render_template('header.html', current_user=flask_login.current_user) +
                            render_template('object.html', form1=form1, form2=form2, sso_form=sso_form, message=message) +
                            render_template('footer.html'))
                else:  #create orbit entry
                    obj_id = obj[0]
                    sso= db.add_elliptical_heliocentric({'object_id': obj_id, 'inclination': sso_form.i.data, 'longascnode_O': sso_form.O.data, 'perihelion_o': sso_form.o.data, 
                                                        'a': sso_form.a.data, 'n': sso_form.n.data, 'e': sso_form.e.data, 'M': sso_form.M.data, 'mjdepoch': Time(str(sso_form.E.data)).mjd, 
                                                        'D': sso_form.D.data, 'M1': sso_form.M1.data, 'M2': sso_form.M2.data})
                    if sso[0] == -1:
                        message = sso[1]
                        return (render_template('header.html', current_user=flask_login.current_user) +
                                render_template('object.html', form1=form1, form2=form2, sso_form=sso_form, message=message) +
                                render_template('footer.html'))
                    else:
                        message = "SSO entry for '%s' added with object_id: %s and heliocentric_elliptical id: %s" % (form1.object_name.data, obj_id, sso[0])
                        # reset the forms as an additional indication of completion
                        form1 = SubmitObjectForm()
                        sso_form = forms.HeliocentricEllipticalForm()
                        return (render_template('header.html', current_user=flask_login.current_user) +
                                render_template('object.html', form1=form1, form2=form2, sso_form=sso_form, period_form=period_form, message=message) +
                                render_template('footer.html'))
            elif form1.typedesig.data == 'h':
                # form validation checked above
                # create base object
                obj =  db.add_object({'name': form1.obj_name.data, 'typedesig': 'h'})
                if obj[0] == -1:
                    message = obj[1]
                    return (render_template('header.html', current_user=flask_login.current_user) +
                            render_template('object.html', form1=form1, form2=form2, sso_form=sso_form, message=message) +
                            render_template('footer.html'))
                else: #create orbit entry
                    obj_id = obj[0]
                    sso= db.add_hyperbolic_heliocentric({'object_id': obj_id, 'inclination': sso_form.i.data, 'longascnode_O': sso_form.O.data, 'perihelion_o': sso_form.o.data, 
                                                        'q': sso_form.q.data, 'e': sso_form.e.data, 'D': sso_form.D.data,  
                                                        'T': str(sso_form.T.data), 'M1': sso_form.M1.data, 'M2': sso_form.M2.data})
                    if sso[0] == -1:
                        message = sso[1]
                        return (render_template('header.html', current_user=flask_login.current_user) +
                                render_template('object.html', form1=form1, form2=form2, sso_form=sso_form, message=message) +
                                render_template('footer.html'))
                    else:
                        message = "SSO entry for '%s' added with object_id: %s and heliocentric_hyperbolic id: %s" % (form1.object_name.data, obj_id, sso[0])
                        # reset the forms as an additional indication of completion
                        form1 = SubmitObjectForm()
                        sso_form = forms.HeliocentricHyperbolicForm()
                        return (render_template('header.html', current_user=flask_login.current_user) +
                                render_template('object.html', form1=form1, form2=form2, sso_form=sso_form, period_form=period_form, message=message) +
                                render_template('footer.html'))
            elif form1.typedesig.data == 'p':
                # form validation checked above
                # create base object
                obj =  db.add_object({'name': form1.obj_name.data, 'typedesig': 'p'})
                if obj[0] == -1:
                    message = obj[1]
                    return (render_template('header.html', current_user=flask_login.current_user) +
                            render_template('object.html', form1=form1, form2=form2, sso_form=sso_form, message=message) +
                            render_template('footer.html'))
                else: #create orbit entry
                    obj_id = obj[0]
                    sso= db.add_parabolic_heliocentric({'object_id': obj_id, 'inclination': sso_form.i.data, 'longascnode_O': sso_form.O.data, 'perihelion_o': sso_form.o.data, 
                                                       'q': sso_form.q.data, 'D': sso_form.D.data,  
                                                       'T': str(sso_form.T.data), 'M1': sso_form.M1.data, 'M2': sso_form.M2.data})
                    if sso[0] == -1:
                        message = sso[1]
                        return (render_template('header.html', current_user=flask_login.current_user) +
                                render_template('object.html', form1=form1, form2=form2, sso_form=sso_form, message=message) +
                                render_template('footer.html'))
                    else:
                        message = "SSO entry for '%s' added with object_id: %s and heliocentric_parabolic id: %s" % (form1.object_name.data, obj_id, per[0])
                        # reset the forms as an additional indication of completion
                        form1 = SubmitObjectForm()
                        sso_form = forms.HeliocentricParabolicForm()
                        return (render_template('header.html', current_user=flask_login.current_user) +
                                render_template('object.html', form1=form1, form2=form2, sso_form=sso_form, period_form=period_form, message=message) +
                                render_template('footer.html'))
            elif form1.typedesig.data == 'E':
                # form validation checked above
                # create base object
                obj =  db.add_object({'name': form1.obj_name.data, 'typedesig': 'E'})
                if obj[0] == -1:
                    message = obj[1]
                    return (render_template('header.html', current_user=flask_login.current_user) +
                            render_template('object.html', form1=form1, form2=form2, sso_form=sso_form, message=message) +
                            render_template('footer.html'))
                else: #create orbit entry
                    obj_id = obj[0]
                    sso= db.add_earth_satellite({'object_id': obj_id, 'inclination': sso_form.i.data, 'ra': sso_form.ra.data, 'pedigree': sso_form.P.data, 
                                                'M': sso_form.M.data, 'e': sso_form.e.data, 'n': sso_form.n.data,  
                                                'T': str(sso_form.T.data), 'decay': sso_form.d.data, 'reforbit': sso_form.r.data})
                    if sso[0] == -1:
                        message = sso[1]
                        return (render_template('header.html', current_user=flask_login.current_user) +
                                render_template('object.html', form1=form1, form2=form2, sso_form=sso_form, message=message) +
                                render_template('footer.html'))
                    else:
                        message = "SSO entry for '%s' added with object_id: %s and earth_satellite id: %s" % (form1.object_name.data, obj_id, per[0])
                        # reset the forms as an additional indication of completion
                        form1 = SubmitObjectForm()
                        sso_form = forms.EarthSatelliteForm()
                        return (render_template('header.html', current_user=flask_login.current_user) +
                                render_template('object.html', form1=form1, form2=form2, sso_form=sso_form, period_form=period_form, message=message) +
                                render_template('footer.html'))
            else:
                print 'How did wtforms get broken, I handled all of the typedesig options already'

    # search for objects
    if form2.submit_obj.data and form2.validate_on_submit():
        if request.method == 'POST' and form2.validate():
            if form2.obj_ra.data and form2.obj_dec.data:
                found = db.get_objects_near(form2.obj_ra.data, form2.obj_dec.data, form2.radius.data)
            elif form2.obj_name.data:
                found = db.get_object_id_from_name(form2.obj_name.data)
            else:
                found = None

            if found is None:
                message = "Not enough information provided"
            elif form2.obj_ra.data and form2.obj_dec.data and not found:
                message = "No objects in database within the coordinate range entered"
            elif form2.obj_name.data and not found:
                message = "No object names in database contain %s" % (form2.obj_name.data,)
            elif found[0] == -1:
                message = found[1]
            elif len(found) > 1:
                urls = [(obj[0], obj[1]) for obj in found]
                message = "multiple objects matching request:"
            else:
                print found[0][0]
                print url_for('show_objects', ident=found[0][0])
                return redirect(url_for('show_objects', ident=found[0][0]))
    # TODO: make an objects "homepage" to have as the base render
    form2.radius.data = 5
    return (render_template('header.html') +#, current_user=flask_login.current_user) +
            render_template('object.html', form1=form1, form2=form2, ssoform=None, object=None,
                            message=message, urls = urls) +
            render_template('footer.html'))


def get_panstars_cutout_tag(ra, dec, optionDictionary = {'width' : '150', 'border' : '0'}):

    '''
    ra: degrees
    dec: degrees
    the color cotout of 100 arcsec wide is stored as /tmp/tmp_name.jpg
    '''

    image_index_url_red = 'http://ps1images.stsci.edu/cgi-bin/ps1filenames.py?ra={0}&dec={1}&filters=i'.format(ra, dec)
    image_index_url_green = 'http://ps1images.stsci.edu/cgi-bin/ps1filenames.py?ra={0}&dec={1}&filters=r'.format(ra, dec)
    image_index_url_blue = 'http://ps1images.stsci.edu/cgi-bin/ps1filenames.py?ra={0}&dec={1}&filters=g'.format(ra, dec)

    ix_red_request = url_requests.get(image_index_url_red)
    ix_green_request = url_requests.get(image_index_url_green)
    ix_blue_request = url_requests.get(image_index_url_blue)

    ix_red_data = ix_red_request.text.split('\n')
    ix_red_data[0] = ix_red_data[0].split(' ')
    ix_red_data[1] = ix_red_data[1].split(' ')
    ix_red_filename = ix_red_data[1][ix_red_data[0].index('filename')]

    ix_green_data = ix_green_request.text.split('\n')
    ix_green_data[0] = ix_green_data[0].split(' ')
    ix_green_data[1] = ix_green_data[1].split(' ')
    ix_green_filename = ix_green_data[1][ix_green_data[0].index('filename')]

    ix_blue_data = ix_blue_request.text.split('\n')
    ix_blue_data[0] = ix_blue_data[0].split(' ')
    ix_blue_data[1] = ix_blue_data[1].split(' ')
    ix_blue_filename = ix_blue_data[1][ix_blue_data[0].index('filename')]

    image_url = "http://ps1images.stsci.edu/cgi-bin/fitscut.cgi?red=%s&green=%s&blue=%s&filetypes=stack&auxiliary=data&size=400&ra=%.6f&dec=%.6f" % (ix_red_filename, ix_green_filename, ix_blue_filename, ra, dec)
    image_request = url_requests.get(image_url)

    image_uri = image_request.content.encode('base64').replace('\n', '')

    image_style = ' '
    for option in optionDictionary:
        image_style += '%s="%s" ' % (option, optionDictionary[option])

    image_tag = '<img src="data:image/jpg;base64,%s"' % image_uri
    image_tag += ' %s>\n' % image_style

    link_url = 'ps1images.stsci.edu/cgi-bin/ps1cutouts?pos=%s+%s&filter=color&filter=g&filter=r&filter=i&filter=z&filter=y&filetypes=stack&auxiliary=data&size=240&output_size=0&verbose=0&autoscale=99.500000&catlist=' % (ra, dec)
    return image_tag, link_url


def plot_visibility(ra, dec, name):
    import base64
    import astropy.units as u
    from astropy.time import Time
    from astropy.coordinates import SkyCoord, EarthLocation, AltAz

    # Get the coordinates of the object:
    obj = SkyCoord(float(ra), float(dec), unit="deg")

    ##############################################################################
    # Use `astropy.coordinates.EarthLocation` to provide the location of Palomar


    palomar_mountain = EarthLocation(lon=243.1361*u.deg, lat=33.3558*u.deg, height=1712*u.m)
    utcoffset = -8*u.hour  # Pacific Daylight Time
    time = Time.now() - utcoffset

    ##############################################################################
    # This is helpful since it turns out M33 is barely above the horizon at this
    # time. It's more informative to find M33's airmass over the course of
    # the night.
    #
    # Find the alt,az coordinates of M33 at 100 times evenly spaced between 10pm
    # and 7am EDT:

    midnight = Time(datetime.datetime(time.datetime.year, time.datetime.month, time.datetime.day)) - utcoffset

    fig = Figure()
    ax = fig.add_subplot(1, 1, 1)

    from astropy.coordinates import get_sun
    delta_midnight = np.linspace(-8, 8, 1000)*u.hour
    times_tonight = midnight + delta_midnight
    frame_tonight = AltAz(obstime=times_tonight, location=palomar_mountain)
    sunaltazs_tonight = get_sun(times_tonight).transform_to(frame_tonight)

    from astropy.coordinates import get_moon
    moon_tonight = get_moon(times_tonight)
    moonaltazs_tonight = moon_tonight.transform_to(frame_tonight)

    ##############################################################################
    # Find the alt,az coordinates of M33 at those same times:

    objaltazs_tonight = obj.transform_to(frame_tonight)
    objairmasss_tonight = objaltazs_tonight.secz
    objairmasss_tonight[objairmasss_tonight < 1] = 1
    objairmasss_tonight[objairmasss_tonight > 6] = 6

    ##############################################################################
    # Make a beautiful figure illustrating nighttime and the altitudes of M33 and
    # the Sun over that time:

    ax.plot(delta_midnight, sunaltazs_tonight.alt, color='r', label='Sun')
    ax.plot(delta_midnight, moonaltazs_tonight.alt, color=[0.75]*3, ls='--', label='Moon')
    scat= ax.scatter(delta_midnight, objaltazs_tonight.alt,
                c=objairmasss_tonight, label=name, lw=0, s=8,
                cmap='gist_stern_r')
    ax.fill_between(delta_midnight.to('hr').value, 0, 90,
                     sunaltazs_tonight.alt < -0*u.deg, color='0.5', zorder=0)
    ax.fill_between(delta_midnight.to('hr').value, 0, 90,
                     sunaltazs_tonight.alt < -18*u.deg, color='k', zorder=0)

    ax.legend(loc='upper left')
    ax.set_xlim(-8, 8)
    #ax.set_xticks(np.arange(9)*2 -12)
    ax.set_ylim(0, 90)
    ax.set_xlabel('Hours from PST Midnight')
    ax.set_ylabel('Altitude [deg]')

    fig.colorbar(scat, label="airmass")
    fig.suptitle("Visibility for %s UTC"%midnight)
    canvas = FigureCanvas(fig)
    output = StringIO.StringIO()
    #canvas.print_png(output)
    #response = make_response(output.getvalue())
    #response.mimetype = 'image/png'
    
    Figure.savefig(fig, output, format='png')
    image_url =  base64.encodestring(output.getvalue())

    image_tag = '<img src="data:image/jpg;base64,%s">' % image_url 

    return image_tag

def plot_visibilities(ra, dec, name, allocations, priority):
    import base64
    import astropy.units as u
    from astropy.time import Time
    from astropy.coordinates import SkyCoord, EarthLocation, AltAz
    from astropy.coordinates import get_sun
    from astropy.coordinates import get_moon

    colors = ['b', 'r', 'g', 'orange', 'gray', 'yellow', 'w', 'purple', 'brown', 'pink', 'lime', 'c']
    allocations = np.array(allocations)
    name = np.array(name)
    priority = np.array(priority, dtype=np.float)
    allocset = list(set(allocations))
    colordic = {}
    for i, a in enumerate(allocset):
        colordic[a] = colors[i%len(colors)]

    # Get the coordinates of the object:
    objs = SkyCoord(np.array(ra, dtype=np.float), np.array(dec, dtype=np.float), unit="deg")

    ##############################################################################
    # Use `astropy.coordinates.EarthLocation` to provide the location of Palomar


    palomar_mountain = EarthLocation(lon=243.1361*u.deg, lat=33.3558*u.deg, height=1712*u.m)
    utcoffset = -8*u.hour  # Pacific Daylight Time
    time = Time.now() - utcoffset

    ##############################################################################
    # Find the alt,az coordinates of the targets at 60 times evenly spaced between -7 and 7 hours

    midnight = Time(datetime.datetime(time.datetime.year, time.datetime.month, time.datetime.day)) - utcoffset

    fig = Figure()
    ax = fig.add_subplot(1, 1, 1)

    delta_midnight = np.linspace(-7, 7, 60)*u.hour
    times_tonight = midnight + delta_midnight
    frame_tonight = AltAz(obstime=times_tonight, location=palomar_mountain)
    sunaltazs_tonight = get_sun(times_tonight).transform_to(frame_tonight)


    moon_tonight = get_moon(times_tonight)
    moonaltazs_tonight = moon_tonight.transform_to(frame_tonight)

    ##############################################################################
    # Make a beautiful figure illustrating nighttime and the altitudes of the Sun over that time:

    ax.plot(delta_midnight, sunaltazs_tonight.alt, color='yellow', label='Sun')
    ax.plot(delta_midnight, moonaltazs_tonight.alt, color=[0.75]*3, ls='--', label='Moon')

    ax.fill_between(delta_midnight.to('hr').value, 0, 90,
                     sunaltazs_tonight.alt < -0*u.deg, color='lightgray', zorder=0)
    ax.fill_between(delta_midnight.to('hr').value, 0, 90,
                     sunaltazs_tonight.alt < -18*u.deg, color='darkgray', zorder=0)
    ax.hlines(30, -7, 7, colors="k", linestyles="dashed", lw=0.5)

    ##############################################################################
    # Find the alt,az coordinates of the objects at those same times:
    # The lines will be grouped by program (color) and priority (line thickness)

    for i, obj in enumerate(objs):

        objaltazs_tonight = obj.transform_to(frame_tonight)

        scat= ax.plot(delta_midnight, objaltazs_tonight.alt, label=name[i], lw=0.3+priority[i]*0.33, c=colordic.get(allocations[i], "m"))

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=6)

    #Set the time edges as 5 degrees above the horizon for the Sun
    xmin = delta_midnight.to('hr').value[sunaltazs_tonight.alt < 5*u.deg][0]
    xmax = delta_midnight.to('hr').value[sunaltazs_tonight.alt < 5*u.deg][-1]
    ax.set_xlim(xmin, xmax)

    ax.set_ylim(0, 90)
    ax.set_xlabel('Hours from PST Midnight')
    ax.set_ylabel('Altitude [deg]')

    fig.suptitle("Visibility for %s UTC"%midnight)
    canvas = FigureCanvas(fig)
    output = StringIO.StringIO()
    #canvas.print_png(output)
    #response = make_response(output.getvalue())
    #response.mimetype = 'image/png'
    
    # Return the figure as a string.
    Figure.savefig(fig, output, format='png')
    image_url =  base64.encodestring(output.getvalue())

    image_tag = '<img src="data:image/jpg;base64,%s">' % image_url 

    return image_tag

@app.route('/objects/<path:ident>')
def show_objects(ident):
    # take name or id and show info about it and images with permission
    try:
        iden = int(ident)
    except ValueError:  # if the name is given instead of id, redirect to the id
        name = str(ident)
        iden = db.get_from_object(['id'], {'name': name})
        if not iden:
            flash("Object with name '%s' not found, please search for object" % (ident,))
            return redirect(url_for('objects'))
        elif iden[0] == -1:
            flash("Searching for name '%s' caused error '%', please report issue" % (ident, iden[1]))
            return redirect(url_for('objects'))
        else:
            return redirect(url_for('show_objects', ident=iden[0][0]))
    # get basic information
    info = db.get_from_object(['name', 'ra', 'dec', 'typedesig', 'epoch', 'id'], {'id': iden})[0]
    info = info + deg2hour(info[1], info[2])
    image_tag, link_url = get_panstars_cutout_tag(info[1], info[2])
    visibility = plot_visibility(info[1], info[2], info[0])
    #visibility = ""
    # TODO: work in SSO objects
    if info[0] == -1:
        flash("Searching for id '%s' caused error '%s', please report issue")
        return redirect(url_for('objects'))
    elif not info:
        flash("Searching for id '%s' produced no results, please search for object")
        return redirect(url_for('objects'))
    # don't check for requests if not logged in
    if not flask_login.current_user.is_authenticated:
        return (render_template('header.html') +
                render_template('object_stats.html', info=info, req_data=[]) +
                render_template('footer.html'))


    # retrieve all of the requests submitted for the object
    request_query = ("SELECT u.username, a.designator, r.inidate, r.enddate, r.priority, r.status "
                     "FROM request r, users u, program p, allocation a WHERE r.user_id = u.id AND a.program_id = p.id "
                     "AND a.id = r.allocation_id AND r.object_id = '%s' AND r.allocation_id IN %s;"
                     % (ident, '(' + str(flask_login.current_user.allocation)[1:-1] +')'))
    req = db.execute_sql(request_query)
    # generate a dataframe, filter for requests that weren't canceled or expired
    req_data = pd.DataFrame(req, columns=['Requester', 'Allocation', 'Start Date', 'End Date', 'Priority', 'Status'])
    req_data = req_data.loc[req_data.Status.isin(['PENDING', 'ACTIVE', 'COMPLETED', 'REDUCED'])]
    req_table = [HTML(req_data.to_html(escape=False, classes='complete', index=False))]
    req_titles = ['', 'Requests']

    # get the time, duration and file for observations of the object
    a_query = ("SELECT obs.mjd, obs.exptime, obs.filter, obs.fitsfile "
               "FROM observation obs, object obj, request r WHERE obj.id = obs.object_id AND "
               "r.id = obs.request_id AND obj.id = '%s' AND r.allocation_id IN %s;" % (ident,'(' + str(flask_login.current_user.allocation)[1:-1] +')'))
    observations = db.execute_sql(a_query)
    obs = []

    def make_image(fitsfile):
        data = fits.open(fitsfile)[0].data
        mydiv = plot([Heatmap(z=data, output_type='div')])
        return mydiv

    for ob in observations:
        obs.append(([ob[0], ob[1], ob[2]], None))
        '''if ob[3]:
            day = ob[3][3:11]
            imagepath = os.path.join(config['path']['path_raw'], day, ob[3])
            obs.append(([ob[0], ob[1], ob[2]], make_image(ob[3])))
        else:'''


    return (render_template('header.html') +  # , current_user=flask_login.current_user) +
            render_template('object_stats.html', info=info, observations=obs, req_table=req_table, req_titles=req_titles, image_tag = image_tag, image_url=link_url, visibility=visibility) +
            render_template('footer.html'))


@app.route('/project_stats')
@flask_login.login_required
def projects_home():  # TODO: list groups and projects.
    # generate a list of (group designator, project id, project designator) 

    if flask_login.current_user.name == 'SEDM_admin':
        user_id = None
    else:
        user_id = flask_login.current_user.id

    alloc_stats = model.get_allocation_stats(user_id)
    stats_plot = stats_web.plot_stats_allocation(alloc_stats)
    script, div = components(stats_plot)

    
    '''alloc_stats_week = model.get_allocation_stats(user_id, datetime.datetime.utcnow()-datetime.timedelta(days=1, hours=8), datetime.datetime.utcnow() - datetime.timedelta(hours=8) )
    stats_plot_week = stats_web.plot_stats_allocation(alloc_stats_week)
    script_week, div_week = components(stats_plot_week)'''

    return (render_template('header.html') + 
            render_template('projects_home.html', script=script, div=div) +  
            #render_template('projects_home.html', script=script_week, div=div_week) +  
            render_template('footer.html'))

    

@app.route('/project_stats/<program>')
@flask_login.login_required
def project_stats(program):
    day = datetime.datetime.now()
    start = Time(day - datetime.timedelta(days=day.weekday())).mjd
    end = start + 7
    alloc_id = None


    permission_query = ("SELECT DISTINCT a.id \
                        FROM groups g, usergroups ug, program p, allocation a \
                        WHERE p.group_id = g.id AND ug.group_id = g.id AND a.program_id = p.id \
                        AND ug.user_id = %d AND a.designator='%s';" % ( long(flask_login.current_user.get_id()), program))
                   
    alloc_id = db.execute_sql(permission_query)
    
    if len(alloc_id)==0:
        flash('allocation %s does not exist for user %s' % (program,flask_login.current_user.name))
        return redirect(url_for('index'))
    else:
        alloc_id = int(alloc_id[0][0])
    
    # TODO: make a form to set start/end dates
    time_used_query = ("SELECT o.time_elapsed, r.status FROM request r, observation o "
                       "WHERE o.request_id = r.id AND r.allocation_id = %s"
                       #" AND r.lastmodified < %s AND r.lastmodified > %s"
                       ";" % (alloc_id,))
    requests = db.execute_sql(time_used_query)
    time_allocated = db.get_from_allocation(['time_allocated'], {'id': alloc_id})[0][0]
    
    if requests: # check whether there are any requests
        reqs = np.array(requests)
        pending = (reqs.T[1] == 'PENDING')
        pending_time = sum(np.nan_to_num(reqs.T[0][pending].astype(float)))/3600.
        # TODO: fix status keywords
        observed = ((reqs.T[1] == 'ACTIVE') | (reqs.T[0] == 'COMPLETED') | (reqs.T[0] == 'REDUCED'))
        observed_time = sum(np.nan_to_num(reqs.T[0][observed].astype(float)))/3600.
    else:
        pending_time = 0.
        observed_time = 0.
    if not time_allocated:
        time_allocated = 0.
    else:
        time_allocated = time_allocated.total_seconds()
    # TODO: make total time show up on graph
    # pygal graphing
    # use different color schemes based on whether there is more time requested than allowed
    if observed_time+pending_time > time_allocated:
        c_style = Style(
            background='transparent',
            colors=('#52e852', '#1a4fba', '#e22626')
            )
        pi_chart = pygal.Pie(style=c_style, height=300, width=300)
        pi_chart.title = "Time allocation (h)"
        pi_chart.add('observed', observed_time)
        pi_chart.add('requested', (time_allocated-observed_time))
        pi_chart.add('requested beyond allocation', (observed_time+pending_time-time_allocated))
    else:
        c_style = Style(
            background='transparent',
            colors=('#52e852', '#1a4fba', '#e5f2de')
            )
        pi_chart = pygal.Pie(style=c_style, height=300, width=300)
        pi_chart.title = "Time allocation (h)"
        pi_chart.add('observed', observed_time)
        pi_chart.add('requested', pending_time)
        pi_chart.add('free', (time_allocated-pending_time-observed_time))
        
    # retrieve all of the requests submitted for the project
    request_query = ("SELECT u.name, u.username, o.name, o.id, r.inidate, r.enddate, r.priority, r.status "
                     "FROM request r, users u, object o WHERE r.user_id = u.id AND r.object_id = o.id "
                     "AND r.allocation_id = %s;" % (alloc_id,))
    req = db.execute_sql(request_query)
    for n, x in enumerate(req):
        req[n] = list(x)
        # decide between name/username
        if not [0]:
            del req[n][0]
        else:
            del req[n][1]
        # decide between name/id
        if not req[n][1]:
            del req[n][1]
        else:
            del req[n][2]

    return (render_template('header.html', current_user=flask_login.current_user) +
            render_template('project_stats.html', img_data=pi_chart.render_data_uri(), req_data=req) +
            render_template('footer.html'))





def get_user_permissions():
    """returns ([groups], [programs], [allocations]) allowed for the current user"""

    
    groups = [-1]
    programs = [-1]
    allocations = [-1]

    if int(flask_login.current_user.get_id()) == 2:  # user_id == 2 is admin
        # use the g_p_a_dict to add all groups, programs and allocations
        # add impossible value to prevent empty IN clauses
        groups = [-1]
        programs = [-1]
        allocations = [-1]
        # set current_user values for groups programs and allocations
        flask_login.current_user.groups = groups
        flask_login.current_user.program = programs
        flask_login.current_user.allocation = allocations
        return
        
    '''group_query = ("SELECT g.id FROM groups g, usergroups ug "
                   "WHERE ug.group_id = g.id "
                   "AND ug.user_id = '%s';" % (flask_login.current_user.get_id()))'''
                   
    permission_query = ("SELECT g.id, p.id, a.id, a.designator \
                        FROM groups g, usergroups ug, program p, allocation a \
                        WHERE p.group_id = g.id AND ug.group_id = g.id AND a.program_id = p.id \
                        AND ug.user_id = '%d';" % ( long(flask_login.current_user.get_id())))
                   
    gr_pr_al = db.execute_sql(permission_query)
    
    if (len (gr_pr_al)>0):
        groups = [int(g[0]) for g in gr_pr_al]  # will be a list of the ids
        programs = [int(g[1]) for g in gr_pr_al]  # will be a list of the ids
        allocations = [int(g[2]) for g in gr_pr_al]  # will be a list of the ids

    else:
        # add impossible value to prevent empty IN clauses
        groups = [-1]
        programs = [-1]
        allocations = [-1]
    
    # set current_user values for groups programs and allocation
    flask_login.current_user.groups = groups
    flask_login.current_user.program = programs
    flask_login.current_user.allocation = allocations
    return


@login_manager.unauthorized_handler
def unauthorized_handler():
    return flask.redirect(flask.url_for('login'))

if __name__ == '__main__':
    app.run()


