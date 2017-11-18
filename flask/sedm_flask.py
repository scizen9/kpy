"""

"""
import sys
import os
sys.path.append(os.path.abspath('../'))

import os
import sys
from datetime import datetime, timedelta
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import pygal
from pygal.style import Style
import json
from flask import Flask, request, flash, redirect, render_template, url_for, make_response, get_flashed_messages
from SEDMDb.SedmDb import SedmDB 
from SEDMDb.SedmDb_tools import DbTools
from werkzeug.security import check_password_hash
import forms
from forms import RequestForm, RedirectForm, FindObjectForm, LoginForm, SubmitObjectForm, is_safe_url
# from flask.ext.appbuilder.charts.views import DirectByChartView
import flask
import flask_login
from wtforms import Form, HiddenField, fields, validators
from astropy.io import fits

from IPython.display import HTML
import matplotlib.pyplot as plt
import pandas as pd
import stats_web
from bokeh.embed import components

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

def search_stats_file():
    '''
    Returns the last stats file that is present in the system according to the present date.
    It also returns a message stating what date that was.
    '''
    return "/scr2/sedm/phot/20171101/stats/stats.log"

@app.route('/')
#@flask_login.login_required
def index():
    sys.stdout.flush()  # send any stdout to the logfile

    if flask_login.current_user.is_authenticated:
        # retrieve all of the requests submitted by the user
        request_query = ("SELECT a.designator, o.name, r.inidate, r.enddate, r.priority, r.status "
                 "FROM request r, object o, allocation a WHERE o.id = r.object_id "
                 "AND a.id = r.allocation_id AND r.user_id = '%s';"
                 % (flask_login.current_user.id,))
        req = db.execute_sql(request_query)
        # organize requests into dataframes by whether they are completed or not
        complete = pd.DataFrame([request for request in req if request[5] in ['COMPLETED', 'REDUCED']], columns=['allocation', 'object', 'start date', 'end date', 'priority','status'])
        active = pd.DataFrame([request for request in req if request[5] in ['PENDING', 'ACTIVE']], columns=['allocation', 'object', 'start date', 'end date', 'priority','status'])
        # retrieve information about the user's allocations
        allocations = db.execute_sql("SELECT a.id, a.designator, a.active, p.designator, g.designator FROM allocation a, program p, groups g WHERE "
                                         "a.program_id=p.id AND p.group_id = g.id AND a.id IN %s;" % ('(' + str(flask_login.current_user.allocation)[1:-1] +')',))
        # create the dataframe and set the allocation names to be linked
        data = pd.DataFrame(allocations, columns=['allocation id', 'allocation', 'active', 'program', 'group'])
        data['allocation'] = data['allocation'].apply(lambda x: '<a href="/project_stats/%s">%s</a>' % (x,x))
        active_alloc = data.loc[data.active == True]  # filter for active allocations
        ac = active_alloc.drop(['allocation id', 'active'], axis=1)
        # generate the html and table titles
        request_tables = [HTML(active.to_html(escape=False, classes='active', index=False)), HTML(complete.to_html(escape=False, classes='complete', index=False))]
        request_titles = ['', 'active requests', 'complete requests']
        alloc_table = [HTML(ac.to_html(escape=False, classes='allocations', index=False))]
        alloc_titles = ['', 'active allocations']
    else:  # if there is not user, set the lists as empty
        request_tables=[]
        request_titles=[]
        alloc_table = []
        alloc_titles = []
    # get the weather stats
    statsfile = search_stats_file()
    stats_plot = stats_web.plot_stats(statsfile)
    message = "Weather statistics for last opened day: %s"%(os.path.basename(os.path.dirname(os.path.dirname(statsfile))))
    script, div = components(stats_plot)

    return render_template('header.html', current_user=flask_login.current_user) + \
            render_template('index.html', req_tables = request_tables, req_titles=request_titles, all_table = alloc_table, all_titles = alloc_titles) + \
            render_template('weather_stats.html', script=script, div=div, message=message) + \
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
    if alloc and alloc[0] == -1:  # TODO: test and remove
        flash(alloc[1])
        return redirect(flask.url_for('/index'))
    elif not alloc:
        choices = [(0, "You have none active!")]
    else:
        choices = alloc
    form1.allocation.choices = choices
    user_id = flask_login.current_user.id
    print Time.now()
    if not form1.inidate.data:
        form1.inidate.data = datetime.today()
    if not form1.enddate.data:    
        form1.enddate.data = datetime.today()+timedelta(2)
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
                if 'ab' in form1 and form1.ab.data:
                    exposures = '{2,'
                else:
                    exposures = '{1,'
            else:
                exposures = '{0,'
            exposures += form1.filters_op.data
            # generate the input dictionary for add_request
            pardic = {'object_id': obj_id, 'exptime': '{120, 2400}',
                      'priority': form1.priority.data, 'inidate': str(form1.inidate.data),
                      'enddate': str(form1.enddate.data), 'user_id': user_id, 'allocation_id': form1.allocation.data,
                      'nexposures': exposures}
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


@app.route('/objects', methods=['GET', 'POST'])
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
                    ra = float(form1.obj_ra.data)
                    coord=SkyCoord(ra, form1.obj_dec.data, unit='deg')
                except:
                    coord=SkyCoord(form1.obj_ra.data, form1.obj_dec.data, unit=('hourangle', 'deg'))
                
                add_obj = db.add_object({'name': form1.obj_name.data, 'ra': coord.ra.value, 'dec': coord.dec.value, 'typedesig': form1.typedesig.data})
                obj_id = add_obj[0]
                if obj_id == -1:
                    message = add_obj[1]
                    return (render_template('header.html', current_user=flask_login.current_user) +
                            render_template('object.html', form1=form1, form2=form2, sso_form=None, message=message) +
                            render_template('footer.html'))
                else:
                    obj_id = id_list[0]
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
        if ob[3]:
            obs.append(([ob[0], ob[1], ob[2]], make_image(ob[3])))
        else:
            obs.append(([ob[0], ob[1], ob[2]], None))
    # TODO: make sure the above works

    return (render_template('header.html') +  # , current_user=flask_login.current_user) +
            render_template('object_stats.html', info=info, observations=obs, req_table=req_table, req_titles=req_titles) +
            render_template('footer.html'))


@app.route('/project_stats')
@flask_login.login_required
def projects_home():  # TODO: list groups and projects.
    # generate a list of (group designator, project id, project designator)    
    allocations = db.execute_sql("SELECT a.id, a.designator, a.active, p.designator, g.designator FROM allocation a, program p, groups g WHERE "
                                 "a.program_id=p.id AND p.group_id = g.id AND a.id IN %s;" % ('(' + str(flask_login.current_user.allocation)[1:-1] +')',))
    data = pd.DataFrame(allocations, columns=['allocation id', 'allocation', 'active', 'program', 'group'])
    data['allocation'] = data['allocation'].apply(lambda x: '<a href="/project_stats/%s">%s</a>' % (x,x))
    active = data.loc[data.active == True]
    inactive = data.loc[data.active == False]
    ac = active.drop(['allocation id', 'active'], axis=1)
    inac = inactive.drop(['allocation id', 'active'], axis=1)
    return (render_template('header.html') + 
            render_template('projects_home.html', allocations = allocations, tables=[HTML(ac.to_html(escape=False, classes='active', index=False)), HTML(inac.to_html(escape=False, classes='inactive', index=False))],
                            titles = ['', 'Active allocations', 'Inactive allocations']) + 
            render_template('footer.html'))


@app.route('/project_stats/<program>')
@flask_login.login_required
def project_stats(program):
    day = datetime.now()
    start = Time(day - timedelta(days=day.weekday())).mjd
    end = start + 7
    alloc_id = None
    for lis in g_p_a_dict.itervalues():
        for prog in lis[1].itervalues():
            for i, allo in prog[1].iteritems():
                if allo == program:
                    alloc_id = i
    if not alloc_id:
        flash('program %s does not exist' % (program,))
        return redirect(url_for('index'))
    if alloc_id not in flask_login.current_user.allocation:
        flash('you do not have permissions to view %s' % (program,))
        return redirect(url_for('index'))
    
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
    if int(flask_login.current_user.get_id()) == 2:  # user_id == 2 is admin
        # use the g_p_a_dict to add all groups, programs and allocations
        groups = [group for group in g_p_a_dict.keys()]
        programs = [p for group in groups for p in g_p_a_dict[group][1].keys()]
        allocations = [a for group in groups for p in g_p_a_dict[group][1].values() for a in p[1].keys()]
        # add impossible value to prevent empty IN clauses
        groups.append(-1)
        programs.append(-1)
        allocations.append(-1)
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
                        AND ug.user_id = '%d';" % (flask_login.current_user.get_id()))
                   
    gr_pr_al = db.execute_sql(permission_query)
    
    if (len (gr_pr_al)>0):
        groups = [long(group[0]) for g in gr_pr_al]  # will be a list of the ids
        # use the g_p_a_dict to generate programs and allocations
        #programs = [p for group in groups for p in g_p_a_dict[group][1].keys()]
        #allocations = [a for group in groups for p in g_p_a_dict[group][1].values() for a in p[1].keys()]
        
        programs = [long(g[1]) for g in gr_pr_al]  # will be a list of the ids
        allocations = [long(g[2]) for g in gr_pr_al]  # will be a list of the ids

    else:
        # add impossible value to prevent empty IN clauses
        groups.append(-1)
        programs.append(-1)
        allocations.append(-1)
    
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


