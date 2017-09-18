"""

"""
import sys
import os
sys.path.append(os.path.abspath('../'))

import os
import datetime
from astropy.time import Time
import numpy as np
import pygal
from pygal.style import Style
import json
from flask import Flask, request, flash, redirect, render_template, url_for, make_response, get_flashed_messages
from SEDMDb.SedmDb import SedmDB
from SEDMDb.SedmDb_tools import DbTools
from werkzeug.security import check_password_hash
from forms import RequestForm, RedirectForm, FindObjectForm, LoginForm, SubmitObjectForm, SSOForm, is_safe_url
# from flask.ext.appbuilder.charts.views import DirectByChartView
import flask
import flask_login
from wtforms import Form, HiddenField
from astropy.io import fits
from plotly.offline import plot
from plotly.graph_objs import Heatmap
import matplotlib.pyplot as plt

# config
SECRET_KEY = 'secret'
USERNAME = 'admin'
PASSWORD = 'default'

db = SedmDB()
tools = DbTools(db)

app = Flask(__name__)
app.config.from_object(__name__)
login_manager = flask_login.LoginManager()
login_manager.init_app(app)

groups = db.execute_sql("SELECT id, designator FROM groups;")
group__dict = {}  # make dict with format {group_id: group designator, ...}
for group in groups:
    group__dict[group[0]] = group[1]

programs = db.execute_sql("SELECT group_id, id, designator from program;")
group_programs_dict = {}  # make dict with format {group_id: (group designator, [(program id, program designator), ...]}
for g_p in programs:
    if g_p[0] not in group_programs_dict.keys():
        group_programs_dict[g_p[0]] = [group__dict[g_p[0]], [(g_p[1], g_p[2])]]
    else:
        group_programs_dict[g_p[0]][1].append((g_p[1], g_p[2]))


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
        return
    user = User()
    #user.id = users[0][1]
    #user.name = users[0][1]
    return user


@app.route('/')
def index():
    # TODO: make it pretty
    return render_template('header.html', current_user=flask_login.current_user) + \
           render_template('index.html') + \
           render_template('footer.html')


@app.route('/login', methods=['GET', 'POST'])
def login():
    form = LoginForm()
    if form.validate_on_submit():
        # Login and validate the user.
        # user should be an instance of your `User` class
        username = form.username.data
        password = form.password.data
        user_pass = db.get_from_users(['username', 'password', 'id'], {'username': username})
        if not user_pass:
            message = "Incorrect username or password!"
            return (render_template('header.html') +
                    render_template('login.html', form=form, message=message) +
                    render_template('footer.html'))

        elif user_pass[0] == -1:
            message = user_pass[1]
            return (render_template('header.html') +
                    render_template('login.html', form=form, message=message) +
                    render_template('footer.html'))
        elif check_password_hash(user_pass[0][1], flask.request.form['password']):
            user = User()
            user.id = user_pass[0][2]
            user.name = username
            flask_login.login_user(user, remember=True)
            group_query = ("SELECT groups.id FROM groups "
                           "JOIN usergroups ON usergroups.group_id = groups.id "
                           "WHERE usergroups.user_id = '%s'" % (flask_login.current_user.id))
            gr = db.execute_sql(group_query)
            groups = tuple([group[0] for group in gr])  # will be a tuple of the ids
            flask_login.current_user.groups = groups
            return redirect(flask.url_for('index'))
        else:
            message = "Incorrect username or password!"
            return (render_template('header.html') +
                    render_template('login.html', form=form, message=message) +
                    render_template('footer.html'))
    return (render_template('header.html') +#, current_user=flask_login.current_user) +
            render_template('login.html', form=form, message=None) +
            render_template('footer.html'))

@app.route('/user_info')
@flask_login.login_required
def user_info():
    return (render_template('header.html') +#, current_user=flask_login.current_user) +
            render_template('footer.html'))


@app.route('/request', methods=['GET', 'POST'])
#@flask_login.login_required
def requests():
    """
    generate forms for request creation
    """
    form1 = RequestForm()
    form2 = FindObjectForm()
    # if form1.validate_on_submit():
    #    print 'form1'
    choices = []
    # TODO: make this restricted to the user's programs
    for g_p in group_programs_dict.values():
        for program in g_p[1]:
            choices.append(program)
    print choices
    form1.program.choices = choices
    form1.user_id.data = 1
    if form1.submit_req.data and form1.validate_on_submit():
        print "request submitted"
        if request.method == 'POST':
            filters = '{'
            for entry in form1.filters.entries:
                filters += str(entry.data['ifu_val']) + ', '
                filters += str(entry.data['u_val']) + ', '
                filters += str(entry.data['g_val']) + ', '
                filters += str(entry.data['r_val']) + ', '
                filters += str(entry.data['i_val']) + '}'
                # TODO: handle 'ordering' in form/parse
            print filters
            pardic = {'object_id': form1.object_id.data, 'exptime': '{120, 2400}',
                                  'priority': form1.priority.data, 'inidate': str(form1.inidate.data),
                                  'enddate': str(form1.enddate.data), 'user_id': 1, 'program_id': form1.program.data,
                                  'nexposures': filters}
            if form1.cadence.data:
                pardic['cadence'] = form1.cadence.data
            if form1.phasesamples.data:
                pardic['phasesamples'] = form1.phasesamples.data
            req = db.add_request(pardic)
            # TODO: set up user_id and program_id=
            message = req[1]
            return (render_template('header.html', current_user=flask_login.current_user) +
                    render_template('request.html', form1=form1, form2=form2, message=message) +
                    render_template('footer.html'))

    if form2.submit_obj.data and form2.validate_on_submit():
        print 'object search'
        if request.method == 'POST' and form2.validate():
            if form2.RA.data and form2.DEC.data:
                ra = form2.RA.data
                dec = form2.DEC.data
                rad = form2.radius.data
                req = db.get_objects_near(ra, dec, rad)
                if req:
                    if req[0] == -1:
                        message = req[1]
                        obj_id = None
                    else:
                        message = 'object(s) within %s arcsec of coordinates: ' % (rad,)
                        for target in req:
                            message += '(name: %s, id: %s)' % (target[1], target[0])
                        obj_id = req[0][0]
                else:
                    message = "No objects in db within 5 arcsec of the coordinates!"
                    obj_id = None
            elif form2.object_name.data:
                obj_name = form2.object_name.data
                req = db.get_object_id_from_name(obj_name)
                if req:
                    message = 'name matches object(s): '
                    for target in req:
                        message += '(name: %s, id: %s)' % (target[1], target[0])
                    obj_id = req[0][0]
                else:
                    message = "No objects found with matching name"
                    obj_id = None
            else:
                message = 'To search for an object_id, either RA and Dec, or a partial name is needed'
                obj_id = None

            form1.object_id.data = obj_id  # set the object_id
            return (render_template('header.html') +#, current_user=flask_login.current_user) +
                    render_template('request.html', form1=form1, form2=form2, message=message) +
                    render_template('footer.html'))
    return (render_template('header.html') +#, current_user=flask_login.current_user) +
            render_template('request.html', form1=form1, form2=form2, message=False) +
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
    ssoform = SSOForm()
    message = None
    urls = []
    if form2.submit_obj.data and form2.validate_on_submit():
        if request.method == 'POST' and form2.validate():
            if form2.RA.data and form2.DEC.data:
                found = db.get_objects_near(form2.RA.data, form2.DEC.data, form2.radius.data)
            elif form2.object_name.data:
                found = db.get_object_id_from_name(form2.object_name.data)
            else:
                found = None

            if found is None:
                message = "Not enough information provided"
            elif form2.RA.data and form2.DEC.data and not found:
                message = "No objects in database within the coordinate range entered"
            elif form2.object_name.data and not found:
                message = "No object names in database contain %s" % (form2.object_name.data,)
            elif found[0] == -1:
                message = found[1]
            elif len(found) > 1:
                urls = [(obj[0], obj[1]) for obj in found]
                message = "multiple objects matching request:"
            else:
                return redirect(url_for('show_objects', ident=found[0][0]))
    # TODO: make an objects "homepage" to have as the base render
    form2.radius.data = 4
    return (render_template('header.html') +#, current_user=flask_login.current_user) +
            render_template('object.html', form1=form1, form2=form2, ssoform=ssoform, object=None,
                            message=message, urls = urls) +
            render_template('footer.html'))


@app.route('/objects/<path:ident>')
def show_objects(ident):
    # take name or id and show info about it and images with permission
    try:
        iden = int(ident)
    except ValueError:
        name = str(ident)
        iden = db.get_from_object(['id'], {'name': name})
        if not iden:
            flash("Object with name '%s' not found, please search for object" % (ident,))
            return redirect(url_for('objects'))
        elif iden[0] == -1:
            flash("Searching for name '%s' caused error '%', please report issue" % (ident, iden[1]))
            return redirect(url_for('objects'))
        else:
            return redirect(url_for('show_objects', iden[0][0]))

    info = db.get_from_object(['name', 'ra', 'dec', 'typedesig', 'epoch', 'id'], {'id': iden})[0]
    # TODO: work in SSO objects
    if info[0] == -1:
        flash("Searching for id '%s' caused error '%s', please report issue")
        return redirect(url_for('objects'))
    elif not info:
        flash("Searching for id '%s' produced no results, please search for object")
        return redirect(url_for('objects'))


    # retrieve all of the requests submitted for the object
    request_query = ("SELECT u.name, u.username, g.designator, r.inidate, r.enddate, r.priority, r.status "
                     "FROM request r, users u, object o, groups g  WHERE r.user_id = u.id AND r.object_id = o.id AND g.id = r.program_id "
                     "AND r.object_id = '%s';" % (ident,))
                    # TODO: add ```AND r.program_id IN %s;" % (,flask_login.current_user.groups)```
    req = db.execute_sql(request_query)
    for n, x in enumerate(req):
        req[n] = list(x)
        # decide between name/username
        if not req[n][0]:
            del req[n][0]
        else:
            del req[n][1]

    # get the time, duration and file for observations of the object
    a_query = ("SELECT obs.mjd, obs.exptime, obs.filter, obs.fitsfile "
               "FROM observation obs, object obj, request r WHERE obj.id = obs.object_id AND "
               "r.id = obs.request_id AND obj.id = '%s';" % (ident,))
               # TODO: add ```AND r.program_id IN %s;" % (, flask_login.current_user.groups)```
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
    # TODO: check again
    # TODO: check even more
    """
       day = datetime.datetime.now()
       start = Time(day - datetime.timedelta(days=day.weekday())).mjd
       end = start + 7
       # TODO: make a form to set start/end dates
       areq_query = ("SELECT atomicrequest.exptime, atomicrequest.status FROM atomicrequest "
                     "JOIN request ON request.id = atomicrequest.request_id "
                     "WHERE request.program_id = '%s'"
                     # " AND  atomicrequest.lastmodified < %s AND atomicrequest.lastmodified > %s"
                     ";" % (group_dict[project],))
       requests = np.array(db.execute_sql(areq_query))
       # TODO: test for whether there are any matching requests
       pending = (requests.T[0] == 'PENDING')
       pending_time = sum(requests.T[1][pending]) / 3600.
       observed = ((requests.T[0] == 'OBSERVED') | (requests.T[0] == 'REDUCED'))
       observed_time = sum(requests.T[1][observed]) / 3600.
       # TODO: make time_allocated part of group definition?
       time_allocated = 14.
       # pygal graphing

       # use different color schemes based on whether there is more time requested than allowed
       if observed_time + pending_time > time_allocated:
           c_style = Style(
               background='transparent',
               colors=('#52e852', '#1a4fba', '#e22626')
           )
           pi_chart = pygal.Pie(style=c_style, height=300, width=300)
           pi_chart.title = "Time allocation (h)"
           pi_chart.add('observed', observed_time)
           pi_chart.add('requested', (time_allocated - observed_time))
           pi_chart.add('requested beyond allocation', (observed_time + pending_time - time_allocated))
       else:
           c_style = Style(
               background='transparent',
               colors=('#52e852', '#1a4fba', '#e5f2de')
           )
           pi_chart = pygal.Pie(style=c_style, height=300, width=300)
           pi_chart.title = "Time allocation (h)"
           pi_chart.add('observed', observed_time)
           pi_chart.add('requested', pending_time)
           pi_chart.add('free', (time_allocated - pending_time - observed_time))
       """

    return (render_template('header.html') +  # , current_user=flask_login.current_user) +
            render_template('object_stats.html', info=info, observations=obs, req_data=req) +
            render_template('footer.html'))


@app.route('/project_stats')
def projects_home():  # TODO: list groups and projects. require login?
    # generate a list of (group designator, project id, project designator)
    projects = []    
    for group in group_programs_dict.values():
        for project in group[1]:
            projects.append((group[0], project[0], project[1]))
    return (render_template('header.html') + 
            render_template('projects_home.html', projects = projects) + 
            render_template('footer.html'))


@app.route('/project_stats/<program>')
@flask_login.login_required
def project_stats(program):
    day = datetime.datetime.now()
    start = Time(day - datetime.timedelta(days=day.weekday())).mjd
    end = start + 7
    prog_id = None
    for lis in group_programs_dict.itervalues():
        for prog in lis[1]:
            if prog[1] == program:
                prog_id = prog[0]
    if not prog_id:
        flash('program %s does not exist' % (program,))
        return redirect(url_for('index'))

    # TODO: make a form to set start/end dates
    time_used_query = ("SELECT o.time_elapsed, r.status FROM request r, observation o "
                       "WHERE o.request_id = r.id AND r.program_id = '%s'"
                       #" AND r.lastmodified < %s AND r.lastmodified > %s"
                       ";" % (prog_id,))
    requests = db.execute_sql(time_used_query)
    time_allocated = db.get_from_program(['time_allocated'], {'id': prog_id})[0][0]
    
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
                     "AND r.program_id = '%s';" % (prog_id,))
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

    return (render_template('header.html') +#, current_user=flask_login.current_user) +
            render_template('project_stats.html', img_data=pi_chart.render_data_uri(), req_data=req) +
            render_template('footer.html'))


@login_manager.unauthorized_handler
def unauthorized_handler():
    return flask.redirect(flask.url_for('login'))

if __name__ == '__main__':
    app.run()


