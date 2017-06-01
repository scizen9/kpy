"""

"""
import sys
import os
sys.path.append(os.path.abspath('../'))

import os
import datetime
from astropy.time import Time
import numpy as np
from flask import Flask, request, flash, redirect, render_template, url_for, make_response
from SEDMDb.SedmDb import SedmDB
from SEDMDb.SedmDb_tools import DbTools
from werkzeug.security import generate_password_hash, check_password_hash
from forms import RequestForm, RedirectForm, FindObjectForm, LoginForm, SubmitObjectForm, SSOForm, is_safe_url
# from flask.ext.appbuilder.charts.views import DirectByChartView
import flask
import flask_login
from wtforms import Form, HiddenField


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
group_dict = {}
for group in groups:
    group_dict[group[1]] = group[0]


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
    user.id = users[0][1]
    user.name = users[0][1]
    return user


@app.route('/')
def index():
    # TODO: make it pretty
    return render_template('header.html', current_user=flask_login.current_user) + render_template('footer.html')


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
            return (render_template('header.html') +#, current_user=flask_login.current_user) +
                    render_template('login.html', form=form, message=message) +
                    render_template('footer.html'))

        elif user_pass[0] == -1:
            message = user_pass[1]
            return (render_template('header.html') +#, current_user=flask_login.current_user) +
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
            return (render_template('header.html') +#, current_user=flask_login.current_user) +
                    render_template('login.html', form=form, message=message) +
                    render_template('footer.html'))
    return (render_template('header.html') +#, current_user=flask_login.current_user) +
            render_template('login.html', form=form, message=None) +
            render_template('footer.html'))


@app.route('/request', methods=['GET', 'POST'])
def requests():
    """
    generate forms for request creation
    """
    form1 = RequestForm()
    form2 = FindObjectForm()
    # if form1.validate_on_submit():
    #    print 'form1'
    form1.project.choices = [(1, 'ZTF'), (2, 'CIT')]
    form1.user_id.data = 1
    if form1.submit_req.data and form1.validate_on_submit():
        print "request submitted"
        if request.method == 'POST':
            req = db.add_request({'object_id': form1.object_id.data, 'exptime': '{120, 2400}', 'priority': form1.priority.data,
                                  'inidate': form1.inidate.data, 'enddate': form1.enddate.data, 'user_id':0, 'program_id': form1.project.data})
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
                message = "multiple objects matching coordinates:"
            else:
                return redirect(url_for('show_objects', ident=found[0][0]))
    # TODO: make an objects "homepage" to have as the base render
    form2.radius.data = 4
    return (render_template('header.html') +#, current_user=flask_login.current_user) +
            render_template('object.html', form1=form1, form2=form2, ssoform=ssoform, object=None,
                            message=message, urls = urls) +
            render_template('footer.html'))


# TODO: make show_objects prettier, add parts of it to 'objects' after search completes
@app.route('/objects/<path:ident>')
def show_objects(ident):
    # take name or id and show info about it and images with permission
    try:
        iden = int(ident)
    except ValueError:
        name = str(ident)
        iden = None
    if iden:
        info = db.get_from_object(['name', 'RA', 'DEC', 'typedesig', 'epoch', 'id'], {'id': iden})[0]
    elif name:
        info = db.get_from_object(['name', 'RA', 'DEC', 'typedesig', 'epoch', 'id'], {'name': name})[0]

    if info[0] == -1:
        # TODO: show a query failed with message ""
        # send report?
        pass
    elif not info:
        # TODO: redirect or show "there is no object %s" % (ident,)
        pass
    elif not flask_login.current_user:
        return (render_template('header.html',  current_user=flask_login.current_user) +
                render_template('show_object.html'))
    else:
        query = ("SELECT atomicrequest.id, atomicrequest.exptime, atomicrequest.status, request.program_id FROM request"
                 " JOIN object ON object.id = request.object_id "
                 "JOIN atomicrequest ON request.id = atomicrequest.request_id"
                 "WHERE object.id = '%s';" % (info[5],))
        obs = db.execute_sql(query)
        # TODO: remove above query?
        # TODO: setup sso object info (below)
        """
        if info[3] == 'e':
            sso_query =
        elif info[3] == 'h':

        elif info[3] == 'p':

        elif info[3] == 'E':
        """
        for areq in obs:
            pass
        # get the time, duration and file for observations of the object
        areq_query = ("SELECT observation.mjd, observation.exptime, observation.fitsfile FROM observation"
                      "JOIN object ON object.id = observation.object_id"
                      "JOIN request ON request.id = observation.request_id"
                      "WHERE object.id = '%s' AND request.program_id IN %s;" % (info[5], flask_login.current_user.groups))
        observations = db.execute_sql(areq_query)




    return (render_template('header.html') +#, current_user=flask_login.current_user) +
            render_template('show_object.html', object=ident) +
            render_template('footer.html'))


@app.route('/project_stats/<path:project>')
def project_stats(project):
    # TODO: show chart with 'total time allocated', 'time requested', "time observed"

    day = datetime.datetime.now()
    start = Time(day - datetime.timedelta(days=day.weekday())).mjd
    end = Time(start + datetime.timedelta(days=7)).mjd

    areq_query = ("SELECT atomicrequest.exptime, atomicrequest.status FROM atomicrequest"
                  "JOIN request ON request.id = atomicrequest.request_id"
                  "WHERE request.project_id = '%s'"
                  #" AND  atomicrequest.lastmodified < %s AND atomicrequest.lastmodified > %s"
                  ";" % (group_dict[project],))
    requests = np.array(db.execute_sql(areq_query))
    # TODO: test for whether there are any matching requests
    pending = (requests.T[0] == 'PENDING')
    pending_time = sum(requests.T[1][pending])/3600.
    observed = ((requests.T[0] == 'OBSERVED') | (requests.T[0] == 'REDUCED'))
    observed_time = sum(requests.T[1][observed])/3600.
    # TODO: make time_allocated part of group definition?
    time_allocated = 14.
    # TODO: use better method than below (pygal or flask-appbuilder?)
    import StringIO
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    from matplotlib.figure import Figure

    fig = Figure()
    ax = fig.add_subplot(111)
    labels = ['observed', 'pending', 'over_budget', 'free']
    if observed_time+pending_time > time_allocated:
        total_time = observed_time+pending_time
        fracs = [observed_time/total_time, (time_allocated-observed_time)/total_time, (total_time-time_allocated)/total_time, 0.]
    else:
        fracs = [observed_time/time_allocated, pending_time/time_allocated, 0., (time_allocated-pending_time-observed_time)]
    ax.pie(fracs, labels=labels, autopct='%1.1f%%', shadow=True)
    canvas = FigureCanvas(fig)
    png_output = StringIO.StringIO()
    canvas.print_png(png_output)
    response = make_response(png_output.getvalue())
    response.headers['Content-Type'] = 'image/png'


# return response
    return (render_template('header.html', current_user=flask_login.current_user) +
            render_template('project_stats', img_data=response) +
            render_template('footer.html'))


if __name__ == '__main__':
    app.run()


