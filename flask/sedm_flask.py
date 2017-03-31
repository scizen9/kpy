"""

"""
import sys
import os
sys.path.append(os.path.abspath('../'))

import os
import datetime
from astropy.time import Time
import numpy as np
from flask import Flask, request, flash, redirect, render_template, url_for
from SEDMDb.SedmDb import SedmDB
from SEDMDb.SedmDb_tools import DbTools
from werkzeug.security import generate_password_hash, check_password_hash
from forms import RequestForm, RedirectForm, FindObjectForm, LoginForm, SubmitObjectForm, SSOForm, is_safe_url
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
    return User.get(user_id)


@app.route('/')
def index():
    # TODO: make it pretty
    return render_template('header.html') + render_template('footer.html')


@app.route('/login', methods=['GET', 'POST'])
def login():
    # Here we use a class of some kind to represent and validate our
    # client-side form data. For example, WTForms is a library that will
    # handle this for us, and we use a custom LoginForm to validate.
    form = LoginForm()
    if form.validate_on_submit():
        # Login and validate the user.
        # user should be an instance of your `User` class
        username = form.username.data
        password = form.password.data
        user_pass = db.get_from_users(['username', 'password'], {'username': username})
        if not user_pass:
            message = "Incorrect username or password"
            return render_template('header.html', current_user=flask_login.current_user) + \
                   render_template('login.html', form=form, message=message) + \
                   render_template('footer.html')

        elif user_pass[0] == -1:
            message = user_pass[1]
            return render_template('header.html', current_user=flask_login.current_user) + \
                   render_template('login.html', form=form, message=message) + \
                   render_template('footer.html')
        elif check_password_hash(user_pass[0][1], flask.request.form['password']):
            user = User()
            user.id = username
            flask_login.login_user(user, remember=True)
            return redirect(flask.url_for('index'))
    return render_template('header.html', current_user=flask_login.current_user) + \
                   render_template('login.html', form=form, message=None) + \
                   render_template('footer.html')


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
            ini = Time(str(form1.inidate.data))
            end = Time(str(form1.enddate.data))
            if ini < end:
                req = (0, 'object_id: %s, priority: %s, inidate: %s, enddate: %s, exptime: {120,2400}' %
                       (form1.object_id.data, form1.priority.data, form1.inidate.data, form1.enddate.data))
            else:
                req = (-1, 'ERROR: inidate: %s after enddate: %s' % (form1.inidate.data, form1.enddate.data))
            # TODO: return to following when able to access database
            #req = db.add_request({'object_id': form1.object_id.data, 'exptime': '{120, 2400}', 'priority': form1.priority.data,
            #                      'inidate': form1.inidate.data, 'enddate': form1.enddate.data, 'user_id':0, 'program_id': 0})
            # TODO: set up user_id and program_id=
            message = req[1]
            return render_template('header.html', current_user=flask_login.current_user) +\
                   render_template('request.html', form1=form1, form2=form2, message=message) + \
                   render_template('footer.html')

    if form2.submit_obj.data and form2.validate_on_submit():
        print 'object search'
        if request.method == 'POST' and form2.validate():
            if form2.RA.data and form2.DEC.data:
                ra = form2.RA.data
                dec = form2.DEC.data
                rad = form2.radius.data
                # req = db.get_objects_near(ra, dec, rad)
                # TODO: uncomment db command
                req = [(1, 'placeholder')]
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
            return render_template('header.html', current_user=flask_login.current_user) + \
                   render_template('request.html', form1=form1, form2=form2, message=message) + \
                   render_template('footer.html')
    return render_template('header.html', current_user=flask_login.current_user) + \
           render_template('request.html', form1=form1, form2=form2, message=False) + \
           render_template('footer.html')
# TODO: set up ability to view and cancel requests
# TODO: allow searching for objects, return link to page with information/observations of it
# TODO: End of Night report as html


@app.route('/schedule')
def schedule():
    # TODO: create a table
    schedule = None  # TODO: decide how schedule will be stored (file/database?)
    return render_template('header.html', current_user=flask_login.current_user) + \
           render_template('schedule.html') + \
           render_template('footer.html')


@app.route('/objects', defaults={'where': ''})
@app.route('/objects/<path:where>')
def objects(where):
    form1 = SubmitObjectForm()
    form2 = FindObjectForm()
    ssoform = SSOForm()
    if where == 'add':
        return render_template('header.html', current_user=flask_login.current_user) + \
               render_template('object.html', form1=form1, form2=form2, ssoform=ssoform, object=None) + \
               render_template('footer.html')
    elif where:
        return render_template('header.html', current_user=flask_login.current_user) + \
               render_template('object.html', object=where) + \
               render_template('footer.html')
    # TODO: make an objects "homepage" to have as the base render
    return render_template('header.html', current_user=flask_login.current_user) + \
           render_template('object.html', form1=form1, form2=form2, ssoform=ssoform, object=None) + \
           render_template('footer.html')


@app.route('/project_stats/<path:project>')
def project_stats(project):
    # TODO: show chart with 'total time allocated', 'time requested', "time observed"
    return render_template('header.html', current_user=flask_login.current_user) + render_template('') + render_template('footer.html')


if __name__ == '__main__':
    app.run()


