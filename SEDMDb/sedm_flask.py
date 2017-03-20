"""

"""


import os
import datetime
from astropy.time import Time
import numpy as np
from flask import Flask, request, flash, redirect, render_template
from SedmDb import SedmDB
from SedmDb_tools import DbTools
import flask


# config
SECRET_KEY = 'secret'
USERNAME = 'admin'
PASSWORD = 'default'

db = SedmDB()
tools = DbTools(db)

app = Flask(__name__)
app.config.from_object(__name__)


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


@app.route('/')
def index():
    return render_template('header.html') + render_template('footer.html')


import forms


@app.route('/request', methods=['GET', 'POST'])
def requests():
    # TODO: make this capable of handling multiple forms (find object by name/ra+dec)
    form1 = forms.RequestForm()
    form2 = forms.FindObjectForm()
    # if form1.validate_on_submit():
    #    print 'form1'
    form1.project.choices = [(1, 'ZTF'), (2, 'CIT')]
    form1.user_id.data = 1
    if form1.submit_req.data and form1.validate_on_submit():
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
            return render_template('header.html') + render_template('request.html', form1=form1, form2=form2,
                                                                    message=message) + render_template('footer.html')

    if form2.submit_obj.data and form2.validate_on_submit():
        print 'form2'
        if request.method == 'POST' and form2.validate():
            if form2.RA.data and form2.DEC.data:
                ra = form2.RA.data
                dec = form2.DEC.data
                # req = db.get_objects_near(ra, dec, 5)
                # TODO: uncomment db command
                req = [(1, 'placeholder')]
                if req:
                    if req[0] == -1:
                        message = req[1]
                        obj_id = None
                    else:
                        message = 'object(s) within 5 arcsec of coordinates: '
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
            return render_template('header.html') + render_template('request.html', form1=form1, form2=form2,
                                                                    message=message) + render_template('footer.html')
    return render_template('header.html') + render_template('request.html', form1=form1, form2=form2,
                                                            message=False) + render_template('footer.html')


@app.route('/schedule')
def schedule():
    pass


@app.route('/objects')
def objects():
    pass

if __name__ == '__main__':
    app.run()


