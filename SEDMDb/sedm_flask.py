"""

"""


import os
import datetime
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

    if form1.submit_req.data and form1.validate_on_submit():
        print 'form1'
        if request.method == 'POST' and form1.validate():
            flash('good request')
            req = db.add_request({'object_id': form1.object_id.data, 'exptime': '{120, 2400}', 'priority': form1.priority.data,
                                  'inidate': form1.inidate.data, 'enddate': form1.enddate.data, 'user_id':0, 'program_id': 0})
            # TODO: set up user_id and program_id
            print req
            if req[0] == 0:
                return redirect(flask.url_for('/'))
            else:
                return render_template('header.html') + render_template('request.html', form=form1, fail=req[1]) \
                       + render_template('find_object.html', form=form2) + render_template('footer.html')

    if form2.submit_obj.data and form2.validate_on_submit():
        print 'form2'
        if request.method == 'POST' and form2.validate():
            ra = form2.RA.data
            dec = form2.DEC.data
            return redirect('/#')
    return render_template('header.html') + render_template('request.html', form=form1, fail=False)\
        + render_template('find_object.html', form=form2) + render_template('footer.html')


@app.route('/schedule')
def schedule():
    pass


@app.route('/objects')
def objects():
    pass

if __name__ == '__main__':
    app.run()


