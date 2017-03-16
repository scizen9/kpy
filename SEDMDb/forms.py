from flask_wtf import FlaskForm
from wtforms import BooleanField, StringField, PasswordField, validators, IntegerField, DateField, FloatField, SubmitField
from flask import request

#form = FlaskForm(csrf_enabled=False)
SECRET_KEY = 'secret'


class RequestForm(FlaskForm):
    object_id = IntegerField('object_id', [validators.data_required()])
    priority = FloatField('priority', [validators.data_required()])
    inidate = DateField('start date (Y-M-D)', [validators.data_required()])
    enddate = DateField('end date (Y-M-D)', [validators.data_required()])
    submit_req = SubmitField('submit request')


class FindObjectForm(FlaskForm):
    object_name = StringField('object name')
    RA = FloatField('Right Ascension (deg)')
    DEC = FloatField('Declination (deg)')
    submit_obj = SubmitField('seach for object')
