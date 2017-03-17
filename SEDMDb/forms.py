from flask_wtf import FlaskForm
from wtforms import fields, validators
from flask import request

#form = FlaskForm(csrf_enabled=False)
SECRET_KEY = 'secret'


class RequestForm(FlaskForm):
    object_id = fields.IntegerField('object_id', [validators.data_required()])
    priority = fields.FloatField('priority', [validators.data_required()])
    inidate = fields.DateField('start date (Y-M-D)', [validators.data_required()])
    enddate = fields.DateField('end date (Y-M-D)', [validators.data_required()])
    user_id = fields.HiddenField('user_id', [validators.data_required()])
    project = fields.SelectField('project', [validators.data_required()], choices=[])
    submit_req = fields.SubmitField('submit request')


class FindObjectForm(FlaskForm):
    object_name = fields.StringField('object name contains')
    RA = fields.FloatField('Right Ascension (deg)')
    DEC = fields.FloatField('Declination (deg)')
    submit_obj = fields.SubmitField('seach for object')
