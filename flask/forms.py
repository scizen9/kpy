from flask_wtf import FlaskForm
from wtforms import fields, validators
from flask import request, redirect, url_for
from urlparse import urlparse, urljoin

SECRET_KEY = 'secret'


class RequestForm(FlaskForm):
    object_id = fields.IntegerField('object_id', [validators.data_required()])
    priority = fields.FloatField('priority', [validators.data_required()])
    inidate = fields.DateField('start date (Y-M-D)', [validators.data_required()])
    enddate = fields.DateField('end date (Y-M-D)', [validators.data_required()])
    user_id = fields.HiddenField('user_id', [validators.data_required()])
    project = fields.SelectField('project', [validators.data_required()], coerce=int, choices=[])
    submit_req = fields.SubmitField('submit request')


class FindObjectForm(FlaskForm):
    object_name = fields.StringField('object name contains')
    RA = fields.FloatField('Right Ascension (deg)')
    DEC = fields.FloatField('Declination (deg)')
    submit_obj = fields.SubmitField('seach for object')


class SubmitObjectForm(FlaskForm):
    object_name = fields.StringField('object Name')
    RA = fields.FloatField('Right Ascension (deg)')
    DEC = fields.FloatField('Declination (deg)')
    add_obj = fields.SubmitField('Add object')


class SSOForm(FlaskForm):
    typedesig = fields.SelectField('Target type', choices=[('fixed', 'f'), ('heliocentric elliptical', 'e'),
                                   ('heliocentric hyperbolic', 'h'), ('heliocentric parabolic', 'p'),
                                   ('geocentric elliptical', 'E')])


class LoginForm(FlaskForm):
    username = fields.StringField('username')
    password = fields.StringField('password')


def is_safe_url(target):
    ref_url = urlparse(request.host_url)
    test_url = urlparse(urljoin(request.host_url, target))
    return test_url.scheme in ('http', 'https') and \
           ref_url.netloc == test_url.netloc


def get_redirect_target():
    for target in request.args.get('next'), request.referrer:
        if not target:
            continue
        if is_safe_url(target):
            return target


class RedirectForm(FlaskForm):
    next = fields.HiddenField()

    def __init__(self, *args, **kwargs):
        FlaskForm.__init__(self, *args, **kwargs)
        if not self.next.data:
            self.next.data = get_redirect_target() or ''

    def redirect(self, endpoint='index', **values):
        if is_safe_url(self.next.data):
            return redirect(self.next.data)
        target = get_redirect_target()
        return redirect(target or url_for(endpoint, **values))