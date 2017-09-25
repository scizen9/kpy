from flask_wtf import FlaskForm
from wtforms import fields, validators, widgets, Field, FormField
from flask import request, redirect, url_for
from urlparse import urlparse, urljoin

SECRET_KEY = 'secret'


#class IntListField(Field):
#    widget = widgets.ListWidget()


class FilterForm(FlaskForm):
    ifu_val = fields.IntegerField('ifu', default=0, validators=(validators.Optional(),))
    u_val = fields.IntegerField('u', default=0, validators=(validators.Optional(),))
    g_val = fields.IntegerField('g', default=0, validators=(validators.Optional(),))
    r_val = fields.IntegerField('r', default=0, validators=(validators.Optional(),))
    i_val = fields.IntegerField('i', default=0, validators=(validators.Optional(),))


class RequestForm(FlaskForm):
    object_id = fields.IntegerField('object_id', [validators.data_required()])
    user_id = fields.HiddenField('user_id', [validators.data_required()])
    program = fields.SelectField('program', [validators.data_required()], coerce=int, choices=[])
    priority = fields.FloatField('priority', [validators.data_required()])
    filters = fields.FieldList(FormField(FilterForm), min_entries=1, label='Obs per filter')
    ab = fields.BooleanField('A B pair', [validators.data_required()])
    cadence = fields.FloatField('cadence', default=None, validators=(validators.Optional(),))
    phasesamples = fields.FloatField('samples per period', default=None, validators=(validators.Optional(),))
    inidate = fields.DateField('start date (Y-M-D)', validators=[validators.data_required()], format='%Y-%m-%d')
    enddate = fields.DateField('end date (Y-M-D)', validators=[validators.data_required()], format='%Y-%m-%d')
    submit_req = fields.SubmitField('submit request')


class FindObjectForm(FlaskForm):
    object_name = fields.StringField('object name contains', validators=(validators.Optional(),))
    RA = fields.FloatField('Right Ascension (deg)', validators=(validators.Optional(),))
    DEC = fields.FloatField('Declination (deg)', validators=(validators.Optional(),))
    radius = fields.FloatField('radius (arcsec)', validators=(validators.Optional(),))
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
    password = fields.PasswordField('password')


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
