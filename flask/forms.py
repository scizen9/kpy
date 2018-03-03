from flask_wtf import FlaskForm
from wtforms import fields, validators, widgets, Field, FormField
from flask import request, redirect, url_for
from urlparse import urlparse, urljoin

SECRET_KEY = 'secret'


class HeliocentricEllipticalForm(FlaskForm):
    class Meta:
        csrf=False
    i = fields.FloatField('inclination (degrees)', validators=[validators.input_required()])
    O = fields.FloatField('longitude of ascending node (degrees)', validators=[validators.input_required()])
    o = fields.FloatField('argument of perihelion (degrees)', validators=[validators.input_required()])
    a = fields.FloatField('mean distance (AU)', validators=[validators.input_required()])
    n = fields.FloatField('mean daily motion (degrees/day)', validators=[validators.input_required()])
    e = fields.FloatField('eccentricity', validators=[validators.input_required()])
    M = fields.FloatField('mean anomaly', validators=[validators.input_required()])
    E = fields.DateField('epoch date (time of mean anomaly) (Y-m-d)', validators=[validators.input_required()], format='%Y-%m-%d')
    D = fields.FloatField('equinox year', validators=[validators.input_required()])
    M1 = fields.FloatField('magnitude model first component', validators=[validators.input_required()])
    M2 = fields.FloatField('magnitude model second component', validators=[validators.input_required()])


class HeliocentricHyperbolicForm(FlaskForm):
    class Meta:
        csrf=False
    T = fields.DateField('date of ephoch of perihelion (Y-m-d)', validators=[validators.input_required()], format='%Y-%m-%d')
    i = fields.FloatField('inclination (degrees)', validators=[validators.input_required()])
    O = fields.FloatField('longitude of ascending node (degrees)', validators=[validators.input_required()])
    o = fields.FloatField('argument of perihelion (degrees)', validators=[validators.input_required()])
    e = fields.FloatField('eccentricity', validators=[validators.input_required()])
    q = fields.FloatField('perihelion distance (AU)', validators=[validators.input_required()])
    D = fields.FloatField('equinox year', validators=[validators.input_required()])
    M1 = fields.FloatField('magnitude model first component', validators=[validators.input_required()])
    M2 = fields.FloatField('magnitude model second component', validators=[validators.input_required()])


class HeliocentricParabolicForm(FlaskForm):
    class Meta:
        csrf=False
    T = fields.DateField('date of ephoch of perihelion (Y-m-d)', validators=[validators.input_required()], format='%Y-%m-%d')
    i = fields.FloatField('inclination (degrees)', validators=[validators.input_required()])
    o = fields.FloatField('argument of perihelion (degrees)', validators=[validators.input_required()])
    q = fields.FloatField('perihelion distance (AU)', validators=[validators.input_required()])
    O = fields.FloatField('longitude of ascending node (degrees)', validators=[validators.input_required()])
    D = fields.FloatField('equinox year', validators=[validators.input_required()])
    M1 = fields.FloatField('magnitude model first component', validators=[validators.input_required()])
    M2 = fields.FloatField('magnitude model second component', validators=[validators.input_required()])


class EarthSatelliteForm(FlaskForm):
    class Meta:
        csrf=False
    T = fields.DateField('Epoch of following fields (Y-m-d)', validators=[validators.input_required()], format='%Y-%m-%d')
    i = fields.FloatField('inclination (degrees)', validators=[validators.input_required()])
    ra = fields.FloatField('RA of ascending node (degrees)', validators=[validators.input_required()])
    e = fields.FloatField('Eccentricity', validators=[validators.input_required()])
    P = fields.FloatField('argument of perigee (degrees)', validators=[validators.input_required()])
    M = fields.FloatField('mean anomaly (degrees)', validators=[validators.input_required()])
    n = fields.FloatField('mean motion (revs/day)', validators=[validators.input_required()])
    d = fields.FloatField('orbit decay rate (revs/day^2)', validators=[validators.input_required()])
    r = fields.FloatField('integral reference orbit no. at Epoch', validators=[validators.input_required()])
    #dr = fields.FloatField('drag coefficient (1/earth radii)', validators=[validators.Optional()])


class PeriodicForm(FlaskForm):
    class Meta:
        csrf=False
    mdj0 = fields.IntegerField('mjd of period start', validators=[validators.input_required()])
    phasedays = fields.FloatField('period length (days)', validators=[validators.input_required()])


class RequestForm(FlaskForm):
    # object data
    obj_name = fields.StringField('Object Name', [validators.input_required()])
    obj_ra = fields.StringField('RA', [validators.Optional()])
    obj_dec = fields.StringField('DEC', [validators.Optional()])
    typedesig = fields.SelectField('object type', [validators.input_required()], coerce=str, 
                                   choices=[('f', 'fixed'), ('v', 'periodic fixed'), ('e', 'heliocentric elliptical'),
                                            ('h', 'heliocentric hyperbolic'), ('p', 'heliocentric parabolic'), ('E', 'geocentric elliptical')])
    # required and optional fields for requests
    allocation = fields.SelectField('allocation', [validators.data_required()], coerce=int, choices=[])
    priority = fields.FloatField('priority', [validators.input_required()])
    filters_op = fields.SelectField('filters', [validators.input_required()], coerce=str, 
                                    choices=[(' 1, 1, 1, 1}', 'u-g-r-i'), (' 0, 1, 1, 1}', 'g-r-i'), \
                                            (' 0, 1, 0, 0}', 'g'), (' 0, 0, 1, 0}', 'r'), ('0, 0, 0, 1}', 'i'), ('0, 0, 0, 0}', '-') ])
    seq_repeats = fields.IntegerField('# of repeats', default = 1)
    ifu = fields.BooleanField('IFU image', [validators.Optional()])
    ab = fields.BooleanField('A B pair', [validators.Optional()])
    cadence = fields.FloatField('cadence', default=None, validators=(validators.Optional(),))
    min_moon_dist = fields.FloatField('minimum moon distance (degrees)',
                                      validators=(validators.Optional(), validators.number_range(0., 180.)))
    max_moon_illum = fields.FloatField('maximum moon illumination (fractional)',
                                       validators=(validators.Optional(), validators.number_range(0., 1.)))
    max_cloud_cover = fields.FloatField('maximum cloud cover (fractional)',
                                        validators=(validators.Optional(), validators.number_range(0., 1.)))
    phasesamples = fields.FloatField('samples per period', default=None, validators=(validators.Optional(),))
    inidate = fields.DateField('start date (Y-m-d)', validators=[validators.input_required()], format='%Y-%m-%d')
    enddate = fields.DateField('end date (Y-m-d)', validators=[validators.input_required()], format='%Y-%m-%d')
    #elliptical = fields.FormField(HeliocentricEllipticalForm, validators=(validators.Optional(),))
    #hyperbolic = fields.FormField(HeliocentricHyperbolicForm, validators=(validators.Optional(),))
    #parabolic = fields.FormField(HeliocentricParabolicForm, validators=(validators.Optional(),))
    #satellite = fields.FormField(EarthSatelliteForm, validators=(validators.Optional(),))
    submit_req = fields.SubmitField('submit request')


class FindObjectForm(FlaskForm):
    obj_name = fields.StringField('object name contains', validators=(validators.Optional(),))
    typedesig = fields.SelectField('object type', [validators.input_required()], coerce=str, 
                                   choices=[('',''), ('f', 'fixed'), ('v', 'periodic fixed'), ('e', 'heliocentric elliptical'),
                                            ('h', 'heliocentric hyperbolic'), ('p', 'heliocentric parabolic'), ('E', 'geocentric elliptical')])
    obj_ra = fields.StringField('Right Ascension (deg or HH:MM:SS)', validators=(validators.Optional(),))
    obj_dec = fields.StringField('Declination (deg or DD:MM:SS)', validators=(validators.Optional(),))
    radius = fields.FloatField('radius (arcsec)', validators=(validators.Optional(),))
    submit_obj = fields.SubmitField('seach for object')


class SubmitObjectForm(FlaskForm):
    obj_name = fields.StringField('object Name')
    typedesig = fields.SelectField('object type', [validators.input_required()], coerce=str, 
                                   choices=[('f', 'fixed'), ('v', 'periodic fixed'), ('e', 'heliocentric elliptical'),
                                            ('h', 'heliocentric hyperbolic'), ('p', 'heliocentric parabolic'), ('E', 'geocentric elliptical')])
    obj_ra = fields.StringField('Right Ascension (deg or HH:MM:SS)', validators=(validators.Optional(),))
    obj_dec = fields.StringField('Declination (deg or DD:MM:SS)', validators=(validators.Optional(),))
    add_obj = fields.SubmitField('Add object')


class LoginForm(FlaskForm):
    username = fields.StringField('username')
    password = fields.PasswordField('password')


class PassChangeForm(FlaskForm):
    password = fields.PasswordField('Old Password', validators=[validators.input_required()])
    pass_new = fields.PasswordField('New Password', validators=[validators.input_required(), validators.EqualTo('pass_conf', message='Passwords must match')] )
    pass_conf = fields.PasswordField('Confirm New Password', validators=[validators.input_required()])


class SearchUserForm(FlaskForm):

    search_string = fields.StringField('Introduce a string to search for your user.')
    search_user = fields.SubmitField('Search User', description=None)

class UsersForm(FlaskForm):

    username = fields.StringField('name', validators=[validators.input_required()])
    name = fields.StringField('name', validators=[validators.input_required()])
    email = fields.StringField('email', validators=[validators.input_required()])
    password = fields.PasswordField('password')
    pass_new = fields.PasswordField('New Password', validators=[validators.EqualTo('pass_conf', message='Passwords must match')] )
    pass_conf = fields.PasswordField('Confirm New Password')
    modify_user = fields.SubmitField('Modify User', description=None)

    old_groups = fields.SelectField('Select group to remove')
    new_groups = fields.SelectField('Select group to add')
    add_group = fields.SubmitField('Add', description=None)
    remove_group = fields.SubmitField('Remove', description=None)

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
