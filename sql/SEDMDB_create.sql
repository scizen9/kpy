-- Created by Vertabelo (http://vertabelo.com)
-- Last modification date: 2016-08-19 00:28:45.832

-- tables
-- Table: classification
CREATE TABLE classification (
    id BIGSERIAL,
    object_id bigint NOT NULL,
    spec_id bigint NOT NULL,
    classification text  NULL,
    redshift decimal(7,5)  NULL,
    redshift_err decimal(7,5)  NULL,
    phase decimal(5,2)  NULL,
    phase_err decimal(5,2)  NULL,
    classifier text  NULL,
    score decimal(5,2)  NULL,
    CONSTRAINT classification_pk PRIMARY KEY (id)
);

CREATE INDEX classification_object_id_key ON classification(
    object_id
);

CREATE INDEX classification_spec_id_key ON classification(
    spec_id
);

CREATE INDEX classification_classification_key ON classification(
    classification
);

CREATE INDEX classification_redshift_key ON classification(
    redshift
);

-- Table: flexure
CREATE TABLE flexure (
    id BIGSERIAL,
    rms decimal(8,4)  NOT NULL,
    spec_id_1 bigint  NOT NULL,
    spec_id_2 bigint  NOT NULL,
    timestamp1 TIMESTAMP NOT NULL,
    timestamp2 TIMESTAMP NOT NULL,
    CONSTRAINT flexure_pk PRIMARY KEY (id)
);

CREATE INDEX flexure_spec_id_1_key ON flexure(
    spec_id_1
);

-- Table: metrics_phot
CREATE TABLE metrics_phot (
    id BIGSERIAL,
    phot_id bigint NOT NULL UNIQUE,
    fwhm decimal(5,2)  NULL,
    background decimal(5,2)  NULL,
    zp decimal(5,2)  NULL,
    zperr decimal(5,2)  NULL,
    ellipticity decimal(5,2)  NULL,
    nsources int  NULL,
    CONSTRAINT metrics_phot_pk PRIMARY KEY (id)
);

-- Table: metrics_spec
CREATE TABLE metrics_spec (
    id BIGSERIAL,
    spec_id bigint NOT NULL UNIQUE,
    fwhm decimal(5,2)  NULL,
    background decimal(5,2)  NULL,
    line_fwhm int  NULL,
    CONSTRAINT metrics_spec_pk PRIMARY KEY (id)
);

-- Table: object
-- f fixed (or at most exhibits constant curvilinear proper motion)
-- e heliocentric elliptical orbit
-- h heliocentric hyperbolic orbit
-- p heliocentric parabolic orbit
-- E geocentric elliptical orbit, i.e., Earth satellite
-- P built-in planet or natural satellite name
--
CREATE TABLE object (
    id BIGSERIAL,
    marshal_id bigint NULL UNIQUE,
    name text  NOT NULL,
    iauname text NULL UNIQUE,
    ra decimal(12,6) NULL,
    dec decimal(12,6) NULL,
    typedesig varchar(1),
    epoch float,
    creationdate timestamp NOT NULL DEFAULT NOW(),
    CONSTRAINT object_pk PRIMARY KEY (id)
);

CREATE INDEX coords_q3c ON object(
    q3c_ang2ipix(ra, dec)
);

CREATE INDEX object_name_key ON object(
    name
);

--parameters from http://www.clearskyinstitute.com/xephem/help/xephem.html#mozTocId468501
CREATE TABLE elliptical_heliocentric (
    id BIGSERIAL,
    object_id bigint not null UNIQUE,
    --inclination
    inclination decimal(10,8),
    -- longitude of ascending node
    longascnode_O decimal(10,8),
    -- argument of perihelion
    perihelion_o decimal(10,8),
    -- mean distance AU
    a decimal(10,8),
    -- mean daily motion: deg/day
    n decimal(10,8),
    -- eccentricity (<1)
    e decimal(10,8),
    -- Mean anomaly
    M decimal(10,8),
    --epoch date. Time of M
    mjdepoch int,
    -- equinox year
    D int,
    -- first abd second components of magnitude model
    M1 decimal(5,4),
    M2 decimal(5,4),
    -- angula size at 1 AU
    s decimal(10,8) NULL,
    CONSTRAINT sso_pk PRIMARY KEY (id)
);

CREATE TABLE hyperbolic_heliocentric (
    id BIGSERIAL,
    object_id bigint not null UNIQUE,
    -- date of the epoch of perihelion
    T timestamp,
    --inclination
    inclination decimal(10,8),
    -- longitude of ascending node
    longascnode_O decimal(10,8),
    -- argument of perihelion
    perihelion_o decimal(10,8),
    -- eccentricity (<1)
    e decimal(10,8),
    -- perihelion distance, AU
    q decimal(10,8),
    -- equinox year
    D int,
    -- first and second components of magnitude model
    M1 decimal(5,4),
    M2 decimal(5,4),
    -- angular size at 1 AU
    s decimal(10,8) NULL,
    CONSTRAINT hyperbolic_heliocentric_pk PRIMARY KEY (id)
);

CREATE TABLE parabolic_heliocentric (
    id BIGSERIAL,
    object_id bigint not null UNIQUE,
    -- date of the epoch of perihelion
    T timestamp,
    --inclination
    inclination decimal(10,8),
    -- argument of perihelion
    perihelion_o decimal(10,8),
    -- perihelion distance, AU
    q decimal(10,8),
    -- longitude of ascending node
    longascnode_O decimal(10,8),
    -- equinox year
    D int,
    -- first and second components of magnitude model
    M1 decimal(5,4),
    M2 decimal(5,4),
    -- angular size at 1 AU
    s decimal(10,8) NULL,
    CONSTRAINT parabolic_heliocentric_pk PRIMARY KEY (id)
);

CREATE TABLE earth_satellite (
    id BIGSERIAL,
    object_id bigint not null UNIQUE,
    -- first date the elements are valid
    T timestamp,
    --inclination
    inclination decimal(10,8),
    -- RA of ascending node
    ra decimal(10,8),
    -- eccentricity (<1)
    e decimal(10,8),
    -- argument of pedigree
    pedigree decimal(10,8),
    -- Mean anomaly
    M decimal(10,8),
    -- mean motion, revs/day
    n decimal(10,8),
    -- orbit decay rate, rev/day^2
    decay decimal(10,8),
    -- integral reference orbit number at epoch
    reforbit int,
    -- drag coefficient, 1/(earth radii)
    drag decimal(10,8),
    CONSTRAINT earth_satellite_pk PRIMARY KEY (id)
);

CREATE TABLE periodic (
    id BIGSERIAL,
    object_id bigint not null,
    mjd0 decimal (10,8),
    phasedays decimal(10,8),
    CONSTRAINT periodic_pk PRIMARY KEY (id)
);

CREATE INDEX periodic_object_id_key ON periodic(
    object_id
);

-- Table: observation
CREATE TABLE observation (
    id BIGSERIAL,
    object_id bigint NOT NULL,
    request_id bigint NOT NULL,
    atomicrequest_id bigint NOT NULL UNIQUE,
    mjd decimal(10,5)  NOT NULL,
    airmass decimal(5,2)  NOT NULL,
    exptime decimal(6,2)  NOT NULL,
    fitsfile text  NOT NULL UNIQUE,
    imtype text  NULL,
    lst text  NOT NULL,
    ra decimal(12,6)  NOT NULL,
    dec decimal(12,6)  NOT NULL,
    tel_ra text  NOT NULL,
    tel_dec text  NOT NULL,
    tel_az decimal(5,2)  NOT NULL,
    tel_el decimal(5,2)  NOT NULL,
    tel_pa decimal(5,2)  NOT NULL,
    ra_off decimal(5,2)  NOT NULL,
    dec_off decimal(5,2)  NOT NULL,
    camera text  NULL,
    CONSTRAINT observation_pk PRIMARY KEY (id)
);

CREATE INDEX observation_request_id_key ON observation(
    request_id
);

CREATE INDEX observation_object_id_key ON observation(
    object_id
);

-- Table: spec
CREATE TABLE spec (
    id BIGSERIAL,
    spec_calib_id bigint NOT NULL,
    observation_id bigint NOT NULL UNIQUE,
    asciifile text  NULL,
    npyfile text NULL,
    fitsfile text NULL,
    imgset text  NULL,
    quality int  NOT NULL,
    cubefile text  NULL,
    standardfile text  NULL,
    marshal_spec_id bigint  NULL,
    skysub boolean  NOT NULL,
    CONSTRAINT spec_pk PRIMARY KEY (id)
);

CREATE INDEX spec_marshal_spec_id_key ON spec(
    marshal_spec_id
);

CREATE TABLE spec_calib (
    id bigint NOT NULL unique,
    dome text NULL,
    bias text NULL,
    flat text NULL,
    cosmic_filter boolean Null,
    drpver text NULL,
    Hg_master text NULL,
    Cd_master text NULL,
    Xe_master text NULL,
    avg_rms float  NULL,
    min_rms float NULL,
    max_rms float NULL
);

-- Table: phot
CREATE TABLE phot (
    id BIGSERIAL,
    phot_calib_id bigint NOT NULL,
    observation_id bigint NOT NULL UNIQUE,
    astrometry boolean  NOT NULL,
    filter text  NOT NULL,
    reducedfile text  NULL,
    sexfile text  NULL,
    maskfile text  NULL,
    pipeline text  NULL,
    marshal_phot_id bigint  NULL,
    CONSTRAINT phot_pk PRIMARY KEY (id)
);

CREATE INDEX phot_marshal_phot_id ON phot(
    marshal_phot_id
);

CREATE TABLE phot_calib (
    id bigint NOT NULL unique,
    bias text NULL,
    flat text NULL
);

-- Table: ref_stars
CREATE TABLE ref_stars (
    id BIGSERIAL,
    phot_id bigint NOT NULL,
    ra decimal(12,6)  NOT NULL,
    dec decimal(12,6)  NOT NULL,
    survey text  NOT NULL,
    filter text  NOT NULL,
    mag decimal(5,2)  NOT NULL,
    magerr decimal(5,2)  NOT NULL,
    instmag decimal(5,2)  NOT NULL,
    istmagerr decimal(5,2)  NOT NULL,
    mjd decimal(8,4)  NOT NULL,
    x int  NOT NULL,
    y int  NOT NULL,
    CONSTRAINT ref_stars_pk PRIMARY KEY (id)
);

CREATE INDEX ref_stars_coordsq3c ON ref_stars(
    q3c_ang2ipix(ra, dec)
);

CREATE INDEX ref_stars_mag_key ON ref_stars(
    mag
);

CREATE INDEX ref_stars_filter_key ON ref_stars(
    filter
);

-- Table: request
CREATE TABLE request (
    id BIGSERIAL,
    object_id bigint NOT NULL,
    user_id bigint NOT NULL,
    program_id bigint  NOT NULL,
    marshal_id bigint NULL,
    exptime integer[]  NOT NULL,
    maxairmass decimal(5,2)  DEFAULT 2.5,
    status text  DEFAULT 'PENDING',
    priority decimal(5,2)  NOT NULL,
    inidate timestamp  NOT NULL,
    enddate timestamp  NOT NULL,
    cadence decimal(5,2) NULL,
    phasesamples decimal(5,2)  NULL,
    sampletolerance decimal(5,2) NULL,
    filters text[] DEFAULT '{ifu,u,g,r,i}',
    nexposures integer[] NULL,
    ordering text[] NULL,
    creationdate timestamp DEFAULT NOW(),
    lastmodified timestamp  DEFAULT NOW(),
    CONSTRAINT request_pk PRIMARY KEY (id)
);

CREATE INDEX request_object_id_key ON request(
    object_id
);

CREATE INDEX request_status_key ON request(
    status
);

CREATE TABLE atomicrequest (
    id BIGSERIAL,
    object_id bigint NOT NULL,
    request_id bigint NOT NULL,
    order_id int NULL,
    exptime float NOT NULL,
    filter text NOT NULL,
    status text DEFAULT 'PENDING',
    priority decimal(5,2)  NOT NULL,
    inidate_plan date  NOT NULL,
    enddate_plan date  NOT NULL,
    inidate timestamp  NOT NULL,
    enddate timestamp  NOT NULL,
    time_elapsed float NULL ,
    creationdate timestamp  DEFAULT NOW(),
    lastmodified timestamp  DEFAULT NOW(),
    CONSTRAINT atomicrequest_pk PRIMARY KEY (id)
);

CREATE INDEX atomicrequest_request_id_key ON atomicrequest(
    request_id
);

CREATE INDEX atomicrequest_object_id_key ON atomicrequest(
    object_id
);

CREATE INDEX atomicrequest_inidate_key ON atomicrequest(
    inidate
);

CREATE INDEX atomicrequest_enddate_key ON atomicrequest(
    enddate
);

-- Table: telescope_stats
CREATE TABLE telescope_stats (
    id BIGSERIAL,
    observation_id bigint NOT NULL UNIQUE,
    date timestamp  NOT NULL,
    dome_status text  NULL,
    in_temp decimal(5,2)  NULL,
    in_humidity decimal(5,2)  NULL,
    in_dew decimal(5,2)  NULL,
    out_temp decimal(5,2)  NULL,
    out_humidity decimal(5,2)  NULL,
    out_dew decimal(5,2)  NULL,
    wind_dir decimal(5,2)  NULL,
    wsp_cur decimal(5,2)  NULL,
    wsp_avg decimal(5,2)  NULL,
    mir_temp decimal(5,2)  NULL,
    top_air decimal(5,2)  NULL,
    pri_temp decimal(5,2)  NULL,
    sec_temp decimal(5,2)  NULL,
    flo_temp decimal(5,2)  NULL,
    bot_temp decimal(5,2)  NULL,
    mid_temp decimal(5,2)  NULL,
    top_temp decimal(5,2)  NULL,
    CONSTRAINT telescope_stats_pk PRIMARY KEY (id)
);

CREATE INDEX telescope_stats_date_key ON telescope_stats(
    date
);

-- Table: users
CREATE TABLE users (
    id BIGSERIAL,
    username text NOT NULL UNIQUE,
    name text  NULL,
    email text  NULL,
    password text NOT NULL,
    CONSTRAINT users_pk PRIMARY KEY (id)
);

-- Table: groups
CREATE TABLE groups (
    id BIGSERIAL,
    designator text NOT NULL UNIQUE,
    CONSTRAINT groups_pk PRIMARY KEY (id)
);

-- Table: groups
CREATE TABLE usergroups (
    user_id bigint NOT NULL,
    group_id bigint NOT NULL,
    CONSTRAINT user_groups PRIMARY KEY (user_id, group_id)
);

-- Table: program
CREATE TABLE program (
    id BIGINT NOT NULL UNIQUE,
    designator text NOT NULL UNIQUE,
    name text NULL ,
    group_id BIGINT NOT NULL,
    PI text NULL,
    priority decimal(5,2) NULL
);

-- Table: allocation
CREATE TABLE allocation (
    pg_designator text NOT NULL,
    inidate DATE NULL,
    enddate DATE NULL,
    time_spent interval NULL,
    time_allocated interval NULL
);

-- foreign keys
-- Reference: atomicrequest_object (table: atomicrequest)
ALTER TABLE usergroups ADD CONSTRAINT usergroups_users
    FOREIGN KEY (user_id)
    REFERENCES users (id)
    NOT DEFERRABLE
    INITIALLY IMMEDIATE
;

ALTER TABLE usergroups ADD CONSTRAINT usergroups_groups
    FOREIGN KEY (group_id)
    REFERENCES groups (id)
    NOT DEFERRABLE
    INITIALLY IMMEDIATE
;

ALTER TABLE program ADD CONSTRAINT program_groups
    FOREIGN KEY (group_id)
    REFERENCES groups (id)
    NOT DEFERRABLE
    INITIALLY IMMEDIATE
;

ALTER TABLE allocation ADD CONSTRAINT allocation_program
    FOREIGN KEY (pg_designator)
    REFERENCES program (designator)
    NOT DEFERRABLE
    INITIALLY IMMEDIATE
;

ALTER TABLE atomicrequest ADD CONSTRAINT atomicrequest_object
    FOREIGN KEY (object_id)
    REFERENCES object (id)
    NOT DEFERRABLE
    INITIALLY IMMEDIATE
;

-- Reference: atomicrequest_request (table: atomicrequest)
ALTER TABLE atomicrequest ADD CONSTRAINT atomicrequest_request
    FOREIGN KEY (request_id)
    REFERENCES request (id)
    NOT DEFERRABLE
    INITIALLY IMMEDIATE
;

-- Reference: classification_spec (table: classification)
ALTER TABLE classification ADD CONSTRAINT classification_spec
    FOREIGN KEY (spec_id)
    REFERENCES spec (id)  
    NOT DEFERRABLE 
    INITIALLY IMMEDIATE
;


-- Reference: classification_spec (table: classification)
ALTER TABLE classification ADD CONSTRAINT classification_object
    FOREIGN KEY (object_id)
    REFERENCES object (id)
    NOT DEFERRABLE
    INITIALLY IMMEDIATE
;

-- Reference: flexure_spec (table: flexure)
ALTER TABLE flexure ADD CONSTRAINT flexure_spec
    FOREIGN KEY (spec_id_1)
    REFERENCES spec (id)  
    NOT DEFERRABLE 
    INITIALLY IMMEDIATE
;

-- Reference: metrics_phot_phot (table: metrics_phot)
ALTER TABLE metrics_phot ADD CONSTRAINT metrics_phot_phot
    FOREIGN KEY (phot_id)
    REFERENCES phot (id)  
    NOT DEFERRABLE 
    INITIALLY IMMEDIATE
;

-- Reference: metrics_spec_spec (table: metrics_spec)
ALTER TABLE metrics_spec ADD CONSTRAINT metrics_spec_spec
    FOREIGN KEY (spec_id)
    REFERENCES spec (id)  
    NOT DEFERRABLE 
    INITIALLY IMMEDIATE
;

-- Reference: phot_phot_calib (table: phot)
ALTER TABLE phot ADD CONSTRAINT phot_phot_calib
    FOREIGN KEY (phot_calib_id)
    REFERENCES phot_calib (id)
    NOT DEFERRABLE
    INITIALLY IMMEDIATE
;

-- Reference: spec_spec_calib (table: spec)
ALTER TABLE spec ADD CONSTRAINT spec_spec_calib
    FOREIGN KEY (spec_calib_id)
    REFERENCES spec_calib (id)
    NOT DEFERRABLE
    INITIALLY IMMEDIATE
;

-- Reference: observation_request (table: observation)
ALTER TABLE observation ADD CONSTRAINT observation_object
    FOREIGN KEY (object_id)
    REFERENCES object (id)
    NOT DEFERRABLE
    INITIALLY IMMEDIATE
;

-- Reference: observation_request (table: observation)
ALTER TABLE observation ADD CONSTRAINT observation_request
    FOREIGN KEY (request_id)
    REFERENCES request (id)
    NOT DEFERRABLE 
    INITIALLY IMMEDIATE
;

-- Reference: observation_atomic_request (table: observation)
ALTER TABLE observation ADD CONSTRAINT observation_atomic_request
    FOREIGN KEY (atomicrequest_id)
    REFERENCES atomicrequest (id)
    NOT DEFERRABLE 
    INITIALLY IMMEDIATE
;

-- Reference: phot_observation (table: phot)
ALTER TABLE phot ADD CONSTRAINT phot_observation
    FOREIGN KEY (observation_id)
    REFERENCES observation (id)  
    NOT DEFERRABLE 
    INITIALLY IMMEDIATE
;

-- Reference: ref_stars_phot (table: ref_stars)
ALTER TABLE ref_stars ADD CONSTRAINT ref_stars_phot
    FOREIGN KEY (phot_id)
    REFERENCES phot (id)  
    NOT DEFERRABLE 
    INITIALLY IMMEDIATE
;

-- Reference: requests_objects (table: request)
ALTER TABLE request ADD CONSTRAINT requests_objects
    FOREIGN KEY (object_id)
    REFERENCES object (id)  
    NOT DEFERRABLE 
    INITIALLY IMMEDIATE
;

-- Reference: requests_program (table: request)
ALTER TABLE request ADD CONSTRAINT program_request
    FOREIGN KEY (program_id)
    REFERENCES program (id)
    NOT DEFERRABLE
    INITIALLY IMMEDIATE
;

-- Reference: spec_flexure (table: flexure)
ALTER TABLE flexure ADD CONSTRAINT spec_flexure
    FOREIGN KEY (spec_id_2)
    REFERENCES spec (id)  
    NOT DEFERRABLE 
    INITIALLY IMMEDIATE
;

-- Reference: spec_observation (table: spec)
ALTER TABLE spec ADD CONSTRAINT spec_observation
    FOREIGN KEY (observation_id)
    REFERENCES observation (id)  
    NOT DEFERRABLE 
    INITIALLY IMMEDIATE
;

-- Reference: telescope_stats_observation (table: telescope_stats)
ALTER TABLE telescope_stats ADD CONSTRAINT telescope_stats_observation
    FOREIGN KEY (observation_id)
    REFERENCES observation (id)  
    NOT DEFERRABLE 
    INITIALLY IMMEDIATE
;

-- Reference: users_request (table: request)
ALTER TABLE request ADD CONSTRAINT users_request
    FOREIGN KEY (user_id)
    REFERENCES users (id)  
    NOT DEFERRABLE 
    INITIALLY IMMEDIATE
;

-- Reference: object_elliptical_heliocentric (table: elliptical_heliocentric)
ALTER TABLE elliptical_heliocentric ADD CONSTRAINT object_elliptical_heliocentric
    FOREIGN KEY (object_id)
    REFERENCES object (id)
    NOT DEFERRABLE
    INITIALLY IMMEDIATE
;

-- Reference: object_earth_satellite (table: earth_satellite)
ALTER TABLE earth_satellite ADD CONSTRAINT object_earth_satellite
    FOREIGN KEY (object_id)
    REFERENCES object (id)
    NOT DEFERRABLE
    INITIALLY IMMEDIATE
;

-- Reference: object_parabolic_heliocentric (table: parabolic_heliocentric)
ALTER TABLE parabolic_heliocentric ADD CONSTRAINT object_parabolic_heliocentric
    FOREIGN KEY (object_id)
    REFERENCES object (id)
    NOT DEFERRABLE
    INITIALLY IMMEDIATE
;

-- Reference: object_hyperbolic_heliocentric (table: hyperbolic_heliocentric)
ALTER TABLE hyperbolic_heliocentric ADD CONSTRAINT object_hyperbolic_heliocentric
    FOREIGN KEY (object_id)
    REFERENCES object (id)
    NOT DEFERRABLE
    INITIALLY IMMEDIATE
;
-- Reference: object_periodic (table: periodic)
ALTER TABLE periodic ADD CONSTRAINT object_periodic
    FOREIGN KEY (object_id)
    REFERENCES object (id)
    NOT DEFERRABLE
    INITIALLY IMMEDIATE
;

-- End of file.

