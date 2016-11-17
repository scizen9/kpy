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
    redshifterr decimal(7,5)  NULL,
    phase decimal(5,2)  NULL,
    phase_err decimal(5,2)  NULL,
    classifier text  NULL,
    score decimal(5,2)  NULL,
    CONSTRAINT classification_pk PRIMARY KEY (id)
);

-- Table: flexure
CREATE TABLE flexure (
    id BIGSERIAL,
    rms decimal(8,4)  NOT NULL,
    spec_id_1 bigint  NOT NULL,
    spec_id_2 bigint  NOT NULL,
    CONSTRAINT flexure_pk PRIMARY KEY (id)
);

-- Table: metrics_phot
CREATE TABLE metrics_phot (
    id BIGSERIAL,
    phot_id bigint NOT NULL,
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
    spec_id bigint NOT NULL,
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
    marshal_id bigint NOT NULL UNIQUE,
    name text  NOT NULL,
    ra decimal(12,6)  NOT NULL,
    dec decimal(12,6)  NOT NULL,
    typedesig varchar(1),
    epoch float,
    CONSTRAINT object_pk PRIMARY KEY (id)
);


--parameters from http://www.clearskyinstitute.com/xephem/help/xephem.html#mozTocId468501
CREATE TABLE elliptical_heliocentric (
    id bigint not null,
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
    s decimal(10,8),
    CONSTRAINT sso_pk PRIMARY KEY (id)

);


CREATE TABLE hyperbolic_heliocentric (
    id bigint not null,
    object_id bigint not null UNIQUE,
    -- date of the epoch of perihelion
    T date,
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
    -- first abd second components of magnitude model
    M1 decimal(5,4),
    M2 decimal(5,4),
    -- angular size at 1 AU
    s decimal(10,8),
    CONSTRAINT hyperbolic_heliocentric_pk PRIMARY KEY (id)

);


CREATE TABLE parabolic_heliocentric (
    id bigint not null,
    object_id bigint not null UNIQUE,
    -- date of the epoch of perihelion
    T date,
    -- eccentricity (<1)
    e decimal(10,8),
    --inclination
    inclination decimal(10,8),
    -- longitude of ascending node
    longascnode_O decimal(10,8),
    -- argument of perihelion
    perihelion_o decimal(10,8),
    -- perihelion distance, AU
    q decimal(10,8),
    -- equinox year
    D int,
    -- first abd second components of magnitude model
    M1 decimal(5,4),
    M2 decimal(5,4),
    -- angular size at 1 AU
    s decimal(10,8),
    CONSTRAINT parabolic_heliocentric_pk PRIMARY KEY (id)
);


CREATE TABLE earth_satellite (
    id bigint not null,
    object_id bigint not null UNIQUE,
    -- first date the elements are valid
    T date,
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
    id bigint not null,
    object_id bigint not null,
    mjd0 decimal (10,8),
    phasedays decimal(10,8),
    CONSTRAINT periodic_pk PRIMARY KEY (id)


);

-- Table: observation
CREATE TABLE observation (
    id BIGSERIAL,
    object_id bigint NOT NULL,
    request_id bigint NOT NULL,
    mjd decimal(10,2)  NOT NULL,
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
    utc date  NOT NULL,
    dec_off decimal(5,2)  NOT NULL,
    camera text  NULL,
    CONSTRAINT observation_pk PRIMARY KEY (id)
);

-- Table: spec
CREATE TABLE spec (
    id BIGSERIAL,
    observation_id bigint NOT NULL UNIQUE,
    reducedfile text  NULL,
    sexfile text  NULL,
    biasfile text  NULL,
    flatfile text  NULL,
    imgset text  NULL,
    quality int  NOT NULL,
    cubefile text  NULL,
    standardfile text  NULL,
    marshal_spec_id bigint  NULL,
    skysub boolean  NOT NULL,
    CONSTRAINT spec_pk PRIMARY KEY (id)
);

-- Table: phot
CREATE TABLE phot (
    id BIGSERIAL,
    observation_id bigint NOT NULL UNIQUE,
    astrometry boolean  NOT NULL,
    filter text  NOT NULL,
    reducedfile text  NULL,
    sexfile text  NULL,
    biasfile text  NULL,
    maskfile text  NULL,
    flatfile text  NULL,
    pipeline text  NULL,
    marshal_phot_id bigint  NULL,
    CONSTRAINT phot_pk PRIMARY KEY (id)
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

-- Table: request
CREATE TABLE request (
    id BIGSERIAL,
    object_id bigint NOT NULL,
    user_id bigint NOT NULL,
    program_id smallint  NOT NULL,
    marshal_id bigint NULL,
    exptime int  NOT NULL,
    maxairmass decimal(5,2)  DEFAULT 2.5,
    status text  DEFAULT 'PENDING',
    priority decimal(5,2)  NOT NULL,
    inidate date  NOT NULL,
    enddate date  NOT NULL,
    cadence decimal(5,2) NULL,
    phasesamples decimal(5,2)  NULL,
    sampletolerance decimal(5,2) NULL,
    filters text[]  NULL,
    nexposures integer[] NULL,
    ordering integer[] NULL,
    creationdate date DEFAULT NOW(),
    lastmodified date  DEFAULT NOW(),
    CONSTRAINT request_pk PRIMARY KEY (id)
);

CREATE TABLE atomicrequest (
    id BIGSERIAL,
    object_id bigint NOT NULL,
    request_id bigint NOT NULL,
    order_id int NULL,
    exptime float  NOT NULL,
    instrument text  NOT NULL,
    status text  NOT NULL,
    active boolean DEFAULT FALSE,
    priority decimal(5,2)  NOT NULL,
    inidate date  NOT NULL,
    enddate date  NOT NULL,
    filter text NULL,
    creationdate date  DEFAULT NOW(),
    lastmodified date  DEFAULT NOW(),
    CONSTRAINT atomicrequest_pk PRIMARY KEY (id)
);

-- Table: telescope_stats
CREATE TABLE telescope_stats (
    id BIGSERIAL,
    observation_id bigint NOT NULL,
    date date  NOT NULL,
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

-- Table: users
CREATE TABLE users (
    id BIGSERIAL,
    group_id bigint NOT NULL,
    name text  NOT NULL,
    email text  NULL,
    CONSTRAINT users_pk PRIMARY KEY (id)
);

-- Table: groups
CREATE TABLE groups (
    id BIGSERIAL,
    designator text NULL UNIQUE,
    CONSTRAINT groups_pk PRIMARY KEY (id)
);

-- Table: groups
CREATE TABLE usergroups (
    user_id bigint NOT NULL,
    group_id bigint NOT NULL,
    CONSTRAINT user_groups PRIMARY KEY (user_id, group_id)
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

-- Reference: observation_request (table: observation)
ALTER TABLE observation ADD CONSTRAINT observation_request
    FOREIGN KEY (request_id)
    REFERENCES request (id)  
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

-- Reference: user_group (table: users)
ALTER TABLE users ADD CONSTRAINT user_group
    FOREIGN KEY (group_id)
    REFERENCES groups (id)
    NOT DEFERRABLE
    INITIALLY IMMEDIATE
;

-- End of file.

