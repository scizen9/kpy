-- Created by Vertabelo (http://vertabelo.com)
-- Last modification date: 2016-08-19 00:28:45.832

-- tables
-- Table: classification
CREATE TABLE classification (
    id bigint  NOT NULL,
    spec_id bigint  NOT NULL,
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
    id int  NOT NULL,
    rms decimal(8,4)  NOT NULL,
    spec_id_1 bigint  NOT NULL,
    spec_id_2 bigint  NOT NULL,
    CONSTRAINT flexure_pk PRIMARY KEY (id)
);

-- Table: metrics_phot
CREATE TABLE metrics_phot (
    id bigint  NOT NULL,
    phot_id bigint  NOT NULL,
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
    id bigint  NOT NULL,
    spec_id bigint  NOT NULL,
    fwhm decimal(5,2)  NULL,
    background decimal(5,2)  NULL,
    line_fwhm int  NULL,
    CONSTRAINT metrics_spec_pk PRIMARY KEY (id)
);

-- Table: object
CREATE TABLE object (
    id bigint  NOT NULL,
    marshal_id bigint  NOT NULL,
    name text  NOT NULL,
    ra decimal(12,6)  NOT NULL,
    dec decimal(12,6)  NOT NULL,
    CONSTRAINT object_pk PRIMARY KEY (id)
);

-- Table: observation
CREATE TABLE observation (
    id bigint  NOT NULL,
    request_id bigint  NOT NULL,
    mjd decimal(10,2)  NOT NULL,
    airmass decimal(5,2)  NOT NULL,
    exptime decimal(6,2)  NOT NULL,
    file text  NOT NULL,
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

-- Table: phot
CREATE TABLE phot (
    id bigint  NOT NULL,
    observation_id bigint  NOT NULL,
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
    id bigint  NOT NULL,
    phot_id bigint  NOT NULL,
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
    id bigint  NOT NULL,
    object_id bigint  NOT NULL,
    user_id bigint  NOT NULL,
    marshal_id bigint  NOT NULL,
    exptime int  NOT NULL,
    maxairmass decimal(5,2)  NOT NULL,
    isspec boolean  NOT NULL,
    status text  NOT NULL,
    priority decimal(5,2)  NOT NULL,
    inidate date  NOT NULL,
    enddate date  NOT NULL,
    jd0 decimal(5,2) NULL,
    phasefreq decimal(12,6) NULL,
    programid smallint  NOT NULL,
    mag decimal(5,2)  NOT NULL,
    filter text  NOT NULL,
    creationdate date  NULL,
    lastmodified date  NULL,
    CONSTRAINT request_pk PRIMARY KEY (id)
);

-- Table: spec
CREATE TABLE spec (
    id bigint  NOT NULL,
    observation_id bigint  NOT NULL,
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

-- Table: telescope_stats
CREATE TABLE telescope_stats (
    id bigint  NOT NULL,
    observation_id bigint  NOT NULL,
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
    id bigint  NOT NULL,
    name text  NOT NULL,
    email text  NULL,
    CONSTRAINT users_pk PRIMARY KEY (id)
);

-- foreign keys
-- Reference: classification_spec (table: classification)
ALTER TABLE classification ADD CONSTRAINT classification_spec
    FOREIGN KEY (spec_id)
    REFERENCES spec (id)  
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

-- End of file.

