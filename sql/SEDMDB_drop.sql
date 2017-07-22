-- Created by Vertabelo (http://vertabelo.com)
-- Last modification date: 2016-08-19 00:28:45.832

-- foreign keys
ALTER TABLE classification
    DROP CONSTRAINT classification_spec;

ALTER TABLE flexure
    DROP CONSTRAINT flexure_spec;

ALTER TABLE metrics_phot
    DROP CONSTRAINT metrics_phot_phot;

ALTER TABLE metrics_spec
    DROP CONSTRAINT metrics_spec_spec;

ALTER TABLE spec_calib
    DROP CONSTRAINT spec_spec_id;

ALTER TABLE phot_calib
    DROP CONSTRAINT phot_phot_id;

ALTER TABLE observation
    DROP CONSTRAINT observation_request;

ALTER TABLE phot
    DROP CONSTRAINT phot_observation;

ALTER TABLE ref_stars
    DROP CONSTRAINT ref_stars_phot;

ALTER TABLE request
    DROP CONSTRAINT requests_objects;

ALTER TABLE flexure
    DROP CONSTRAINT spec_flexure;

ALTER TABLE spec
    DROP CONSTRAINT spec_observation;

ALTER TABLE telescope_stats
    DROP CONSTRAINT telescope_stats_observation;

ALTER TABLE request
    DROP CONSTRAINT users_request;

ALTER TABLE program
    DROP CONSTRAINT program_groups;

ALTER TABLE allocation
    DROP CONSTRAINT allocation_program;


-- tables

DROP TABLE usergroups;

DROP TABLE classification;

DROP TABLE flexure;

DROP TABLE metrics_phot;

DROP TABLE metrics_spec;

DROP TABLE spec_calib;

DROP TABLE phot_calib;

DROP TABLE observation;

DROP TABLE phot;

DROP TABLE ref_stars;

DROP TABLE atomicrequest;

DROP TABLE spec;

DROP TABLE telescope_stats;

DROP TABLE allocation;

DROP TABLE program;

DROP TABLE users;

DROP TABLE elliptical_heliocentric;

DROP TABLE hyperbolic_heliocentric;

DROP TABLE parabolic_heliocentric;

DROP TABLE earth_satellite;

DROP TABLE periodic;

DROP TABLE object;

DROP TABLE request;

DROP TABLE groups;

-- End of file.

