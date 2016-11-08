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

-- tables
DROP TABLE classification;

DROP TABLE flexure;

DROP TABLE metrics_phot;

DROP TABLE metrics_spec;

DROP TABLE object;

DROP TABLE observation;

DROP TABLE phot;

DROP TABLE ref_stars;

DROP TABLE request;

DROP TABLE schedule;

DROP TABLE spec;

DROP TABLE telescope_stats;

DROP TABLE users;

-- End of file.

