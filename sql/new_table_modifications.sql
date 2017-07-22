CREATE TABLE spec_calib (
    spec_id bigint NOT NULL unique,
    dome text NULL,
    bias text NULL,
    flat text NULL,
    cosmic_filter boolean Null,
    DRPVER text NULL,
    Hg_master text NULL,
    Cd_master text NULL,
    Xe_master text NULL,
    avg_rms text NULL,
    min_rms decimal(5,2) NULL,
    max_rms decimal(5,2) NULL
);

CREATE TABLE phot_calib (
    phot_id bigint NOT NULL unique,
    bias text NULL,
    flat text NULL
);

-- Table: program
CREATE TABLE program (
    id BIGINT NOT NULL UNIQUE,
    designator text NOT NULL UNIQUE,
    name text NULL ,
    group_id BIGINT NOT NULL,
    PI text NULL,
    time_allocated interval NULL,
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

ALTER TABLE phot_calib ADD CONSTRAINT phot_phot_id
    FOREIGN KEY (phot_id)
    REFERENCES phot (id)
    NOT DEFERRABLE
    INITIALLY IMMEDIATE
;

ALTER TABLE spec_calib ADD CONSTRAINT spec_spec_id
    FOREIGN KEY (spec_id)
    REFERENCES spec (id)
    NOT DEFERRABLE
    INITIALLY IMMEDIATE
;

