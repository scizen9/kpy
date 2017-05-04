#In yupana
psql -h localhost -d ptftransient -U ptftransient -c "\copy (SELECT id, name, iauname, ra, dec, 'f' as typedesig, 2000 as epoch, creationdate  from sources s limit 3) To '/tmp/test.csv' With CSV"

#Copy the file to pharos

#In pharos

psql -h localhost -d sedmdbtest -U sedmadmin -c "\copy object(marshal_id,name,iauname,ra,dec,typedesig,epoch,creationdate) from '/tmp/test.csv' DELIMITER ',' CSV HEADER"

