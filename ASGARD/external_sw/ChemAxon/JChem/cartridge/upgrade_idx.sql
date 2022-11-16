spool spoolfile
set serveroutput on
declare
  ityp_owner VARCHAR2(30) := 'PKOVACS_TRUNK';
  upgrade_only VARCHAR2(300) := ' parameters(''upgradeOnly=y'')';
  cursor c1 (itowner VARCHAR2) is
    select owner, index_name, partitioned
          from dba_indexes
          where ityp_owner = itowner and ityp_name = 'JC_IDXTYPE';
  cursor c2 (owner VARCHAR2, idxname VARCHAR2) is
    select partition_name
          from dba_ind_partitions
          where index_owner = owner and index_name = idxname;
  v_sql VARCHAR2(4000);
  idxtbl_entry integer;
begin
  execute immediate 'call ' || ityp_owner || '.jchem_core_pkg.init()';
  for idx in c1(ityp_owner)
  loop
    if idx.partitioned = 'YES' then
      for part in c2(idx.owner, idx.index_name)
      loop
        execute immediate 'select count(prop_value) from ' || idx.owner || '.jc_idx_property '
                || 'where prop_name = ''' || idx.owner || '.' || idx.index_name
                || '_' || part.partition_name || '.idxTable'''
                into idxtbl_entry;
        if idxtbl_entry > 0 then -- regular structure table
          v_sql := 'alter index '
                || idx.owner || '.' || idx.index_name
                || ' rebuild partition ' || part.partition_name
                || upgrade_only;
          dbms_output.put_line(v_sql);
          execute immediate v_sql;
        end if;
      end loop;
    else
      execute immediate 'select count(prop_value) from ' || idx.owner || '.jc_idx_property '
              || 'where prop_name = ''' || idx.owner || '.' || idx.index_name || '.idxTable'''
              into idxtbl_entry;
      if (idxtbl_entry > 0) then -- regular structure table
        v_sql := 'alter index '
              || idx.owner || '.' || idx.index_name || ' rebuild'
              || upgrade_only;
        dbms_output.put_line(v_sql);
        execute immediate v_sql;
      end if;
    end if;
  end loop;
end;
/
quit
