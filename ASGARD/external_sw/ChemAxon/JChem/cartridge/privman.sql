create or replace package privman_pkg authid current_user as
  procedure grants_on_jcobjs(jchemowner varchar2, user_name varchar2);
  procedure grants_on_jcobjs(jchemowner varchar2, user_name varchar2,
                             jcproptab_name varchar2,
                             can_create_index number);
  procedure syns_for_jcobjs(jchemowner varchar2);
  procedure syns_for_jcobjs(jchemowner varchar2, user_name varchar2);
  procedure public_syns_for_jcobjs(jchemowner varchar2);
  procedure drop_syns_for_jcobjs(jchemowner varchar2);
  procedure drop_syns_for_jcobjs(jchemowner varchar2, user_name varchar2);
  procedure drop_public_syns_for_jcobjs(jchemowner varchar2);
  function get_synctrl_sql(jchemowner varchar2, dro number := 0,
                           pub varchar2 := '', user_name varchar := '') return varchar2;
  procedure grants_on_jcidx(user_name varchar,
                            idx_schema varchar2, idx_name varchar2,
                            with_search number, with_insert number,
                            with_update number, with_delete number);
  procedure set_jcpt_apart(boo number);
  procedure set_logtable_name(ltbl_name varchar2);
  procedure set_dryrun(drun number);
end privman_pkg;
/
show errors

create or replace package body privman_pkg is
  jcpt_apart integer := 0;
  logtable_name varchar2(61);
  dry_run integer := 0;

  procedure set_jcpt_apart(boo number) is
  begin
    jcpt_apart := boo;
  end;

  procedure set_logtable_name(ltbl_name varchar2) is
  begin
    logtable_name := ltbl_name;
  end;
  
  procedure set_dryrun(drun number) is
  begin
    dry_run := drun;
  end;

  procedure log(stmt varchar2) is
  begin
      execute immediate 'insert into '
                        || logtable_name 
                        || ' (stmt) values(:a)' using stmt;
      commit;
  end;

  procedure exec(stmt varchar2) is
  begin
    if dry_run > 0 then
      dbms_output.put_line('>>>>> ' || stmt);
    end if;

    if logtable_name is not null then
      log(stmt);
    end if;

    if dry_run = 0 then
      execute immediate stmt;
    end if;
  end;

  procedure grants_on_jcobjs(jchemowner varchar2, user_name varchar2) as
  begin
    grants_on_jcobjs(jchemowner, user_name, null, null);
  end;

  procedure grants_on_jcobjs(jchemowner varchar2, user_name varchar2,
                             jcproptab_name varchar2,
                             can_create_index number) as
    v_sql varchar(4000);
    cursor c1 is
    select 'grant execute on ' 
            || a.OWNER || '."' || a.OBJECT_NAME 
            || '" to ' || user_name as grant_text
        from all_objects a
        where lower(a.OWNER) = lower(jchemowner)
            and a.OBJECT_TYPE in ('PACKAGE',
                                  'FUNCTION',
                                  'OPERATOR',
                                  'TYPE',
                                  'INDEXTYPE',
                                  'PROCEDURE',
                                  'PACKAGE BODY',
                                  'JAVA CLASS',
                                  'TYPE BODY')
    union all
    select 'grant select on '
              || a.OWNER || '.' || a.OBJECT_NAME
              || ' to ' || user_name as grant_text
        from all_objects a
        where lower(a.OWNER) = lower(jchemowner)
            and a.OBJECT_TYPE in ('SEQUENCE', 'TABLE')
            and lower(a.OBJECT_NAME) NOT LIKE 'bin%';
  begin
    exec('grant select on ' || jchemowner || '.jc_idx_property to ' || user_name);
    if (jcproptab_name is not null) then
      begin
        exec('grant select on ' || jcproptab_name || ' to ' || user_name);
        if (can_create_index > 0) then
          exec('grant insert on ' || jcproptab_name || ' to ' || user_name);
          exec('grant update on ' || jcproptab_name || ' to ' || user_name);
          exec('grant delete on ' || jcproptab_name || ' to ' || user_name);
        end if;
      exception
      when others then
        if sqlcode <> -942 then
          raise;
        end if;
      end;
    end if;

    for t in c1 loop
      v_sql := t.grant_text;
      exec(v_sql);
    end loop;
  end;

  function get_synctrl_sql(jchemowner varchar2, dro number := 0,
                           pub varchar2 := '', user_name varchar := '') return varchar2 as
    v_sql varchar2(4000);

    cuser varchar2(61);
    cuserprefix varchar2(62) := '';
    action varchar2(100) := 'create or replace';

    clau varchar2(4000) := 'from all_objects a '
        || 'where lower(a.OWNER) = lower(''' || jchemowner || ''') '
        ||    'and a.OBJECT_TYPE in (''PACKAGE'', '
        ||                          '''FUNCTION'', '
        ||                          '''OPERATOR'', '
        ||                          '''TYPE'', '
        ||                          '''INDEXTYPE'', '
        ||                          '''PROCEDURE'') '
        ||    'and not ((a.OBJECT_TYPE = ''PACKAGE'' '
        ||             'or a.OBJECT_TYPE = ''PACKAGE BODY'') '
        ||             'and a.OBJECT_NAME = ''PRIVMAN_PKG'') ';
  begin
    if user_name is not null then
	  cuserprefix := user_name || '.';
    elsif pub is null then 
      execute immediate 'SELECT USER FROM DUAL' into cuser;
      cuserprefix := cuser || '.';
    end if;

    if dro > 0 then
      action := 'drop';
    end if;

    v_sql := 'select ''' || action || ' ' || pub || ' synonym '
          || cuserprefix || ''' || a.OBJECT_NAME ';

    if dro = 0 then
      v_sql := v_sql || '|| '' for '' || a.OWNER || ''.'' || a.OBJECT_NAME ';
    end if;

    v_sql := v_sql || ' as syn_text ';
    return v_sql || clau;
  end; -- get_synctrl_sql

  procedure synctrl_for_jcobjs(jchemowner varchar2, dro number := 0,
                               pub varchar2 := '', user_name varchar2 := '') as
    v_sql varchar2(4000);
    v_sql1 varchar2(4000);
    TYPE cursor_ref IS REF CURSOR;
    c1 cursor_ref;

  begin
    v_sql := get_synctrl_sql(jchemowner, dro, pub, user_name);

    if dry_run > 0 then
      dbms_output.put_line('>>>>> ' || v_sql);
    end if;

    open c1 for v_sql;
    loop
        fetch c1 into v_sql1;
        exit when c1%notfound;

        dbms_output.put_line('>>>>> ' || v_sql1);
        begin
          exec(v_sql1);
        exception
        when others then
          if dro > 0 then
            -- Ignore synonym not found error
            if sqlcode not in (-1432,-1434) then
              raise;
            end if;
          end if;
        end;
    end loop;
    close c1;
  end;

  procedure syns_for_jcobjs(jchemowner varchar2) as
  begin
    synctrl_for_jcobjs(jchemowner);
  end;

  procedure syns_for_jcobjs(jchemowner varchar2,user_name varchar2) as
  begin
    synctrl_for_jcobjs(jchemowner=>jchemowner,user_name=>user_name);
  end;

  procedure public_syns_for_jcobjs(jchemowner varchar2) as
  begin
    synctrl_for_jcobjs(jchemowner, 0, 'public');
  end;

  procedure drop_syns_for_jcobjs(jchemowner varchar2) as
  begin
    synctrl_for_jcobjs(jchemowner, 1);
  end;

  procedure drop_syns_for_jcobjs(jchemowner varchar2,user_name varchar2) as
  begin
    synctrl_for_jcobjs(jchemowner=>jchemowner,user_name=>user_name,dro=>1);
  end;

  procedure drop_public_syns_for_jcobjs(jchemowner varchar2) as
  begin
    synctrl_for_jcobjs(jchemowner, 1, 'public');
  end;

  procedure grants_on_jcidx(user_name varchar,
                            idx_schema varchar2, idx_name varchar2,
                            with_search number, with_insert number,
                            with_update number, with_delete number) as
    table_owner varchar(30);
    table_name varchar(62);
    jcpt_name varchar(62);
    idxtbl_qname varchar(62);
    v_sql varchar2(4000);
    err_num number;
    prop_name varchar2(4000);
  begin
      dbms_output.put_line('>>>>> ' || idx_schema || ', ' || idx_name);
      v_sql := 'select table_owner, table_name '
                || 'from all_indexes where owner = :a and index_name = :b';

      begin
      dbms_output.put_line('>>>>> Executing ' || v_sql);
      execute immediate v_sql into table_owner, table_name
                              using upper(idx_schema), upper(idx_name);
      exception
        when others then
          if sqlcode = 100 then
            raise_application_error(-20102, 'No data found for '''
                                    || v_sql || '''');
          else
            raise;
          end if;
      end;

      dbms_output.put_line('>>>>> table_owner: ' || table_owner
                          || ', table_name: ' || table_name);

      v_sql := 'select upper(prop_value) from ' 
                || idx_schema || '.jc_idx_property where upper(prop_name) = :a';
      begin
        prop_name := upper(idx_schema) || '.' || upper(idx_name) || '.JCHEMPROPERTIES';
        dbms_output.put_line('>>>>> Executing ' || v_sql || ' with '
                             || prop_name || ' as :a');
        execute immediate v_sql into jcpt_name using prop_name;

        prop_name := upper(idx_schema) || '.' || upper(idx_name) || '.IDXTABLE';
        dbms_output.put_line('>>>>> Executing ' || v_sql || ' with '
                             || prop_name || ' as :a');
        execute immediate v_sql into idxtbl_qname using prop_name;
      exception
        when others then
          err_num := sqlcode;
          if (err_num = 1403 or err_num = 100) then
            idxtbl_qname := table_owner || '.' || table_name;
          else
            raise;
          end if;
      end;

      if jcpt_apart = 0 then 
        exec('grant select on ' || jcpt_name || ' to ' || user_name);
        -- FS#12761
        exec('grant update on ' || jcpt_name || ' to ' || user_name);
      end if;

      exec('grant select on ' || table_owner || '.JC_IDX_UDOP to ' || user_name);

      exec('grant select on ' || table_owner || '.jc_idx_property to ' || user_name);

      if (with_search > 0) then
        exec('grant select on ' || table_owner || '.' || table_name || ' to ' || user_name);

        exec('grant select on ' || idxtbl_qname || ' to ' || user_name);
        -- for reading update logs
        exec('grant select on ' || idxtbl_qname || '_UL to ' || user_name);
        -- for updating update logs
        exec('grant update on ' || idxtbl_qname || '_UL to ' || user_name);
        -- for deleting update logs:
        exec('grant delete on ' || idxtbl_qname || '_UL to ' || user_name);

        -- FIXME: we may not need SELECT privilege for this case
        exec('grant select on ' || jcpt_name || '_CR to ' || user_name);
        -- StructureCache will want to update the cache registration timestamp
        exec('grant update on ' || jcpt_name || '_CR to ' || user_name);
        -- StructureCache will want to delete expired cache registrations
        exec('grant delete on ' || jcpt_name || '_CR to ' || user_name);
      end if;

      if (with_insert > 0) then
        exec('grant insert on ' || table_owner || '.' || table_name || ' to ' || user_name);
        -- needed for LOB structure columns
        exec('grant select on ' || table_owner || '.' || table_name || ' to ' || user_name);
        -- needed for LOB structure columns
        exec('grant update on ' || table_owner || '.' || table_name || ' to ' || user_name);

        exec('grant select on ' || idxtbl_qname || '_UL to ' || user_name);
        exec('grant insert on ' || idxtbl_qname || '_UL to ' || user_name);
        exec('grant update on ' || idxtbl_qname || '_UL to ' || user_name);
        exec('grant delete on ' || idxtbl_qname || '_UL to ' || user_name);

        -- We well want to now which cache ids are on line and requiring
        -- notifications of INSERTs
        exec('grant select on ' || jcpt_name || '_CR to ' || user_name);

        if (upper(idxtbl_qname) = upper(table_owner || '.' || table_name)) then
           exec('grant select on ' || idxtbl_qname || '_sq to ' || user_name);
           exec('grant select on ' || idxtbl_qname || '_usq to ' || user_name);
        end if;
      end if;

      if (with_update > 0) then
        exec('grant select on ' || table_owner || '.' || table_name || ' to ' || user_name);
        exec('grant update on ' || table_owner || '.' || table_name || ' to ' || user_name);

        exec('grant select on ' || idxtbl_qname || '_UL to ' || user_name);
        exec('grant insert on ' || idxtbl_qname || '_UL to ' || user_name);
        exec('grant delete on ' || idxtbl_qname || '_UL to ' || user_name);

        -- We well want to now which cache ids are on line and requiring
        -- notifications of UPDATEs
        exec('grant select on ' || jcpt_name || '_CR to ' || user_name);

        if (upper(idxtbl_qname) = upper(table_owner || '.' || table_name)) then
           exec('grant select on ' || idxtbl_qname || '_usq to ' || user_name);
        end if;
      end if;

      if (with_delete > 0) then
        exec('grant delete on ' || table_owner || '.' || table_name || ' to ' || user_name);

        exec('grant select on ' || idxtbl_qname || '_UL to ' || user_name);
        exec('grant insert on ' || idxtbl_qname || '_UL to ' || user_name);
        exec('grant delete on ' || idxtbl_qname || '_UL to ' || user_name);

        -- We well want to now which cache ids are on line and requiring
        -- notifications of DELETEs 
        exec('grant select on ' || jcpt_name || '_CR to ' || user_name);

        if (upper(idxtbl_qname) = upper(table_owner || '.' || table_name)) then
           exec('grant select on ' || idxtbl_qname || '_usq to ' || user_name);
        end if;
      end if;
  end;

end;
/
show errors

quit
