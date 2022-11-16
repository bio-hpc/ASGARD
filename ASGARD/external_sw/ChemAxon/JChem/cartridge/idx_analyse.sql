drop type jcc_nonidxed_array;
drop type jcc_name_value_pair_array;

whenever sqlerror exit 1

create or replace type jcc_nonidxed_record as object (
    rid varchar2(61), smiles varchar2(32767)
    );
/
show errors;

create type jcc_nonidxed_array is varray(20000000) of jcc_nonidxed_record;
/
show errors;

create or replace type jcc_name_value_pair as object (
    name varchar2(32767), value varchar2(32767)
    );
/
show errors;

create or replace type jcc_name_value_pair_array is varray(20000000) of jcc_name_value_pair;
/
show errors;

create or replace package jcc_idx_stats authid current_user as
  function get_idxname(tabname varchar2, colname varchar2) return varchar2;
  function get_prim_idx_attr(idxname varchar2, attrname varchar2)
      return varchar2;
  function get_idxtabname(idxname varchar2) return varchar2;
  function get_jcptname(idxname varchar2) return varchar2;
  function get_g_jcprop(idxname varchar2, propname varchar2) return varchar2;
  function get_jcprop(idxname varchar2, propname varchar2) return varchar2;
  function get_all_jcprop(idxname varchar2)
      return jcc_name_value_pair_array pipelined;
  function get_cacheregtabname(idxname varchar2) return varchar2;
  procedure list_cacheids(idxname varchar2);
  function t_nonidxed(tabname varchar2, colname varchar2)
      return jcc_nonidxed_array pipelined;
  function get_avg_cdsmiles_len(idxname varchar2) return number;
end jcc_idx_stats;
/
show errors


create or replace package body jcc_idx_stats as

  function get_idxname(tabname varchar2, colname varchar2) return varchar2 as
    idxname varchar2(61);
  begin
    for rec in (select a.index_name
                from all_ind_columns a, all_indexes b
                where a.table_owner = user
                  and a.table_name = upper(tabname)
                  and a.column_name = upper(colname)
                  and a.index_name = b.index_name
                  and b.ityp_name = 'JC_IDXTYPE')
    loop
      idxname := rec.index_name;
    end loop;
    if idxname is null then
      raise_application_error(-20001, 'Index name not found for '
          || tabname || '.' || colname);
    end if;
    idxname := user || '.' || idxname;
    dbms_output.put_line('=====> Index name: ' || idxname);
    return idxname;
  end;

  function get_prim_idx_attr(idxname varchar2, attrname varchar2)
      return varchar2 as
    type mycurtype is ref cursor;
    mycur mycurtype;
    s varchar(32767);

    attrval varchar(32767) := null;
  begin
    s := 'select prop_value
                from jc_idx_property
                where prop_name = '''
                || idxname || '.' || attrname || '''';
    dbms_output.put_line('=====> Opening cursor for : ' || s);
    open mycur for s;
    loop
      fetch mycur into attrval;
      exit when mycur%notfound;
    end loop;
    dbms_output.put_line('=====> Master index attribute for '
        || idxname || '.' || attrname || ': ' || attrval);
    return attrval;
  end;

  function get_idxtabname(idxname varchar2) return varchar2 as
    idxtabname varchar2(61);
  begin
    idxtabname := get_prim_idx_attr(idxname, 'idxTable');
    if idxtabname is null then
      raise_application_error(-20001, 'Index table not found for '
          || idxname);
    end if;

    dbms_output.put_line('=====> Index table name: ' || idxtabname);
    return idxtabname;
  end;

  function get_jcptname(idxname varchar2) return varchar2 as
    jcptname varchar2(61);
    notfound exception;
  begin
    jcptname := get_prim_idx_attr(idxname, 'JChemProperties');
    if jcptname is null then
      raise_application_error(-20001, 'JChem Properties table not found for '
          || idxname);
    end if;
    dbms_output.put_line('=====> JChem Properties table name: ' || jcptname);
    return jcptname;
  end;

  function get_g_jcprop(idxname varchar2, propname varchar2) return varchar2 as
    jcptname varchar2(61);
    propval varchar2(32767);
    sqlstr varchar2(32767);
    fullpropname varchar2(1000);
  begin
    jcptname := get_jcptname(idxname);
    sqlstr := 'select prop_value from ' || jcptname || ' where prop_name = :a';
    dbms_output.put_line('===> Executing ' || sqlstr
            || ' [:a: ' || propname || ']');
    execute immediate sqlstr into propval using propname;
    return propval;
  end;

  function get_jcprop(idxname varchar2, propname varchar2) return varchar2 as
    jcptname varchar2(61);
    propval varchar2(32767);
    sqlstr varchar2(32767);
    fullpropname varchar2(1000);
  begin
    jcptname := get_jcptname(idxname);
    sqlstr := 'select prop_value from ' || jcptname || ' where prop_name = :a';
    fullpropname := 'idxtable.' || get_idxtabname(idxname) || '.' || propname;
    dbms_output.put_line('===> Executing ' || sqlstr
            || ' [:a: ' || fullpropname || ']');
    execute immediate sqlstr into propval using fullpropname;
    return propval;
  end;

  function get_all_jcprop(idxname varchar2) return
      jcc_name_value_pair_array pipelined as
    type mycurtype is ref cursor;
    mycur mycurtype;

    s varchar(32000);
    jcptname varchar(61);
    idxtabname varchar(61);

    name varchar(4000);
    value varchar(4000);
    value_ext long raw;
    bufext varchar2(32767);
    extblsize integer := 4000;
    extblcount integer;
    substr_start integer;
    substr_len integer;
  begin
    jcptname := get_jcptname(idxname);
    idxtabname := get_idxtabname(idxname);
    s := 'select prop_name, prop_value, prop_value_ext '
          || ' from ' || jcptname
          || ' where lower(prop_name) like ''%table.'
          || lower(idxtabname) || '.%''';
    dbms_output.put_line('=====> Opening cursor for : ' || s);
    open mycur for s;
    loop
      fetch mycur into name, value, value_ext;
      exit when mycur%notfound;
      if value is null then
        dbms_output.put_line('==========> length(value_ext): '
            || length(value_ext));
        if length(value_ext) < (32767 - extblsize) then
          dbms_output.put_line('==========> length(value_ext)/4: '
              || length(value_ext)/extblsize);
          extblcount := length(value_ext)/extblsize;
        else
          extblcount := (32767 - extblsize)/extblsize;
        end if;
        dbms_output.put_line('=========> extblcount: ' || extblcount);
        for n in 0..extblcount-1 loop
          substr_start := 1+n*extblsize;
          substr_len := extblsize;
          bufext := utl_raw.cast_to_varchar2(substrb(value_ext, substr_start, substr_len));
          pipe row(jcc_name_value_pair(name, bufext));
        end loop;
        substr_start := substr_start + extblsize;
        substr_len := mod(length(value_ext), extblsize);
        dbms_output.put_line('==========> substr_start: ' || substr_start || ', substr_len: ' || substr_len);
        bufext := utl_raw.cast_to_varchar2(substrb(value_ext, substr_start, substr_len));
        pipe row(jcc_name_value_pair(name, bufext));
      else
        pipe row(jcc_name_value_pair(name, value));
      end if;
    end loop;

    return;
  end;

  function get_cacheregtabname(idxname varchar2) return varchar2 as
    jcptname varchar2(61);
    cacheregtabname varchar2(61);
  begin
    jcptname := get_jcptname(idxname);
    cacheregtabname := get_g_jcprop(idxname, 'cache.registration_table');
    dbms_output.put_line('=====> Cache registraton table name: '
          || cacheregtabname);
    return cacheregtabname;
  end;

  procedure list_cacheids(idxname varchar2) as
    type mycurtype is ref cursor;
    mycur mycurtype;
    cachereg_tabname varchar2(61);
    s varchar2(32767);
    cache_id varchar2(4000);
    registration_time varchar2(4000);
    is_protected integer;
  begin
    cachereg_tabname := jcc_idx_stats.get_cacheregtabname(user
                                      || '.' || upper(idxname));
    dbms_output.put_line('====================');
    dbms_output.put_line('Registered cache ids:');
    dbms_output.put_line('>>>>>>>>>>>>>>>>>>>>');
    s := 'select * from ' || cachereg_tabname;
    open mycur for s;
    loop
      fetch mycur into cache_id, registration_time, is_protected;
      exit when mycur%notfound;
      dbms_output.put_line(cache_id
          || ', ' || registration_time || ', ' || is_protected);
    end loop;
    dbms_output.put_line('<<<<<<<<<<<<<<<<<<<<');
  end;

  function t_nonidxed(tabname varchar2, colname varchar2)
      return jcc_nonidxed_array pipelined as
    type mycurtype is ref cursor;
    mycur mycurtype;
    
    s varchar(32000);
    idxtabname varchar(61);

    rid varchar(200);
    smiles varchar(4000);
  begin
    idxtabname := get_idxtabname(get_idxname(tabname, colname));
    s := 'select rowid, ' || colname
          || ' from ' || tabname
          || ' where ' || tabname || '.rowid not in (select rid from '
                || idxtabname || ')';
    dbms_output.put_line('=====> Opening cursor for : ' || s);
    open mycur for s;
    loop
      fetch mycur into rid, smiles;
      exit when mycur%notfound;
      pipe row(jcc_nonidxed_record(rid, smiles));
    end loop;
    return;
  end;

  function get_avg_cdsmiles_len(idxname varchar2) return number as
    idxtabname varchar2(61);
    avlen number;
  begin
    idxtabname := get_idxtabname(idxname);  
    execute immediate 'select avg(length(cd_smiles)) from '
              || idxtabname into avlen;
    return avlen;
  end;
end jcc_idx_stats;
/
show errors

set serveroutput on size unlimited
set line 300
column name format a80
column value format a80
-- select user from dual;
-- define table_name = 'nci_1k';
-- define column_name = 'structure';
-- select * from table(jcc_idx_stats.get_all_jcprop(jcc_idx_stats.get_idxname('&table_name', '&column_name')));  
-- call jcc_idx_stats.list_cacheids('&index_name');
quit
