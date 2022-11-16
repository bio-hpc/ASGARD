/*
 * Copyright (c) 1998-2015 ChemAxon Ltd. All Rights Reserved.
 *
 * This software is the confidential and proprietary information of
 * ChemAxon. You shall not disclose such Confidential Information
 * and shall use it only in accordance with the terms of the agreements
 * you entered into with ChemAxon.
 *
 */

whenever sqlerror exit 1

/*
 * Calibration script for JChem Oracle Cartridge. Contains methods for setting session level and permanent cost factors
 * for index and function mode of JChem substructure search to help the Oracle CBO (cost based optimizer) to select
 * better execution plans.
 *
 * author: Zsuzsanna Szabo
 * since: 31/03/2015
 */
create or replace package jchem_calibration_pkg authid current_user as

  /**
   * Calibrates the index and function mode costs of the JC_COMPARE() method for the current session.
   * Prints UPDATE statement required to be executed as JChem owner for permanent storage of the calibrated costs.
   *
   * Index mode is calibrated with a SELECT statement of the following type:
   * SELECT COUNT(*) FROM <table_name> WHERE JC_COMPARE(<column_name>, <molecule>, t:s)
   * Index mode cost is set up from the execution time of the above statement. Cost is calculated to have
   * the same time calculated in the ORACLE explain plan as the real execution time.
   *
   * Function mode cost is calibrated based on the index mode cost and the following statement:
   * SELECT COUNT(*) FROM <table_name> WHERE JC_COMPARE(<column_name>, <molecule>, t:s) and <key_name> < ?   
   *
   * For a lower limit in the part <key_name> < ? function mode is expected to perform better. For a higher
   * limit index mode should execute faster. Calibration searches for the limit where function mode and index
   * mode execution time is roughly the same and sets function mode cost to match index mode cost at this point.
   * This way below this point ORACLE will choose function mode, above this it will choose index mode execution
   * of JC_COMPARE().
   *
   * Parameters:
   * table_name   name of structural table to be used for calibration
   * column_name  structural column name
   * key_name     name of key column used for function mode calibration
   * molecule     query molecule for search performed at calibration
   */
  procedure calibrate(
        table_name in varchar2,
        column_name in varchar2,
        key_name in varchar2,
        molecule in varchar2);
 
  /**
   * Calibrates the index and function mode costs of the JC_COMPARE() method for the current session.
   * Prints UPDATE statement required to be executed as JChem owner for permanent storage of the calibrated costs.
   *
   * Given SELECT SQL statement should contain a JC_COMPARE call and a limitation on another column of the
   * same table or another table joint to the structural table. The limitation should be an integer
   * value marked with a '?' character in the statement. This value will be varied during the calibration.
   * Example:
   * SELECT COUNT(*) FROM STRUCTURE_TABLE S, ANOTHER_TABLE A WHERE
   *    JC_COMPARE(S.STRUCTURE, 'c1ccccc1', 't:s') = 1 AND
   *    A.NUMBER < ? AND
   *    S.STRUCT_KEY = A.STRUCT_KEY
   *
   * Index mode is calibrated with a SELECT statement of the following type:
   * SELECT COUNT(*) FROM <table_name> WHERE JC_COMPARE(<column_name>, <molecule>, 't:s') = 1,
   * where table name, column name and molecule are extracted from the given SELECT statement.
   * Index mode cost is set up from the execution time of the above statement. Cost is calculated to have
   * the same time calculated in the ORACLE explain plan as the real execution time.
   *
   * Function mode cost is calibrated based on the statement given as parameter and the index mode cost.
   *
   * For a lower limit value replaced into the '?' character function mode is expected to perform better.
   * For a higher limit index mode should execute faster. Calibration searches for the limit where function
   * mode and index  mode execution time is roughly the same and sets function mode cost to match index mode
   * cost at this point. This way below this point ORACLE will choose function mode, above this it will
   * choose index mode execution of JC_COMPARE().
   *
   * Parameter:
   * select_statement  SQL select statement tp calibrate for
   */
  procedure calibrate(select_statement in varchar2);

  /**
   * Returns that value for the structural key < ? below which function mode and above which
   * index mode is chosen.
   */
  function get_changing_key_limit return int;      

  /**
   * Returns the update statement that should be executed as Jchem owner to set
   * the Jchem cost factors permanently.
   */
  function get_update_statement return varchar2;

end jchem_calibration_pkg;
/
show errors;


create or replace package body jchem_calibration_pkg as

  /* The first check of constraint key < ? starts with this value.
   * It is assumed that function mode of JC_COMPARE() performs better with this limit than the index mode.
   */
  key_absolute_lower_limit constant int := 2;

  /* The first value for checking whether index mode performs better than function mode. */
  key_start_upper_limit constant int := 5000;
  check_factor constant int := 1E8;

  /* Minimal required size for table used for calibration. */
  min_table_size constant int := 10000;

  /* This CPU cost is set for function mode if already at 'key_absolute_lower_limit' index mode performs better than function mode
   * (it means that index mode is always better than function mode).
   */
  extra_large_func_cost constant int := 1E15;
 
  /* Line feed. */
  crlf constant varchar2(2) := chr(13) || chr(10);

  /* 0.1 second, times below that are considered the same. */
  time_epsilon constant int := 10;

  /* Will contain the limit at which the function mode and the index mode performs the same. */
  key_turning_limit int := -1;

  /* Update statement to be executed for permanent cost factor settings. */
  update_statement varchar2(300) := NULL;

  /* Measures the time taken to execute given SQL statement. */
  function measure_time(sql_statement in varchar2) return number as
  start_time number;
  result number;
  begin
    start_time := dbms_utility.get_time;
    execute immediate sql_statement into result;
    return dbms_utility.get_time - start_time;
  end;

  /**
   * Returns the simple not join statement for index mode calibration.
   */
  function get_simple_index_statement(
        table_name in varchar2,
        index_name in varchar2,
        column_name in varchar2,
        molecule in varchar2) return varchar2 as
  begin
    return 'select /*+ INDEX(' || table_name || ' ' || index_name || ')*/ count(*) from ' || table_name ||
            ' where jc_compare(' || column_name || ',' || ' ''' || molecule || ''',''t:s'')=1';
  end;

  /**
   * Adds word 'select' and hint for index mode execution to sql statement and
   * returns the statement with the hint added.
   */
  function prepare_index_statement(
        table_name in varchar2,
        index_name in varchar2,
        statement_skeleton varchar2) return varchar2 as
  begin
    return 'select /*+ INDEX(' || table_name || ' ' || index_name || ')*/' || statement_skeleton;
  end;
 
  /**
   * Adds word 'select' and hint for function mode execution to sql statement and
   * returns the statement with the hint added.
   */
  function prepare_no_index_statement(
        table_name in varchar2,
        index_name in varchar2,
        statement_skeleton in varchar2) return varchar2 as
  begin
    return 'select /*+ NO_INDEX(' || table_name || ' ' || index_name || ')*/' || statement_skeleton;
  end;

  /**
   * Prepares an returns statement skeleton with  JC_COMPARE(...) and key <
   */
  function prepare_statement_skeleton(
        table_name in varchar2,
        column_name in varchar2,
        molecule in varchar2,
        key in varchar2) return varchar2 as
  begin
    return ' count(*) from ' || table_name || ' where jc_compare(' || column_name || ',' ||
     ' ''' || molecule || ''',''t:s'')=1 and ' || key || ' < ?';
  end;

  /**
   * Adds given limit to statement as key < {given limit}.
   */
  function add_key_limit_to_statement(statement_skeleton in varchar2, key_limit in int) return varchar2 as
  begin
    if instr(statement_skeleton, '?') = 0 then
      raise_application_error(-20008, 'Statement should contain a ''?'' character which marks the calibrated limit.');
    end if;
    return replace(statement_skeleton, '?', to_char(key_limit));
  end;
 
  /**
   * Finds a lower and upper bound for the limit where index mode and function mode execution time is about the same.
   * At the lower bound function mode should be the winner, at the upper bound the index mode.
   */
  procedure find_key_limits(
        index_statement_skeleton in varchar2,
        function_statement_skeleton in varchar2,
        key_lower_limit out int,
        key_upper_limit out int) as
  index_time number;
  function_time number;
  begin
    key_lower_limit := key_absolute_lower_limit;
    key_upper_limit := key_start_upper_limit;
    index_time := measure_time(add_key_limit_to_statement(index_statement_skeleton, key_upper_limit));
    function_time := measure_time(add_key_limit_to_statement(function_statement_skeleton, key_upper_limit));
    while index_time - function_time > time_epsilon loop
      key_lower_limit := key_upper_limit;
      key_upper_limit := key_lower_limit * 2;
      index_time := measure_time(add_key_limit_to_statement(index_statement_skeleton, key_upper_limit));
      function_time := measure_time(add_key_limit_to_statement(function_statement_skeleton, key_upper_limit));
    end loop;
  end;
 
   /**
    * Finds the value for 'key < ?' limit where index mode and function mode execution time is about the same.
    */
  function find_turning_point(
        index_statement_skeleton in varchar2,
        function_statement_skeleton in varchar2) return int as
  key_upper_limit int;
  key_lower_limit int;
  key_current_limit int;
  index_time number;
  function_time number;
  begin
    find_key_limits(index_statement_skeleton, function_statement_skeleton, key_lower_limit, key_upper_limit);
    while key_upper_limit > key_lower_limit * 1.1 loop
      key_current_limit := (key_lower_limit + key_upper_limit) / 2;
      index_time := measure_time(add_key_limit_to_statement(index_statement_skeleton, key_current_limit));
      function_time := measure_time(add_key_limit_to_statement(function_statement_skeleton, key_current_limit));
      if function_time < index_time then
        key_lower_limit := key_current_limit;
      else
        key_upper_limit := key_current_limit;
      end if;
    end loop;
    return (key_lower_limit + key_upper_limit) / 2;
  end;

  /**
   * Executes given statement so that it becomes cached and this execution time gets stabilized.
   */
  procedure warm_up(sql_statement in varchar2) as
  result number;
  begin
    execute immediate sql_statement into result;
  end;

  /**
   * Gets row containing full, IO and CPU costs from the explain plan of the given statement.
   */
  function get_costs(statement in  varchar2) return plan_table%ROWTYPE as
  plan_table_row plan_table%ROWTYPE;
  begin
    delete from plan_table;
    execute immediate 'explain plan for ' || statement;
    select * into plan_table_row from plan_table where operation = 'SELECT STATEMENT';
    return plan_table_row;
  end;

  /**
   * Gets row containing full, IO and CPU costs from the explain plan of the given statement with
   * function mode CPU cost of JC_COMPARE() set to given value. Index mode costs remain as set before.
   */
  function get_function_costs(function_statement varchar2, function_factor int) return plan_table%ROWTYPE as
  index_cpu_cost number;
  index_io_cost number;
  index_net_cost number;
  begin
    index_cpu_cost := jchem_opti_pkg.get_cost_factor('index', 'compare', 't:s', 'cpu', NULL, NULL);
    index_io_cost := jchem_opti_pkg.get_cost_factor('index', 'compare', 't:s', 'io', NULL, NULL);
    index_net_cost := jchem_opti_pkg.get_cost_factor('index', 'compare', 't:s', 'net', NULL, NULL);
    jchem_opti_pkg.set_volatile_cost_factors('compare', index_cpu_cost, index_io_cost, index_net_cost, function_factor, 0, 0);
    return get_costs(function_statement);
  end;

  /**
   * Gets row containing full, IO and CPU costs from the explain plan of the given statement with
   * index mode CPU cost of JC_COMPARE() set to given value.
   */
  function get_index_costs(statement in varchar2, index_factor in number) return plan_table%ROWTYPE as
  begin
    jchem_opti_pkg.set_volatile_cost_factors('compare', index_factor, 0, 0, 0, 0, 0);
    return get_costs(statement);
  end;

  /**
   * Calculates the required CPU cost for JC_COMPARE() based on the given required full cost
   * and the full costs calculated by the ORACLE CBO with other CPU cost factor settings
   * (0 and a number large enough for precise calculation.
   */
  function calculate_factor_from_costs(
        required_cost in number,
        check_factor in number,
        zero_costs in plan_table%ROWTYPE,
        n_costs in plan_table%ROWTYPE) return number as        
  schema_name varchar2(100);      
  begin
    return round(check_factor *
        (1.0 * (required_cost - zero_costs.io_cost) / (n_costs.cost - zero_costs.cost) -
         1.0 * zero_costs.cpu_cost / (n_costs.cpu_cost - zero_costs.cpu_cost)));
    exception     
    when ZERO_DIVIDE then    
      select sys_context('userenv','current_schema') into schema_name from dual;
      raise_application_error(-20005, 'JChem costs can not be calibrated. Possible causes:' || crlf ||
          '- Table chosen for calibration is too small. (Choose a table that contains at least ' ||
          to_char(min_table_size * 10) || ' structures.)' || crlf ||
          '- Statistics are not gathered for the schema yet. (Please, execute ''CALL DBMS_STATS.GATHER_SCHEMA_STATS(''' || schema_name || ''')''.)' || crlf ||
          'For more details, please refer to the online Jchem Cartridge cost estimation documentation.');
  end;

  /**
   * Adjusts the CPU cost of JC_COMPARE() so that the full cost of given <fucntion_statement> will be
   * equal the full cost of the given <index_statement>.
   */
  function adjust_function_cost(
        index_statement in varchar2,
        function_statement in varchar2) return number as
  required_cost number;
  zero_costs plan_table%ROWTYPE;
  n_costs plan_table%ROWTYPE;
  begin
    delete from plan_table;
    execute immediate 'explain plan for ' || index_statement;
    select cost into required_cost from plan_table where operation = 'SELECT STATEMENT';

    zero_costs := get_function_costs(function_statement, 0);
    n_costs := get_function_costs(function_statement, check_factor);
    return calculate_factor_from_costs(required_cost, check_factor, zero_costs, n_costs);
  end;

  /**
   * Adjusts index mode CPU cost of JC_COMPARE() so that its cost in the explain plan
   * will have the given cost value.
   */
  function adjust_index_cost_to_cost(statement in varchar2, required_cost in number) return number as
  n_costs plan_table%ROWTYPE;
  zero_costs plan_table%ROWTYPE;
  begin
    zero_costs := get_index_costs(statement, 0);
    n_costs := get_index_costs(statement, check_factor);
    return calculate_factor_from_costs(required_cost, check_factor, zero_costs, n_costs);
  end;

  /**
   * Calibrates function mode cost of JC_COMPARE() so that it has the same cost as
   * index mode execution when the two modes have approximately the same execution time.
   */
  procedure calibrate_function_mode(
        alias_name in varchar2,
        index_name in varchar2,
        statement_skeleton in varchar2) as
  index_statement_skeleton varchar2(300);
  function_statement_skeleton varchar2(300);
  index_statement varchar2(300);
  function_statement varchar2(300);
  index_time number;
  function_time number;
  index_cpu_cost number;
  function_cpu_cost number;
  required_index_cost number;
  begin
    index_statement_skeleton := prepare_index_statement(alias_name, index_name, statement_skeleton);
    function_statement_skeleton := prepare_no_index_statement(alias_name, index_name, statement_skeleton);
    dbms_output.put_line('Function mode calibration index statement: ' || crlf || index_statement_skeleton);        
    dbms_output.put_line('Function mode calibration function statement: ' || crlf || function_statement_skeleton);        

    index_statement := add_key_limit_to_statement(index_statement_skeleton, key_absolute_lower_limit);
    function_statement := add_key_limit_to_statement(function_statement_skeleton, key_absolute_lower_limit);
    warm_up(index_statement);
    warm_up(function_statement);

    index_time := measure_time(index_statement);
    function_time := measure_time(function_statement);
    index_cpu_cost := jchem_opti_pkg.get_cost_factor('index', 'compare', 't:s', 'cpu', NULL, NULL);

    if index_time - function_time < time_epsilon then
      key_turning_limit := key_absolute_lower_limit;
      function_cpu_cost := extra_large_func_cost;
    else
      key_turning_limit := find_turning_point(index_statement_skeleton, function_statement_skeleton);
      index_statement := add_key_limit_to_statement(index_statement_skeleton, key_turning_limit);
      function_statement := add_key_limit_to_statement(function_statement_skeleton, key_turning_limit);
      function_cpu_cost := adjust_function_cost(index_statement, function_statement);

      -- If Oracle CBO is buggy and the value of function cost can not reach the value of
      -- index cost then raise index cost.
      if function_cpu_cost < 0 then
        function_cpu_cost := 0;
        required_index_cost := get_function_costs(function_statement, function_cpu_cost).cost;
        index_cpu_cost := adjust_index_cost_to_cost(index_statement, required_index_cost);
      end if;
    end if;    
    jchem_opti_pkg.set_volatile_cost_factors('compare', index_cpu_cost, 0, 0, function_cpu_cost, 0, 0);
  end;

  /**
   * Adjusts index mode CPU cost of JC_COMPARE() so that its time in the explain plan
   * will have the given time value.
   */
  function adjust_index_cost(statement in varchar2, exec_time in number) return number as
  n_costs plan_table%ROWTYPE;
  zero_costs plan_table%ROWTYPE;
  required_cost number;
  begin
    zero_costs := get_index_costs(statement, 0);
    n_costs := get_index_costs(statement, check_factor);
    -- dbms_utility.get_time returns time in 1/100th of a second
    required_cost := 0.01 * exec_time * n_costs.cost / n_costs.time;
    return calculate_factor_from_costs(required_cost, check_factor, zero_costs, n_costs);
  end;

  /**
   * Calibrates index mode CPU cost of JC_COMPARE() based on real execution time
   * of a simple SELECT statement with substructure search.
   */
  procedure precalibrate_index_mode(
        table_name in varchar2,
        index_name in varchar2,
        column_name in varchar2,
        molecule in varchar2) as
  statement varchar2(300);
  exec_time number;
  index_cpu_cost number;
  begin
    statement := get_simple_index_statement(table_name, index_name, column_name, molecule);
    dbms_output.put_line('Index mode calibration statement: ' || crlf || statement);        
    warm_up(statement);
    exec_time := measure_time(statement);
    index_cpu_cost := adjust_index_cost(statement, exec_time);
    jchem_opti_pkg.set_volatile_cost_factors('compare', index_cpu_cost, 0, 0, 0, 0, 0);
  end;

  /**
   * Finds JChem index name of given structural column of given table. Throws exception if not found.
   */
  function find_index_name(table_name in varchar2, column_name in varchar2) return varchar2 as
  index_name varchar2(100);
  begin
    execute immediate 'select col.index_name from user_ind_columns col, user_indexes ind where col.table_name =  ''' || table_name ||
        ''' and col.column_name = ''' || column_name || ''' and col.index_name = ind.index_name and ind.ityp_name = ''JC_IDXTYPE''' into index_name;
    return index_name;    
    exception
      when NO_DATA_FOUND then
      raise_application_error(-20001, 'No domain index of type ''JC_IDXTYPE'' is defined for ''' || table_name || '.' || column_name || '''.' || crlf ||
      'Please check table name, column name and the existence of JChem domain index on the given column.');
  end;

  /**
   * Checks whether statistics essential for cost estimation are available.
   */
  procedure check_index_statistics(index_name varchar2) as
  num_rows int;
  schema_name varchar2(100);      
  begin
    execute immediate 'select num_rows from all_tables where table_name = ''' || index_name || '_JCX''' into num_rows;
    if num_rows is NULL then
      select sys_context('userenv','current_schema') into schema_name from dual;
      raise_application_error(-20002, 'No statistics are available for Jchem index.' || crlf ||
        'Please gather statistics by executing ''CALL DBMS_STATS.GATHER_SCHEMA_STATS(''' ||
        schema_name || ''')'' before calibrating.');
    end if;
  end;

  /**
   * Checks whether statistics have been associated with JChem functions and are working properly.
   */
  procedure check_statistics_association(table_name varchar2, column_name varchar2) as
  the_cpu_cost int;
  table_size int;
  compare_result int;
  begin
    delete from plan_table;
    jchem_opti_pkg.set_volatile_cost_factors('compare', 100, 0, 0, 100000, 0, 0);
    execute immediate 'explain plan for select count(*) from ' || table_name ||
      ' where jc_compare(' || column_name || ', ''Brc1ccccc1'', ''t:s'') = 1';
    select cpu_cost into the_cpu_cost from plan_table where operation like 'DOMAIN_INDEX%';
    if the_cpu_cost is NULL then
      execute immediate 'select count(*) from ' || table_name into table_size;
      if table_size < min_table_size then
        raise_application_error(-20004, 'Table ' || table_name || ' is too small for calibration (' ||
          to_char(table_size) || ' rows).' || crlf ||
          'Try calibrating with a table that contains at least ' || to_char(min_table_size) || ' rows.');
      else
        -- Test whether connection is alive and use_password() is set.
        execute immediate 'select count(*) from ' || table_name ||
            ' where jc_compare(' || column_name || ', ''BrNSClc1ccccc1'', ''t:s'') = 1';
        raise_application_error(-20003, 'The ''assoc_stats.sql'' script has not been executed yet.' || crlf ||
          'Please run the script as the JCHEM owner before calibrating.');
      end if;
    end if;
  end;

  /**
   * Check prerequisites and perform calibration.
   */
  procedure calibrate(
        table_name in varchar2,
        alias_name in varchar2,
        column_name in varchar2,
        statement_skeleton in varchar2,
        molecule in varchar2) as
  index_cpu_cost number;
  function_cpu_cost number;
  index_name varchar2(50);
  begin
    index_name := find_index_name(table_name, column_name);
    check_index_statistics(index_name);
    check_statistics_association(table_name, column_name);

    precalibrate_index_mode(table_name, index_name, column_name, molecule);
    calibrate_function_mode(alias_name, index_name, statement_skeleton);

    index_cpu_cost := jchem_opti_pkg.get_cost_factor('index', 'compare', 't:s', 'cpu', NULL, NULL);
    function_cpu_cost := jchem_opti_pkg.get_cost_factor('func', 'compare', 't:s', 'cpu', NULL, NULL);
    dbms_output.put_line('Successfully calibrated for session. Cost factors are (index: ' || to_char(index_cpu_cost) || ', function: ' ||  to_char(function_cpu_cost) || ').');

    dbms_output.put_line('For permanent calibration execute the following statement as the JChem user:');
    update_statement := 'UPDATE JC_IDX_PROPERTY SET PROP_VALUE = '''  || to_char(index_cpu_cost) ||
        '; 0; 0; ' || to_char(function_cpu_cost) || '; 0; 0'' WHERE PROP_NAME = ''cost.factors.COMPARE.default''';
    dbms_output.put_line(update_statement);
  end;

  /**
   * Calibrates the index and function mode costs of the JC_COMPARE() method for the current session.
   * Prints UPDATE statement required to be executed as JChem owner for permanent storage of the calibrated costs.
   *
   * Index mode is calibrated with a SELECT statement of the following type:
   * SELECT COUNT(*) FROM <table_name> WHERE JC_COMPARE(<column_name>, <molecule>, 't:s') = 1
   * Index mode cost is set up from the execution time of the above statement. Cost is calculated to have
   * the same time calculated in the ORACLE explain plan as the real execution time.
   *
   * Function mode cost is calibrated based on the index mode cost and the following statement:
   * SELECT COUNT(*) FROM <table_name> WHERE JC_COMPARE(<column_name>, <molecule>, 't:s') = 1 and <key_name> < ?
   *
   * For a lower limit in the part <key_name> < ? function mode is expected to perform better. For a higher
   * limit index mode should execute faster. Calibration searches for the limit where function mode and index
   * mode execution time is roughly the same and sets function mode cost to match index mode cost at this point.
   * This way below this point ORACLE will choose function mode, above this it will choose index mode execution
   * of JC_COMPARE().
   *
   * Parameters:
   * table_name   name of structural table to be used for calibration
   * column_name  structural column name
   * key_name     name of key column used for function mode calibration
   * molecule     query molecule for search performed at calibration
   */
  procedure calibrate(
        table_name in varchar2,
        column_name in varchar2,
        key_name in varchar2,
        molecule in varchar2) as
  tab_name varchar2(50);
  col_name varchar2(50);
  k_name varchar2(50);
  statement_skeleton varchar2(300);
  begin
    tab_name := upper(table_name);
    col_name := upper(column_name);
    k_name := upper(key_name);
    statement_skeleton := prepare_statement_skeleton(tab_name, col_name, molecule, k_name);
    calibrate(tab_name, tab_name, col_name, statement_skeleton, molecule);
  end;

  /**
   * Returns the first word in the given string at ar after the given position. Returns NULL
   * if not found.
   */
  function get_next_word(
    string in varchar2,
    start_pos in int,
    after_pos out int) return varchar2 as
  word_start_pos int;
  begin
    word_start_pos := regexp_instr(string, '\w', start_pos);
    if start_pos = 0 then
      return NULL;
    else
      after_pos := regexp_instr(string, '\W', word_start_pos);
      return substr(string, word_start_pos, after_pos - word_start_pos);
    end if;
  end;

  /**
   * Checks whether the word at or after the given position (only whitespaces are allowed) is
   * an alias name or an SLQ reserved word that is allowed after a table name.
   */
  function check_for_alias_name(string in varchar2, start_pos in int) return varchar2 as
  after_pos int;
  non_whitespace_start_pos int;
  begin
    non_whitespace_start_pos := regexp_instr(string, '\S', start_pos);
    after_pos := regexp_instr(string,
      'join|full join|inner join|left join|right join|where|on|\W', non_whitespace_start_pos, 1, 0, 'i');
    if after_pos > non_whitespace_start_pos then
      return substr(string, non_whitespace_start_pos, after_pos - non_whitespace_start_pos);
    else
      return NULL;
    end if;
  end;

  /**
   * Finds name of the structural table and its alias (if present) in the given statement based on
   * found table-or-alias name or by testing each word in the statement whether it is the table name
   * if table-or-alias name is not present.
   */
  procedure find_table_and_alias_name(
    statement in varchar2,
    name_or_alias in varchar2,
    column_name in varchar2,
    table_name out varchar2,
    alias_name out varchar2) as
  pos int;
  end_pos int;
  word_pos int;
  word_after_pos int;
  table_count int;
  word varchar2(50);
  begin
    table_name := NULL;
    if name_or_alias is not NULL then
      execute immediate 'select count(*) from user_tables where table_name = ''' ||
          upper(name_or_alias) || '''' into table_count;
      if table_count = 1 then
        table_name := upper(name_or_alias);
        alias_name := table_name;
      else
        alias_name := name_or_alias;
        pos := regexp_instr(statement, '(\w+)(\W+)' || alias_name || '(\W+)', 1, 1, 0, 'i');
        end_pos := regexp_instr(statement, '\W', pos);
        table_name := upper(substr(statement, pos, end_pos - pos));
      end if;
    else
      word_pos := 1;
      word := get_next_word(statement, word_pos, word_after_pos);
      while table_name is NULL and word is not null loop
        word := upper(word);
        execute immediate 'select count(*) from user_tab_columns where table_name = ''' ||
            word || ''' and column_name = ''' || column_name || '''' into table_count;
        if table_count = 1 then
          table_name := word;
          alias_name := check_for_alias_name(statement, word_after_pos);
          if alias_name is NULL then
            alias_name := table_name;
          end if;
        else
          word_pos := word_after_pos;
          word := get_next_word(statement, word_pos, word_after_pos);
        end if;
      end loop;
      if table_name is NULL then
        raise_application_error(-20008, 'Error while parsing SELECT statement.');
      end if;
    end if;
  end;

  /**
   * Parses parameters required for index mode calibration statement and for adding hints to function
   * mode calibration statement from the given statement.
   */
  procedure parse_parameters(
        statement in varchar2,
        table_name out varchar2,
        alias_name out varchar2,
        column_name out varchar2,
        molecule out varchar2) as
  pos int;
  comma_pos int;
  dot_pos int;
  mol_end_pos int;
  name_or_alias varchar2(50) := NULL;
  begin
    pos := regexp_instr(statement, 'jc_compare', 1, 1, 1, 'i');
    if pos = 0 then
      raise_application_error(-20007, 'Statement should contain a call to function JC_COMPARE().');
    end if;
    pos := instr(statement, '(', pos);
    comma_pos := instr(statement, ',', pos);
    column_name := regexp_replace(substr(statement, pos + 1, comma_pos - pos - 1), '\s');
    dot_pos := instr(column_name, '.');
    if dot_pos <> 0 then
      name_or_alias := substr(column_name, 1, dot_pos - 1);
      column_name := substr(column_name, dot_pos + 1);
    end if;
    column_name := upper(column_name);

    pos := instr(statement, '''', comma_pos + 1);
    mol_end_pos := instr(statement, '''', pos + 1);
    molecule := substr(statement, pos + 1, mol_end_pos - pos - 1);

    find_table_and_alias_name(statement, name_or_alias, column_name, table_name, alias_name);
  end;

  /**
   * Calibrates the index and function mode costs of the JC_COMPARE() method for the current session.
   * Prints UPDATE statement required to be executed as JChem owner for permanent storage of the calibrated costs.
   *
   * Given SELECT SQL statement should contain a JC_COMPARE call of type
   * JC_COMPARE(<column_name>, <molecule>, 't:s') = 1 and a limitation on another column of the
   * same table or another table joint to the structural table. The limitation should be an integer
   * value marked with a '?' character in the statement. This value will be varied during the calibration.
   * Example:
   * SELECT COUNT(*) FROM STRUCTURE_TABLE S, ANOTHER_TABLE A WHERE
   *    JC_COMPARE(S.STRUCTURE, 'c1ccccc1', 't:s') = 1 AND
   *    A.NUMBER < ? AND
   *    S.STRUCT_KEY = A.STRUCT_KEY
   *
   * Index mode is calibrated with a SELECT statement of the following type:
   * SELECT COUNT(*) FROM <table_name> WHERE JC_COMPARE(<column_name>, <molecule>, 't:s') = 1,
   * where table name, column name and molecule are extracted from the given SELECT statement.
   * Index mode cost is set up from the execution time of the above statement. Cost is calculated to have
   * the same time calculated in the ORACLE explain plan as the real execution time.
   *
   * Function mode cost is calibrated based on the statement given as parameter and the index mode cost.
   *
   * For a lower limit value replaced into the '?' character function mode is expected to perform better.
   * For a higher limit index mode should execute faster. Calibration searches for the limit where function
   * mode and index  mode execution time is roughly the same and sets function mode cost to match index mode
   * cost at this point. This way below this point ORACLE will choose function mode, above this it will
   * choose index mode execution of JC_COMPARE().
   *
   * Parameter:
   * select_statement  SQL select statement tp calibrate for
   */
  procedure calibrate(select_statement in varchar2) as
  table_name varchar2(50);
  alias_name varchar2(50);
  column_name varchar2(50);
  molecule varchar2(50);
  statement_skeleton varchar2(300);
  begin
    if upper(substr(select_statement, 1, 6)) <> 'SELECT' then
      raise_application_error(-20006, 'Statement should begin with the word SELECT.');
    end if;
    -- Check correctness of the statement.
    execute immediate add_key_limit_to_statement(select_statement, 1);
    statement_skeleton := substr(select_statement, 7);
    parse_parameters(statement_skeleton, table_name, alias_name, column_name, molecule);
    calibrate(table_name, alias_name, column_name, statement_skeleton, molecule);
  end;

  /**
   * Returns that value for the structural key < ? below which function mode and above which
   * index mode is chosen.
   */
  function get_changing_key_limit return int as
  begin
    return key_turning_limit;
  end;  

  /**
   * Returns the update statement that should be executed as Jchem owner to set
   * the Jchem cost factors permanently.
   */
  function get_update_statement return varchar2 as
  begin
    return update_statement;
  end;  

end;
/
show errors;

quit




