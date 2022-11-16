Rem
Rem Copyright (c) 1998-2013 ChemAxon Ltd. All Rights Reserved.
Rem This software is the confidential and proprietary information of
Rem ChemAxon. You shall not disclose such Confidential Information
Rem and shall use it only in accordance with the terms of the agreements
Rem you entered into with ChemAxon.
Rem

whenever sqlerror exit 1

---------------------------------------------------------------------
--    			JChem Data Cartridge			   --
---------------------------------------------------------------------

-- Type used for storing the results of screenning
--CREATE OR REPLACE TYPE SMI_ARRAY is VARRAY(1000000000) OF VARCHAR2(32767)
--/

CREATE OR REPLACE TYPE CD_ID_ARRAY is VARRAY(1000000000) OF INTEGER;
/

-- Type used for storing the results of sss
CREATE OR REPLACE TYPE RESARRAY is VARRAY(1000000000) OF VARCHAR2(256)
/

-- Type used for storing the result rowids
CREATE OR REPLACE TYPE RIDARRAY is VARRAY(1000000000) OF VARCHAR2(20)
/

CREATE OR REPLACE TYPE MOLPROPS_ARRAY is VARRAY(1000000000) OF VARCHAR(2000)
/

CREATE OR REPLACE TYPE MOLPROPS_ARRAY_ARRAY is VARRAY(1000000000) OF MOLPROPS_ARRAY
/
show errors;

CREATE OR REPLACE TYPE CHAR_ARRAY is VARRAY(2000000000) OF VARCHAR2(32767)
/
show errors;

CREATE OR REPLACE TYPE COMPOSITE_CHAR as OBJECT (c VARCHAR2(32767));
/
show errors;

CREATE OR REPLACE TYPE COMP_CHAR_ARRAY is VARRAY(20000000) OF COMPOSITE_CHAR;
/
show errors;

CREATE OR REPLACE TYPE CLOB_ARRAY is VARRAY(2000000000) OF CLOB
/
show errors;

CREATE OR REPLACE TYPE COMPOSITE_CLOB as OBJECT (b CLOB);
/
show errors;

CREATE OR REPLACE TYPE COMP_CLOB_ARRAY is varray(20000000) OF COMPOSITE_CLOB;
/
show errors;

CREATE OR REPLACE TYPE BLOB_ARRAY is VARRAY(2000000000) OF BLOB
/
show errors;

CREATE OR REPLACE TYPE COMPOSITE_BLOB as OBJECT (b BLOB);
/
show errors;

CREATE OR REPLACE TYPE COMP_BLOB_ARRAY is varray(20000000) OF COMPOSITE_BLOB;
/
show errors;

CREATE OR REPLACE TYPE CHAR_PRODUCT_RECORD AS OBJECT (product VARCHAR2(32767), synthesis_code VARCHAR2(32767));
/
show errors;

CREATE OR REPLACE TYPE CHAR_PRODUCT_ARRAY IS VARRAY(20000000) OF CHAR_PRODUCT_RECORD;
/
show errors;

CREATE OR REPLACE TYPE CLOB_PRODUCT_RECORD AS OBJECT (product CLOB, synthesis_code VARCHAR2(32767));
/
show errors;

CREATE OR REPLACE TYPE CLOB_PRODUCT_ARRAY IS VARRAY(20000000) OF CLOB_PRODUCT_RECORD;
/
show errors;

CREATE OR REPLACE TYPE BLOB_PRODUCT_RECORD AS OBJECT (product BLOB, synthesis_code VARCHAR2(32767));
/
show errors;

CREATE OR REPLACE TYPE BLOB_PRODUCT_ARRAY IS VARRAY(20000000) OF BLOB_PRODUCT_RECORD;
/
show errors;

create or replace type SCAN_ARGUMENTS is varray(6) of varchar2(32767);
/
show errors;

create or replace type SCAN_CONTEXT as object (
  scanId number,
  baseTableSchema varchar(30),
  baseTableName varchar(30),
  indexedColumnName varchar(30),
  indexSchema varchar(30),
  indexName varchar(30),
  usrOpId varchar2(200),
  startDate date,
  operatorName varchar(30),
  arguments SCAN_ARGUMENTS,
  note varchar2(32767)
)
/
show errors;

create or replace type SCAN_CONTEXT_ARRAY is varray(20000000) of SCAN_CONTEXT;
/
show errors;

create or replace type ERROR_RECORD as object(
  scanId number,
  rid varchar2(100),
  errorMessage varchar2(32767),
  note varchar2(32767)
)
/
show errors;

create or replace type ERROR_RECORD_ARRAY is varray(20000000) of ERROR_RECORD;
/
show errors;

create or replace type search_exec_info as object(
  screened_count_unique integer,
  screened_count_total integer
)
/
show errors;

create or replace type gmem_util_record as object (description varchar2(300), util_mb number);
/
show errors;

create or replace type gmem_util_array is VARRAY(20000000) of gmem_util_record;
/
show errors;

create or replace type taskinfo_record as object (
  id number,
  user_name varchar2(100),
  description varchar2(4000),
  estim_memory_use number,
  timeout number,
  jcc_session_id varchar2(255),
  user_assigned_op_id varchar2(255),
  start_time timestamp,
  last_rescheduled timestamp 
  );
/
show errors;

create or replace type taskinfo_array is varray(20000000) of taskinfo_record;
/
show errors;

create or replace type license_record as object (
  software varchar2(4000),
  license_term varchar2(4000),
  licensee varchar2(4000),
  licensed_unit_count varchar2(4000),
  expiration_date varchar2(4000),
  support_expiration_date varchar2(4000),
  license_comment varchar2(4000),
  restrictions VARCHAR2(4000)
  );
/
show errors;

create or replace type license_table is table of license_record;
/
show errors;

create or replace package jchem_refcur_pkg is
  type refcur_t is ref cursor;
end jchem_refcur_pkg;
/
show errors;

-- CREATE INDEXTYPE IMPLEMENTATION TYPE
CREATE OR REPLACE TYPE jc_idxtype_im AUTHID CURRENT_USER AS OBJECT
(
  ia sys.odciindexinfo,
  scan_num number,
  cnum number,

  STATIC FUNCTION ODCIGetInterfaces(ifclist OUT sys.ODCIObjectList) 
								RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE,
  STATIC FUNCTION ODCIIndexCreate(ia sys.odciindexinfo, parms VARCHAR2, 
						env sys.ODCIEnv) RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE,
  STATIC FUNCTION ODCIIndexAlter(ia SYS.ODCIIndexInfo, parms IN OUT VARCHAR2, 
				 alter_option NUMBER, env sys.ODCIEnv)
      RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE,

  STATIC FUNCTION ODCIIndexDrop(ia sys.odciindexinfo, 
				env sys.ODCIEnv) RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE,

  STATIC FUNCTION ODCIIndexTruncate(ia sys.ODCIIndexInfo,
                                env sys.ODCIEnv) RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE,
                                
  STATIC FUNCTION ODCIIndexStart(sctx in out jc_idxtype_im, 
			 ia sys.odciindexinfo,
                         op sys.odciPredInfo, 
                         qi sys.ODCIQueryInfo, 
                         strt NUMBER, 
                         stop NUMBER,
                         env sys.ODCIEnv) RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE,

  STATIC FUNCTION ODCIIndexStart(sctx in out NOCOPY jc_idxtype_im, 
			 ia sys.odciindexinfo,
                         op sys.odciPredInfo, 
                         qi sys.ODCIQueryInfo, 
                         strt NUMBER, 
                         stop NUMBER,
                         query CLOB, env sys.ODCIEnv)
            RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE,

  STATIC FUNCTION ODCIIndexStart(sctx in out jc_idxtype_im, 
			 ia sys.odciindexinfo,
                         op sys.odciPredInfo, 
                         qi sys.ODCIQueryInfo, 
                         strt NUMBER, 
                         stop NUMBER,
                         query CLOB, 
			 param VARCHAR2, env sys.ODCIEnv)
      RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE,

  STATIC FUNCTION ODCIIndexStart(sctx in out jc_idxtype_im, 
			 ia sys.odciindexinfo,
                         op sys.odciPredInfo, 
                         qi sys.ODCIQueryInfo, 
                         strt NUMBER, 
                         stop NUMBER,
                         query CLOB, 
			 param1 number,
                         param2 number,
                         env sys.ODCIEnv)
      RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE,

  STATIC FUNCTION ODCIIndexStart(sctx in out NOCOPY jc_idxtype_im, 
        ia sys.odciindexinfo,
        op sys.odciPredInfo, 
        qi sys.ODCIQueryInfo, 
        strt NUMBER, 
        stop NUMBER,
        query BLOB, env sys.ODCIEnv) RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE,

  STATIC FUNCTION ODCIIndexStart(sctx in out jc_idxtype_im, 
        ia sys.odciindexinfo,
        op sys.odciPredInfo, 
        qi sys.ODCIQueryInfo, 
        strt NUMBER, 
        stop NUMBER,
        query BLOB, 
        param VARCHAR2, env sys.ODCIEnv) RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE,

  STATIC FUNCTION ODCIIndexStart(sctx in out jc_idxtype_im, 
        ia sys.odciindexinfo,
        op sys.odciPredInfo, 
        qi sys.ODCIQueryInfo, 
        strt NUMBER, 
        stop NUMBER,
        query BLOB, 
        param1 number,
        param2 number,
        env sys.ODCIEnv) RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE,

  MEMBER FUNCTION ODCIIndexFetch(self IN OUT NOCOPY jc_idxtype_im, 
	nrows NUMBER, rids OUT sys.odciridlist, env sys.ODCIEnv)
    RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE,

  MEMBER FUNCTION ODCIIndexClose(self IN OUT jc_idxtype_im, 
						env sys.ODCIEnv) RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE,

  STATIC FUNCTION ODCIIndexInsert(ia sys.odciindexinfo, rid VARCHAR2, 
	               	       newval VARCHAR2, env sys.ODCIEnv)
      RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE,
  
  STATIC FUNCTION ODCIIndexDelete(ia sys.odciindexinfo, rid VARCHAR2, 
            		       oldval VARCHAR2, env sys.ODCIEnv)
      RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE,
  
  STATIC FUNCTION ODCIIndexUpdate(ia sys.odciindexinfo, rid VARCHAR2, 
       	      oldval VARCHAR2, newval VARCHAR2, env sys.ODCIEnv)
      RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE,

  STATIC FUNCTION ODCIIndexInsert(ia sys.odciindexinfo, rid VARCHAR2, 
             newval CLOB, env sys.ODCIEnv)
      RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE,
          
  STATIC FUNCTION ODCIIndexDelete(ia sys.odciindexinfo, rid VARCHAR2, 
             oldval CLOB, env sys.ODCIEnv)
      RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE,

  STATIC FUNCTION ODCIIndexUpdate(ia sys.odciindexinfo, rid VARCHAR2, 
       	      oldval CLOB, newval CLOB, env sys.ODCIEnv)
      RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE,

  STATIC FUNCTION ODCIIndexInsert(ia sys.odciindexinfo, rid VARCHAR2, 
        newval BLOB, env sys.ODCIEnv) RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE,
  
  STATIC FUNCTION ODCIIndexDelete(ia sys.odciindexinfo, rid VARCHAR2, 
        oldval BLOB, env sys.ODCIEnv) RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE,
  
  STATIC FUNCTION ODCIIndexUpdate(ia sys.odciindexinfo, rid VARCHAR2,
        oldval BLOB, newval BLOB, env sys.ODCIEnv) RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE,

  STATIC FUNCTION ODCIIndexExchangePartition(ia sys.ODCIIndexInfo, 
        ia1 sys.ODCIIndexInfo, env sys.ODCIEnv) RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE,
  STATIC FUNCTION ODCIIndexMergePartition(ia sys.ODCIIndexInfo, 
        part_name1 sys.ODCIPartInfo, part_name2 sys.ODCIPartInfo, 
        parms VARCHAR2, env sys.ODCIEnv) RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE,
  STATIC FUNCTION ODCIIndexSplitPartition(ia sys.ODCIIndexInfo, 
        part_name1 sys.ODCIPartInfo, part_name2 sys.ODCIPartInfo, 
        parms VARCHAR2, env sys.ODCIEnv) RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE
);
/
show errors


-------------------------------------------------------------------
-- Create package used by ODCIIndexGetMetadata and jc_idxtype_im 
-------------------------------------------------------------------
CREATE OR REPLACE PACKAGE jchem_core_pkg AUTHID CURRENT_USER AS

  lastError varchar2(32767);

  simCalcSeparator varchar2(1) := ';';
  
  useArrays boolean:= TRUE;
  
  procedure set_use_arrays(use_arrays varchar2);
  
  function test(password varchar2) return varchar2;

  function trim_error_messages return number as language java
      name 'chemaxon.jchem.cartridge.oraresident.JFunctions.trimErrorMessages() return int';
  
  procedure handle_java_error(errm varchar2);

  FUNCTION get_cartowner_schema RETURN VARCHAR2 DETERMINISTIC PARALLEL_ENABLE;

  FUNCTION getJChemVersion RETURN VARCHAR2 DETERMINISTIC PARALLEL_ENABLE;

  FUNCTION getTableVersion RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;

  PROCEDURE init PARALLEL_ENABLE AS LANGUAGE JAVA NAME
      'chemaxon.jchem.cartridge.oraresident.dbsession.JavaStoredProcSession.instance()';

  PROCEDURE set_jccowner_schema(jcc_owner_schema varchar2)
      PARALLEL_ENABLE AS LANGUAGE JAVA NAME
      'chemaxon.jchem.cartridge.oraresident.dbsession.JavaStoredProcSession.instance(java.lang.String)';

  PROCEDURE check_master_table PARALLEL_ENABLE AS LANGUAGE JAVA NAME
      'chemaxon.jchem.cartridge.oraresident.JFunctions.checkMasterTable()';

  PROCEDURE set_master_property(idx_schema varchar2, prop_name varchar2, prop_value varchar2)
      PARALLEL_ENABLE AS LANGUAGE JAVA NAME
      'chemaxon.jchem.cartridge.oraresident.JFunctions.setMasterProperty(
                                                      java.lang.String,
                                                      java.lang.String,
                                                      java.lang.String)';

  PROCEDURE upgr_from_pre24 AS LANGUAGE JAVA NAME
        'chemaxon.jchem.cartridge.oraresident.JFunctions.upgradeAllFromPre_2_4()';

  FUNCTION get_remote_environment(result out varchar2)
      RETURN VARCHAR2 DETERMINISTIC PARALLEL_ENABLE AS LANGUAGE JAVA NAME
	'chemaxon.jchem.cartridge.oraresident.JFunctions.getEnvironment(java.lang.String[])
            return java.lang.String';

  FUNCTION getEnvironment RETURN VARCHAR2 DETERMINISTIC PARALLEL_ENABLE;
  
  FUNCTION getLicenses RETURN license_table DETERMINISTIC PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
    'chemaxon.jchem.cartridge.oraresident.LicenseType.getLicenses()
            return oracle.sql.ARRAY';

  FUNCTION jctf_autocalccts(index_schema VARCHAR2,
                            index_name VARCHAR2)
  RETURN COMP_CHAR_ARRAY PIPELINED;
  FUNCTION jctf_autocalccts(index_schema VARCHAR2,
                            index_name VARCHAR2,
                            index_part VARCHAR2)
  RETURN COMP_CHAR_ARRAY PIPELINED;
  FUNCTION jctf_autocalccts_bycol(table_schema VARCHAR2,
                                  table_name VARCHAR2,
                                  table_col VARCHAR2)
  RETURN COMP_CHAR_ARRAY PIPELINED;

  function get_jcc_session_id return varchar2 as language java name
      'chemaxon.jchem.cartridge.oraresident.dbsession.JavaStoredProcSession.jccSessionId() return java.lang.String';

  function use_password0(login varchar2, password VARCHAR2) return varchar2
      PARALLEL_ENABLE AS LANGUAGE JAVA NAME
        'chemaxon.jchem.cartridge.oraresident.JFunctions.usePassword(
          java.lang.String, java.lang.String)
          return java.lang.String';

  procedure use_password(password varchar2) parallel_enable;

  procedure use_password(login varchar2, password varchar2) parallel_enable;

  function use_default_account0 return varchar2
      PARALLEL_ENABLE AS LANGUAGE JAVA NAME
        'chemaxon.jchem.cartridge.oraresident.JFunctions.useDefaultAccount()
          return java.lang.String';

  procedure use_default_account;

  function set_password(password varchar2) return number parallel_enable;

  PROCEDURE checkTableVersion(
        indexSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2,
        tblSchema VARCHAR2, tblName VARCHAR2, colName VARCHAR2, colType VARCHAR2) PARALLEL_ENABLE 
        AS LANGUAGE JAVA NAME
        'chemaxon.jchem.cartridge.oraresident.JFunctions.checkTableVersion(
        java.lang.String, java.lang.String, java.lang.String,
        java.lang.String, java.lang.String, java.lang.String, java.lang.String)';

  PROCEDURE checkTableVersionEx(func_name VARCHAR2,
                                indexctx IN SYS.ODCIIndexCtx) PARALLEL_ENABLE;

  PROCEDURE regenIndexes(indexSchema VARCHAR2, idxName VARCHAR2,
              qualifTblName VARCHAR2) PARALLEL_ENABLE
    AS LANGUAGE JAVA NAME 'chemaxon.jchem.cartridge.oraresident.JFunctions.checkTableVersion(
          java.lang.String, java.lang.String, java.lang.String)';

  FUNCTION indexCreate(indexSchema VARCHAR2, indexName VARCHAR2,
    indexPartition VARCHAR2, tblSchema VARCHAR2, tblName VARCHAR2,tblPartition VARCHAR2,
    colName VARCHAR2, colType VARCHAR2, parms VARCHAR2) RETURN VARCHAR2
    PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
	'chemaxon.jchem.cartridge.oraresident.Indexing.indexCreate(java.lang.String, 
		java.lang.String, java.lang.String, java.lang.String, java.lang.String,
		java.lang.String, java.lang.String, java.lang.String, java.lang.String)
            return java.lang.String';
			
  FUNCTION indexAlter(indexSchema VARCHAR2, idxName VARCHAR2,
                      indexPartition VARCHAR2, tblSchema VARCHAR2,
                      tblName VARCHAR2, tablePartition VARCHAR2,
                      colName VARCHAR2, colType VARCHAR2, parms VARCHAR2, alter_option NUMBER)
    RETURN NUMBER PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
	'chemaxon.jchem.cartridge.oraresident.Indexing.indexAlter(java.lang.String, 
        java.lang.String, java.lang.String, java.lang.String, java.lang.String,
        java.lang.String, java.lang.String, java.lang.String, java.lang.String, int) return int';

  FUNCTION indexDrop(schemaName VARCHAR2, idxName VARCHAR2,
        idxPartition VARCHAR2, tblSchema VARCHAR2, tblName VARCHAR2)
    RETURN NUMBER PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
	'chemaxon.jchem.cartridge.oraresident.Indexing.indexDrop(java.lang.String,
        java.lang.String, java.lang.String, java.lang.String,
        java.lang.String) return int';

  FUNCTION indexTruncate(schemaName VARCHAR2, idxName VARCHAR2,
        idxPartition VARCHAR2, tblSchema VARCHAR2, tblName VARCHAR2)
    RETURN NUMBER PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
	'chemaxon.jchem.cartridge.oraresident.Indexing.indexTruncate(java.lang.String,
        java.lang.String, java.lang.String, java.lang.String,
        java.lang.String) return int';

  function get_idxtable_qname(idxSchema VARCHAR2,
                              idxName VARCHAR2,
                              partition VARCHAR2,
                              jctable number) return varchar2
    as language java name
    'chemaxon.jchem.cartridge.oraresident.JFunctions.getIdxTableQName(java.lang.String,
                                    java.lang.String,java.lang.String, int)
     return java.lang.String';

  FUNCTION exchangePartitions(localIdxSchema VARCHAR2,
        localIdxName VARCHAR2, localIdxPartition VARCHAR2,
        globalIdxSchema VARCHAR2, globalIdxName VARCHAR2)
    RETURN NUMBER PARALLEL_ENABLE AS LANGUAGE JAVA NAME
	'chemaxon.jchem.cartridge.oraresident.JFunctions.exchangePartitions(
        java.lang.String, java.lang.String, java.lang.String, 
        java.lang.String, java.lang.String) return int';

  FUNCTION evaluate_arr(target VARCHAR2, options VARCHAR2, rid VARCHAR2,
        idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2,
        tblSchema VARCHAR2, tblName VARCHAR2)
    RETURN CHAR_ARRAY PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
    'chemaxon.jchem.cartridge.oraresident.JFunctions.evaluateArr(java.lang.String,
      java.lang.String, java.lang.String, java.lang.String,
      java.lang.String, java.lang.String, java.lang.String,
      java.lang.String) return java.lang.String[]';

  FUNCTION autocalccts(index_schema VARCHAR2,
                       index_name VARCHAR2,
                       index_partition VARCHAR2)
  RETURN CHAR_ARRAY PARALLEL_ENABLE AS LANGUAGE JAVA NAME
    'chemaxon.jchem.cartridge.oraresident.JFunctions.getAutoCalcCtsARR(
    java.lang.String, java.lang.String, java.lang.String) return
    java.lang.String[]';

  FUNCTION autocalccts_bycol(table_schema VARCHAR2,
                             table_name VARCHAR2,
                             col_name VARCHAR2)
  RETURN CHAR_ARRAY PARALLEL_ENABLE AS LANGUAGE JAVA NAME
    'chemaxon.jchem.cartridge.oraresident.JFunctions.getAutoCalcCtsByColARR(
    java.lang.String, java.lang.String, java.lang.String) return
    java.lang.String[]';

  FUNCTION index_scan(indexSchema VARCHAR2, idxName VARCHAR2,
                       indexPartition VARCHAR2, tblSchema VARCHAR2, tblname VARCHAR2, 
					   colName VARCHAR2, colType VARCHAR2, optype VARCHAR2,
                       opflavor VARCHAR2, strt NUMBER, stop NUMBER, opFlag NUMBER, 
					   query VARCHAR2, options VARCHAR2, scanId number)
              RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME
	'chemaxon.jchem.cartridge.oraresident.JFunctions.indexScan(java.lang.String, java.lang.String, 
			java.lang.String, java.lang.String, java.lang.String,
			java.lang.String, java.lang.String, java.lang.String, 
                        java.lang.String, java.lang.Double, java.lang.Double, int, 
						java.lang.String, java.lang.String, long)
                        return java.lang.String'; 

  function get_search_info(usrOpId varchar2) return search_exec_info parallel_enable
          as language java name
        'chemaxon.jchem.cartridge.oraresident.resultset.SearchExInfoSupport.getSearchExInfo(
            java.lang.String) return oracle.sql.STRUCT';

  FUNCTION get_hit_count(tblSchema VARCHAR2, tblname VARCHAR2, colName VARCHAR2, 
                         query VARCHAR2, options VARCHAR2)
        RETURN NUMBER PARALLEL_ENABLE AS LANGUAGE JAVA NAME
	'chemaxon.jchem.cartridge.oraresident.JFunctions.getHitCount(java.lang.String, 
			java.lang.String, java.lang.String, java.lang.String,  
                        java.lang.String) return int'; 

  FUNCTION similaritySearch(query VARCHAR2, stype VARCHAR2,
        strt NUMBER, stop NUMBER, opFlag number, searchOptions VARCHAR2,
        indexSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2,
        tblSchema VARCHAR2, tblname VARCHAR2, colName VARCHAR2, colType VARCHAR2, scanId number)
  RETURN VARCHAR2
      PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
	'chemaxon.jchem.cartridge.oraresident.JcSimilarity.getSimilarity(
                java.lang.String, java.lang.String,
				java.lang.Double, java.lang.Double, int,
                java.lang.String, java.lang.String, java.lang.String,
                java.lang.String, java.lang.String, java.lang.String,
                java.lang.String, java.lang.String, long)
          return java.lang.String';
  
  FUNCTION get_rowids_array(scanId NUMBER, nrows NUMBER,
            result_set out RIDARRAY)
      return VARCHAR2 PARALLEL_ENABLE is language java name
      'chemaxon.jchem.cartridge.oraresident.JFunctions.getRowidsArray(long, int, oracle.sql.ARRAY[])
      return java.lang.String';

  FUNCTION get_rowids_resultset(scanId NUMBER, nrows NUMBER,
            result_set out jchem_refcur_pkg.refcur_t)
      return VARCHAR2 PARALLEL_ENABLE is language java name
      'chemaxon.jchem.cartridge.oraresident.JFunctions.getRowidsResultSet(long, int, java.sql.ResultSet[])
      return java.lang.String';

  FUNCTION nr_remaining_rowids(scanId NUMBER, nr out NUMBER) return VARCHAR2
      PARALLEL_ENABLE is language java name
      'chemaxon.jchem.cartridge.oraresident.JFunctions.nrRemainingRowids(long, long[])
          return java.lang.String';

  FUNCTION close_scan_resultset(scanId NUMBER) RETURN VARCHAR2
      PARALLEL_ENABLE is language java name
      'chemaxon.jchem.cartridge.oraresident.JFunctions.closeScanResultSet(long) return java.lang.String';

  FUNCTION exec_function(sqlOperator VARCHAR2, target VARCHAR2,
              query VARCHar2, options VARCHAR2, rid VARCHAR2,
              idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2,
              tblSchema VARCHAR2, tblname VARCHAR2, colName VARCHAR2, colType VARCHAR2,
              scanId number, scanFlg number, result out varchar2)
    RETURN VARCHAR2 PARALLEL_ENABLE
	AS LANGUAGE JAVA NAME 'chemaxon.jchem.cartridge.oraresident.JFunctions.execFunction(
		java.lang.String, java.lang.String, java.lang.String,
		java.lang.String, java.lang.String, java.lang.String, 
                java.lang.String, java.lang.String, java.lang.String, java.lang.String,
                java.lang.String, java.lang.String, int, int, java.lang.String[])
                return java.lang.String';

  FUNCTION insert_mol_into_idxtable(str VARCHAR2, indexSchema VARCHAR2,
                idxName VARCHAR2, idxPartition VARCHAR2,
                tblSchema VARCHAR2, tblname VARCHAR2, colName VARCHAR2, colType VARCHAR2,
                rid VARCHAR2) RETURN VARCHAR2
        PARALLEL_ENABLE AS LANGUAGE JAVA NAME
		'chemaxon.jchem.cartridge.oraresident.JCartDml.insertMolIntoIndexTable(
		java.lang.String, java.lang.String, java.lang.String, java.lang.String,
		java.lang.String, java.lang.String, java.lang.String, java.lang.String, java.lang.String)
                return java.lang.String';

  PROCEDURE delete_mol_from_idxtable(idxSchema VARCHAR2, idxName VARCHAR2,
		idxPartition VARCHAR2, tblSchema VARCHAR2, tblname VARCHAR2, 
		rid VARCHAR2) PARALLEL_ENABLE AS LANGUAGE JAVA NAME
		'chemaxon.jchem.cartridge.oraresident.JCartDml.deleteMolFromIndexTable(java.lang.String, 
		java.lang.String, java.lang.String, java.lang.String,
		java.lang.String, java.lang.String)';

  FUNCTION update_mol_idxtable(oldval VARCHAR2, str VARCHAR2,
                indexSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2,
                tblSchema VARCHAR2, tblname VARCHAR2, colName VARCHAR2, colType VARCHAR2,
                rid VARCHAR2)
        RETURN VARCHAR2
        PARALLEL_ENABLE AS LANGUAGE JAVA NAME
		'chemaxon.jchem.cartridge.oraresident.JCartDml.updateMolIndexTable(java.lang.String,
		java.lang.String, java.lang.String, java.lang.String, java.lang.String,
		java.lang.String, java.lang.String, java.lang.String, java.lang.String, java.lang.String)
                return java.lang.String';

  FUNCTION calc_molProp(query VARCHAR2, propName VARCHAR2, rid VARCHAR2, 
		idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2,
        tblSchema VARCHAR2, tableName VARCHAR2)
        RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
		'chemaxon.jchem.cartridge.oraresident.JFunctions.calcMolProp(
                java.lang.String, java.lang.String, java.lang.String,
		java.lang.String, java.lang.String, java.lang.String,
		java.lang.String, java.lang.String) return java.lang.String';

  FUNCTION calc_molPropNum(query VARCHAR2, propName VARCHAR2, rid VARCHAR2, 
		idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2,
        tblSchema VARCHAR2, tableName VARCHAR2)
        RETURN NUMBER PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
		'chemaxon.jchem.cartridge.oraresident.JFunctions.calcMolPropNum(
                java.lang.String, java.lang.String, java.lang.String,
		java.lang.String, java.lang.String, java.lang.String,
		java.lang.String, java.lang.String) return java.lang.Double';

  FUNCTION calc_molProp_from_idx(rid VARCHAR2, propName VARCHAR2, 
	idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2,
    tblSchema VARCHAR2, tableName VARCHAR2)
    RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
	'chemaxon.jchem.cartridge.oraresident.JFunctions.calcMolPropFromRowid(
    java.lang.String, java.lang.String, java.lang.String, java.lang.String,
    java.lang.String, java.lang.String, java.lang.String)
    return java.lang.String';

  FUNCTION calc_molPropNum_from_idx(rid VARCHAR2, propName VARCHAR2, 
	idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2,
    tblSchema VARCHAR2, tableName VARCHAR2)
    RETURN NUMBER PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
	'chemaxon.jchem.cartridge.oraresident.JFunctions.calcMolPropNumFromRowid(
    java.lang.String, java.lang.String, java.lang.String, java.lang.String,
    java.lang.String, java.lang.String, java.lang.String)
    return java.lang.Double';

  FUNCTION molconvertv(query VARCHAR2, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2, rid VARCHAR2,
        idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2,
        tblSchema VARCHAR2, tableName VARCHAR2)
        RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME
        'chemaxon.jchem.cartridge.oraresident.JFunctions.molConvertVV(java.lang.String,
        java.lang.String, java.lang.String, java.lang.String, java.lang.String,
        java.lang.String, java.lang.String, java.lang.String, java.lang.String,
        java.lang.String)
        return java.lang.String';

  FUNCTION molconvertv(query CLOB, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2, rid VARCHAR2,
        idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2,
        tblSchema VARCHAR2, tableName VARCHAR2)
        RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME
        'chemaxon.jchem.cartridge.oraresident.JFunctions.molConvertVC(oracle.sql.CLOB,
        java.lang.String, java.lang.String, java.lang.String, java.lang.String,
        java.lang.String, java.lang.String, java.lang.String, java.lang.String,
        java.lang.String)
        return java.lang.String';

  FUNCTION molconvertv(query BLOB, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2, rid VARCHAR2,
        idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2,
        tblSchema VARCHAR2, tableName VARCHAR2)
        RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME
        'chemaxon.jchem.cartridge.oraresident.JFunctions.molConvertVB(oracle.sql.BLOB,
        java.lang.String, java.lang.String, java.lang.String, java.lang.String,
        java.lang.String, java.lang.String, java.lang.String, java.lang.String,
        java.lang.String)
        return java.lang.String';

  FUNCTION molconvertv_from_idx(rid VARCHAR2, type VARCHAR2,
        idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2,
        tblSchema VARCHAR2, tableName VARCHAR2)
        RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME
        'chemaxon.jchem.cartridge.oraresident.JFunctions.molConvertFromRowId(
        java.lang.String, java.lang.String, java.lang.String,
        java.lang.String, java.lang.String, java.lang.String)
        return java.lang.String';
  
  FUNCTION send_user_func(name VARCHAR2, delim VARCHAR2, 
	params VARCHAR2) RETURN VARCHAR2 AS LANGUAGE JAVA NAME 
		'chemaxon.jchem.cartridge.oraresident.JFunctions.sendUserFunc(java.lang.String, 
		java.lang.String, java.lang.String) return java.lang.String';

  PROCEDURE send_user_func_batch(opName VARCHAR2, fieldName VARCHAR2, 
		params VARCHAR2, indexSchema VARCHAR2, tableSchema VARCHAR2,
		tableName VARCHAR2, idxSchema VARCHAR2, scanId NUMBER)
      AS LANGUAGE JAVA NAME
	'chemaxon.jchem.cartridge.oraresident.JFunctions.sendUserFuncBatch(java.lang.String,
		java.lang.String, java.lang.String, java.lang.String,
		java.lang.String, java.lang.String, java.lang.String, long)';

  PROCEDURE register_user_func(opname VARCHAR2, params VARCHAR2, class VARCHAR2)
	AS LANGUAGE JAVA NAME 'chemaxon.jchem.cartridge.oraresident.JFunctions.registerUserOps(java.lang.String, 
				java.lang.String, java.lang.String)';

  FUNCTION is_jchem_table(indexSchema VARCHAR2, idxName VARCHAR2,
	tableSchema VARCHAR2, tblName VARCHAR2) RETURN VARCHAR2 AS
	LANGUAGE JAVA NAME 'chemaxon.jchem.cartridge.oraresident.JFunctions.isJChemTable(
        java.lang.String, java.lang.String, java.lang.String, java.lang.String)
		return java.lang.String';

  PROCEDURE invert_result(res IN OUT NOCOPY CD_ID_ARRAY, tableName VARCHAR2);

  FUNCTION getSqlForFormulaScan(idxSchema VARCHAR2, idxName VARCHAR2,
                 idxPartition VARCHAR2, query VARCHAR2, strt VARCHAR2)
      RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME
    'chemaxon.jchem.cartridge.oraresident.JFunctions.getSqlForFormulaScan(java.lang.String, java.lang.String, 
        java.lang.String, java.lang.String, java.lang.String)
        return java.lang.String';

  FUNCTION getSqlForMolweightScan(idxSchema VARCHAR2, idxName VARCHAR2,
                 idxPartition VARCHAR2, strt VARCHAR2, stop VARCHAR2, opFlag NUMBER)
      RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME
    'chemaxon.jchem.cartridge.oraresident.JFunctions.getSqlForMolweightScan(java.lang.String,
        java.lang.String, java.lang.String, java.lang.String, java.lang.String,
        int) return java.lang.String';

  FUNCTION react(reaction VARCHAR2, reactant1 VARCHAR2, reactant2 VARCHAR2,
                 reactant3 VARCHAR2, reactant4 VARCHAR2, options VARCHAR2,
                 old_jc_react VARCHAR)
      RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
        'chemaxon.jchem.cartridge.oraresident.JFunctions.react(java.lang.String,
        java.lang.String, java.lang.String, java.lang.String, java.lang.String,
        java.lang.String, java.lang.String) return java.lang.String';

  FUNCTION react_arr(reaction VARCHAR2, reactant1 VARCHAR2, reactant2 VARCHAR2,
                 reactant3 VARCHAR2, reactant4 VARCHAR2, options VARCHAR2,
                 old_jc_react VARCHAR)
      RETURN CHAR_PRODUCT_ARRAY PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
        'chemaxon.jchem.cartridge.oraresident.JFunctions.reactArr(java.lang.String,
        java.lang.String, java.lang.String, java.lang.String, java.lang.String,
        java.lang.String, java.lang.String) return oracle.sql.ARRAY';

  FUNCTION standardize(structure VARCHAR2, param VARCHAR2)
      RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME
              'chemaxon.jchem.cartridge.oraresident.JFunctions.standardize(
            java.lang.String, java.lang.String) return java.lang.String';

  FUNCTION formula_search(target VARCHAR2, query VARCHAR2, searchType VARCHAR2,
                          result out VARCHAR2)
      RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME
              'chemaxon.jchem.cartridge.oraresident.JccFormulaSearch.isMatching(
               java.lang.String, java.lang.String, java.lang.String,
               java.lang.String[]) return java.lang.String';

  PROCEDURE insert_into_idxtbl(idxtbl_name IN OUT NOCOPY VARCHAR2,
                               cdid_seq_name IN OUT NOCOPY VARCHAR2,
                               molprops IN OUT NOCOPY MOLPROPS_ARRAY_ARRAY);

  FUNCTION create_parse_insidxtbl_cur(idxtbl_name IN OUT NOCOPY VARCHAR2,
                                      cdid_seq_name IN OUT NOCOPY VARCHAR2,
                                      molprops IN OUT NOCOPY MOLPROPS_ARRAY) RETURN INTEGER;

  PROCEDURE insert_into_idxtbl_single(idxtbl_name IN OUT NOCOPY VARCHAR2,
                                     cdid_seq_name IN OUT NOCOPY VARCHAR2,
                                     molprops IN OUT NOCOPY MOLPROPS_ARRAY,
                                     cn INTEGER);

  function get_idx_stats(idx_schema varchar2, idx_name varchar2,
                         idx_partition varchar2) return varchar2
    as language java name
    'chemaxon.jchem.cartridge.oraresident.Indexing.getIndexStatistics(java.lang.String,
     java.lang.String, java.lang.String) return java.lang.String';
     
   function get_idx_stats(idx_schema varchar2, idx_name varchar2,
                         idx_partition varchar2, usrAssOpId varchar2) return varchar2
    as language java name
    'chemaxon.jchem.cartridge.oraresident.Indexing.getIndexStatistics(java.lang.String,
     java.lang.String, java.lang.String, java.lang.String) return java.lang.String';

  FUNCTION suspend_idx_update(idx_schema VARCHAR2, idx_name VARCHAR2,
                              idx_partition VARCHAR2, options VARCHAR2)
      RETURN VARCHAR2 AS LANGUAGE JAVA NAME
      'chemaxon.jchem.cartridge.oraresident.Indexing.suspendIndexUpdate(java.lang.String,
      java.lang.String, java.lang.String, java.lang.String) return java.lang.String';

  PROCEDURE resume_idx_update(idx_schema VARCHAR2, idx_name VARCHAR2,
                                 options VARCHAR2) AS LANGUAGE JAVA NAME
      'chemaxon.jchem.cartridge.oraresident.Indexing.resumeIndexUpdate(java.lang.String,
      java.lang.String, java.lang.String)';

  procedure purge_connection_cache as language java name
      'chemaxon.jchem.cartridge.oraresident.JFunctions.purgeConnectionCache()';

  function get_gmemutil_arr return gmem_util_array
      as language java name
      'chemaxon.jchem.cartridge.oraresident.JFunctions.getGlobalMemUtil()
        return oracle.sql.ARRAY';

  function get_taskinfo_arr return taskinfo_array
      as language java name
      'chemaxon.jchem.cartridge.oraresident.JFunctions.getTaskInfos()
        return oracle.sql.ARRAY';

  procedure load_cache(idxSchema IN VARCHAR2, idxName IN VARCHAR2, partition IN VARCHAR2, options IN VARCHAR2) as language java name
      'chemaxon.jchem.cartridge.oraresident.AdminSp.loadCache(java.lang.String, java.lang.String, java.lang.String, java.lang.String)';	 

  procedure unload_cache(idxSchema IN VARCHAR2, idxName IN VARCHAR2, partition IN VARCHAR2) as language java name
      'chemaxon.jchem.cartridge.oraresident.AdminSp.unLoadCache(java.lang.String, java.lang.String, java.lang.String)';	 
  
  function get_error_table_name(idx_schema varchar2, idx_name varchar2, idx_partition varchar2) return varchar2
      as language java name
      'chemaxon.jchem.cartridge.oraresident.Indexing.getErrorTableName(java.lang.String,
      java.lang.String, java.lang.String) return java.lang.String';
	  
  FUNCTION fuse(first_structure IN VARCHAR2, second_structure IN VARCHAR2, inputformat IN VARCHAR2)
      RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME
              'chemaxon.jchem.cartridge.oraresident.JFunctions.fuse(
            java.lang.String, java.lang.String, java.lang.String) return java.lang.String';	  

END jchem_core_pkg;
/
show errors;

------------------
-- jchem_clob_pkg
------------------
CREATE OR REPLACE PACKAGE jchem_clob_pkg AUTHID CURRENT_USER AS

  FUNCTION exec_function_cv(sqlOperator VARCHAR2, target CLOB, query VARCHAR2, 
        options VARCHAR2, rid VARCHAR2, idxSchema VARCHAR2, idxName VARCHAR2,
        idxPartition VARCHAR2, tblSchema VARCHAR2, tblName VARCHAR2, colName VARCHAR2, colType VARCHAR2,
        scanId number, scanFlg in number, result out varchar2)
    RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
    'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.execFunction(
        java.lang.String, oracle.sql.CLOB, java.lang.String, java.lang.String,
        java.lang.String, java.lang.String, java.lang.String,  java.lang.String, java.lang.String,
        java.lang.String, java.lang.String, java.lang.String, int, int, java.lang.String[])
      return java.lang.String';

  FUNCTION exec_function_vc(sqlOperator VARCHAR2, target VARCHAR2, query CLOB, 
        options VARCHAR2, rid VARCHAR2, idxSchema VARCHAR2, idxName VARCHAR2,
        idxPartition VARCHAR2, tblSchema VARCHAR2, tblName VARCHAR2, colName VARCHAR2, colType VARCHAR2,
        scanId number, scanFlg in number, result out varchar2)
    RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
    'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.execFunction(java.lang.String,
    java.lang.String, oracle.sql.CLOB, java.lang.String, java.lang.String,
    java.lang.String, java.lang.String,  java.lang.String, java.lang.String, java.lang.String, 
    java.lang.String, java.lang.String, int, int, java.lang.String[]) return java.lang.String';

  FUNCTION exec_function_cb(sqlOperator VARCHAR2, target CLOB, query BLOB, 
        options VARCHAR2, rid VARCHAR2, idxSchema VARCHAR2, idxName VARCHAR2,
        idxPartition VARCHAR2, tblSchema VARCHAR2, tblName VARCHAR2, colName VARCHAR2, colType VARCHAR2, 
        scanId number, scanFlg in number, result out varchar2)
    RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
    'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.execFunction(java.lang.String,
    oracle.sql.CLOB, oracle.sql.BLOB, java.lang.String, java.lang.String,
    java.lang.String, java.lang.String,  java.lang.String, java.lang.String, java.lang.String,
    java.lang.String, java.lang.String, int, int, java.lang.String[]) return java.lang.String';

  FUNCTION exec_function(sqlOperator VARCHAR2, target CLOB, query CLOB, 
        options VARCHAR2, rid VARCHAR2, idxSchema VARCHAR2, idxName VARCHAR2,
        idxPartition VARCHAR2, tblSchema VARCHAR2, tblName VARCHAR2, colName VARCHAR2, colType VARCHAR2,
        scanId number, scanFlg in number, result out varchar2)
    RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
    'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.execFunction(java.lang.String,
    oracle.sql.CLOB, oracle.sql.CLOB, java.lang.String, java.lang.String,
    java.lang.String, java.lang.String,  java.lang.String, java.lang.String, java.lang.String, 
    java.lang.String, java.lang.String, int, int, java.lang.String[]) return java.lang.String';

  FUNCTION exec_function__c(sqlOperator VARCHAR2, query CLOB, target CLOB, 
        options VARCHAR2, rid VARCHAR2, idxSchema VARCHAR2, idxName VARCHAR2,
        idxPartition VARCHAR2, tblSchema VARCHAR2, tblName VARCHAR2, colName VARCHAR2, colType VARCHAR2,
        scanId number, scanFlg in number, result out CLOB)
    RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
    'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.execFunctionC(java.lang.String,
    oracle.sql.CLOB, oracle.sql.CLOB, java.lang.String, java.lang.String,
    java.lang.String, java.lang.String,  java.lang.String, java.lang.String, java.lang.String, 
    java.lang.String, java.lang.String, int, int, oracle.sql.CLOB[]) return java.lang.String';

  FUNCTION evaluate_arr(target CLOB, options VARCHAR2, rid VARCHAR2,
        idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2,
        tblSchema VARCHAR2, tblName VARCHAR2)
    RETURN CLOB_ARRAY PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
    'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.evaluateArr(oracle.sql.CLOB,
      java.lang.String, java.lang.String, java.lang.String,
      java.lang.String, java.lang.String, java.lang.String,
      java.lang.String) return oracle.sql.CLOB[]';

  FUNCTION hitColorAndAlign(tblSchema VARCHAR2, tblName VARCHAR2,
              colName VARCHAR2, query CLOB, rowids VARCHAR2,
              options VARCHAR2, hitColorAndAlignOptions VARCHAR2,
              scanId number)
              RETURN CLOB_ARRAY PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
    'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.hitColorAndAlignOptions(
      java.lang.String, java.lang.String, java.lang.String,
      oracle.sql.CLOB, java.lang.String, java.lang.String, java.lang.String, long)
      return oracle.sql.CLOB[]';

  FUNCTION index_scan(indexSchema VARCHAR2, idxName VARCHAR2,
                       indexPartition VARCHAR2, tblSchema VARCHAR2,
                       tblname VARCHAR2, colName VARCHAR2, colType VARCHAR2, optype VARCHAR2,
                       opflavor VARCHAR2, strt NUMBER, stop NUMBER, opFlag
                       number, query clob, options VARCHAR2, scanId number)
  RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME
	'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.indexScan(java.lang.String, 
			java.lang.String, java.lang.String, java.lang.String, java.lang.String,
			java.lang.String, java.lang.String, java.lang.String, 
                        java.lang.String, java.lang.Double, java.lang.Double,
                        int, oracle.sql.CLOB,
                        java.lang.String, long) return java.lang.String'; 

  FUNCTION get_hit_count(tblSchema VARCHAR2, tblname VARCHAR2, colName VARCHAR2, 
                         query CLOB, options VARCHAR2)
        RETURN NUMBER PARALLEL_ENABLE AS LANGUAGE JAVA NAME
	'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.getHitCount(java.lang.String,  
			java.lang.String, java.lang.String, oracle.sql.CLOB,
                        java.lang.String) return int'; 

  FUNCTION similaritySearch(query CLOB, stype VARCHAR2,
            strt NUMBER, stop NUMBER, opFlag number, searchOptions
            VARCHAR2, indexSchema VARCHAR2, idxName VARCHAR2, idxPartition
            VARCHAR2, tblSchema VARCHAR2, tblname VARCHAR2, colName VARCHAR2,
            colType VARCHAR2, scanId NUMBER) RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA
            NAME 
	'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.getSimilarity(
        oracle.sql.CLOB, java.lang.String, java.lang.Double, java.lang.Double, int,
        java.lang.String, java.lang.String, java.lang.String,
        java.lang.String, java.lang.String,
        java.lang.String, java.lang.String, java.lang.String, long)
          return java.lang.String';
  
  FUNCTION insert_mol_into_idxtable(str CLOB, 
		indexSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2,
                tblSchema VARCHAR2, tblname VARCHAR2, colName VARCHAR2, colType VARCHAR2, rid VARCHAR2)
        RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME
		'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.insertMolIntoIndexTable(
		oracle.sql.CLOB, java.lang.String, java.lang.String, java.lang.String,
		java.lang.String, java.lang.String, java.lang.String, java.lang.String, java.lang.String)
                return java.lang.String';

  FUNCTION update_mol_idxtable(oldval CLOB, str CLOB, 
                  indexSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2,
                  tblSchema VARCHAR2, tblname VARCHAR2, colName VARCHAR2, colType VARCHAR2, 
                  rid VARCHAR2)
        RETURN VARCHAR2
        PARALLEL_ENABLE AS LANGUAGE JAVA NAME
		'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.updateMolIndexTable(
                 oracle.sql.CLOB,
                 oracle.sql.CLOB, java.lang.String, java.lang.String, java.lang.String,
		 java.lang.String, java.lang.String, java.lang.String, java.lang.String, java.lang.String)
                return java.lang.String';

  FUNCTION calc_molProp(query CLOB, propName VARCHAR2, rid VARCHAR2, 
		idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2, 
        tblSchema VARCHAR2, tableName VARCHAR2) RETURN CLOB
        PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
		'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.calcMolProp(
        oracle.sql.CLOB, java.lang.String, java.lang.String, 
		java.lang.String, java.lang.String, java.lang.String,
		java.lang.String, java.lang.String) return oracle.sql.CLOB';

  FUNCTION calc_molPropNum(query CLOB, propName VARCHAR2, rid VARCHAR2, 
		idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2, 
        tblSchema VARCHAR2, tableName VARCHAR2) RETURN NUMBER
        PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
		'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.calcMolPropNum(
        oracle.sql.CLOB, java.lang.String, java.lang.String, 
		java.lang.String, java.lang.String, java.lang.String,
		java.lang.String, java.lang.String) return java.lang.Double';

  FUNCTION molconvertc(query VARCHAR2, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2, rid VARCHAR2,
        idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2, 
        tblSchema VARCHAR2, tableName VARCHAR2, tmpBlob CLOB)
        RETURN CLOB PARALLEL_ENABLE AS LANGUAGE JAVA NAME
        'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.molConvertCV(
        java.lang.String, java.lang.String,java.lang.String,java.lang.String,
        java.lang.String, java.lang.String, java.lang.String, java.lang.String,
        java.lang.String, java.lang.String, oracle.sql.CLOB)
        return oracle.sql.CLOB';

  FUNCTION molconvertc(query CLOB, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2, rid VARCHAR2,
        idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2, 
        tblSchema VARCHAR2, tableName VARCHAR2, tmpBlob CLOB)
        RETURN CLOB PARALLEL_ENABLE AS LANGUAGE JAVA NAME
        'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.molConvertCC(
        oracle.sql.CLOB, java.lang.String,java.lang.String,java.lang.String,
        java.lang.String, java.lang.String, java.lang.String, java.lang.String,
        java.lang.String, java.lang.String, oracle.sql.CLOB)
        return oracle.sql.CLOB';

  FUNCTION molconvertc(query BLOB, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2, rid VARCHAR2,
        idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2, 
        tblSchema VARCHAR2, tableName VARCHAR2, tmpBlob CLOB)
        RETURN CLOB PARALLEL_ENABLE AS LANGUAGE JAVA NAME
        'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.molConvertCB(
        oracle.sql.BLOB, java.lang.String,java.lang.String,java.lang.String,
        java.lang.String, java.lang.String, java.lang.String, java.lang.String,
        java.lang.String, java.lang.String, oracle.sql.CLOB)
        return oracle.sql.CLOB';

  FUNCTION molconvertc_from_idx(rid VARCHAR2, type VARCHAR2,
        idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2, 
        tblSchema VARCHAR2, tableName VARCHAR2, tmpBlob CLOB)
        RETURN CLOB PARALLEL_ENABLE AS LANGUAGE JAVA NAME
        'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.molConvertCFromRowId(
        java.lang.String, java.lang.String, java.lang.String, java.lang.String,
        java.lang.String, java.lang.String, java.lang.String, oracle.sql.CLOB)
        return oracle.sql.CLOB';
  
  FUNCTION get_molweight(query CLOB, rid VARCHAR2, 
		idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2,
        tblSchema VARCHAR2, tableName VARCHAR2) RETURN NUMBER
        PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
		'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.getMolweight(
        oracle.sql.CLOB, java.lang.String,java.lang.String,
        java.lang.String, java.lang.String, java.lang.String,
        java.lang.String) return java.lang.Double';

  FUNCTION get_molformula(query CLOB, rid VARCHAR2, 
		idxSchema VARCHAR2, idxName VARCHAR2,  idxPartition VARCHAR2,
        tblSchema VARCHAR2, tableName VARCHAR2) RETURN VARCHAR2
        PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
		'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.getMolformula(
        oracle.sql.CLOB, java.lang.String, java.lang.String,
        java.lang.String, java.lang.String, java.lang.String,
        java.lang.String) return java.lang.String';

  FUNCTION calc_molProp_from_idx(rid VARCHAR2, propName VARCHAR2, 
	idxSchema VARCHAR2, idxName VARCHAR2,  idxPartition VARCHAR2,
    tblSchema VARCHAR2, tableName VARCHAR2) RETURN CLOB
    PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
	'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.calcMolPropFromRowid(
    java.lang.String, java.lang.String, java.lang.String, java.lang.String,
    java.lang.String, java.lang.String, java.lang.String)
	return oracle.sql.CLOB';

  FUNCTION calc_molPropNum_from_idx(rid VARCHAR2, propName VARCHAR2, 
	idxSchema VARCHAR2, idxName VARCHAR2,  idxPartition VARCHAR2,
    tblSchema VARCHAR2, tableName VARCHAR2) RETURN NUMBER
    PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
	'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.calcMolPropNumFromRowid(
    java.lang.String, java.lang.String, java.lang.String, java.lang.String,
    java.lang.String, java.lang.String, java.lang.String)
	return java.lang.Double';

  FUNCTION send_user_func(name VARCHAR2, delim VARCHAR2, 
	params CLOB) RETURN CLOB AS LANGUAGE JAVA NAME 
		'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.sendUserFunc(
        java.lang.String, java.lang.String, oracle.sql.CLOB)
        return oracle.sql.CLOB';

  PROCEDURE send_user_func_batch(opName VARCHAR2, fieldName VARCHAR2, 
		operator_str VARCHAR2, params CLOB, indexSchema VARCHAR2,
        tableSchema VARCHAR2, tableName VARCHAR2,
        idxSchema VARCHAR2, scanId NUMBER) AS LANGUAGE JAVA NAME
	'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.sendUserFuncBatch(
        java.lang.String, java.lang.String,
		oracle.sql.CLOB, java.lang.String, java.lang.String,
		java.lang.String, java.lang.String, long)';

  FUNCTION getSqlForFormulaScan(idxSchema VARCHAR2, idxName VARCHAR2,
                 idxPartition VARCHAR2, query CLOB, pred VARCHAR2)
      RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME
    'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.getSqlForFormulaScan(java.lang.String, java.lang.String, 
    java.lang.String, oracle.sql.CLOB, java.lang.String)
    return java.lang.String';

  FUNCTION react(reaction CLOB, reactant1 CLOB, reactant2 CLOB, reactant3 CLOB,
          reactant4 CLOB, options VARCHAR2, tempBlob CLOB, old_jc_react VARCHAR)
          RETURN CLOB PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
        'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.react( oracle.sql.CLOB,
        oracle.sql.CLOB, oracle.sql.CLOB, oracle.sql.CLOB, oracle.sql.CLOB,
        java.lang.String, oracle.sql.CLOB, java.lang.String)
        return oracle.sql.CLOB';

  FUNCTION react_arr(reaction CLOB, reactant1 CLOB, reactant2 CLOB, reactant3 CLOB,
          reactant4 CLOB, options VARCHAR2, tempBlob CLOB)
          RETURN CLOB_PRODUCT_ARRAY PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
        'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.reactArr(oracle.sql.CLOB,
        oracle.sql.CLOB, oracle.sql.CLOB, oracle.sql.CLOB, oracle.sql.CLOB,
        java.lang.String, oracle.sql.CLOB)
        return oracle.sql.ARRAY';

  FUNCTION standardize(structure CLOB, param VARCHAR2, temp_blob CLOB)
      RETURN CLOB PARALLEL_ENABLE AS LANGUAGE JAVA NAME
      'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.standardize(
      oracle.sql.CLOB, java.lang.String, oracle.sql.CLOB)
      return oracle.sql.CLOB';

  FUNCTION fuse(first_structure IN CLOB, second_structure IN CLOB, inputformat IN VARCHAR2, tmpClob CLOB)
      RETURN CLOB PARALLEL_ENABLE AS LANGUAGE JAVA NAME
              'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.fuse(
            oracle.sql.CLOB, oracle.sql.CLOB, java.lang.String, oracle.sql.CLOB) return oracle.sql.CLOB';	  
		
  FUNCTION equals(b1 CLOB, b2 CLOB) RETURN NUMBER PARALLEL_ENABLE
        AS LANGUAGE JAVA NAME 
		'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.equals(
        oracle.sql.CLOB, oracle.sql.CLOB) return int';
		
END jchem_clob_pkg;
/
show errors;

------------------
-- jchem_blob_pkg
------------------
CREATE OR REPLACE PACKAGE jchem_blob_pkg AUTHID CURRENT_USER AS

  FUNCTION exec_function_vb(sqlOperator VARCHAR2, target VARCHAR2, query BLOB,
        options VARCHAR2, rid VARCHAR2, idxSchema VARCHAR2,
        idxName VARCHAR2, idxPartition VARCHAR2, 
        tblSchema VARCHAR2, tblName VARCHAR2, colName VARCHAR2, colType VARCHAR2, 
        scanId number, scanFlg in number, result out varchar2)
    RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
    'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.execFunction(
      java.lang.String, java.lang.String, oracle.sql.BLOB, java.lang.String,
      java.lang.String, java.lang.String, java.lang.String,  java.lang.String, java.lang.String, 
      java.lang.String, java.lang.String, java.lang.String, int, int, java.lang.String[])
      return java.lang.String';

  FUNCTION exec_function_bv(sqlOperator VARCHAR2, target BLOB, query VARCHAR2, 
        options VARCHAR2, rid VARCHAR2, idxSchema VARCHAR2, idxName VARCHAR2,
        idxPartition VARCHAR2, tblSchema VARCHAR2, tblName VARCHAR2, colName VARCHAR2, colType VARCHAR2,
        scanId number, scanFlg in number, result out varchar2)
    RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
    'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.execFunction(java.lang.String,
    oracle.sql.BLOB, java.lang.String, java.lang.String, java.lang.String,
    java.lang.String, java.lang.String, java.lang.String, java.lang.String, java.lang.String, 
    java.lang.String, java.lang.String, int, int,
    java.lang.String[]) return java.lang.String';

  FUNCTION exec_function(sqlOperator VARCHAR2, query BLOB, target BLOB, 
        options VARCHAR2, rid VARCHAR2, idxSchema VARCHAR2, idxName VARCHAR2,
        idxPartition VARCHAR2, tblSchema VARCHAR2, tblName VARCHAR2, colName VARCHAR2, colType VARCHAR2,
        scanId number, scanFlg in number, result out varchar2)
    RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
    'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.execFunction(java.lang.String,
    oracle.sql.BLOB, oracle.sql.BLOB, java.lang.String, java.lang.String,
    java.lang.String, java.lang.String,  java.lang.String, java.lang.String, java.lang.String, 
    java.lang.String, java.lang.String, int, int, java.lang.String[]) return java.lang.String';

  FUNCTION exec_function__b(sqlOperator VARCHAR2, query BLOB, target BLOB, 
        options VARCHAR2, rid VARCHAR2, idxSchema VARCHAR2, idxName VARCHAR2,
        idxPartition VARCHAR2, tblSchema VARCHAR2, tblName VARCHAR2, colName VARCHAR2, colType VARCHAR2,
        scanId number, scanFlg in number, result out BLOB)
    RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
    'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.execFunctionB(java.lang.String,
    oracle.sql.BLOB, oracle.sql.BLOB, java.lang.String, java.lang.String,
    java.lang.String, java.lang.String,  java.lang.String, java.lang.String, java.lang.String, 
    java.lang.String, java.lang.String, int, int, oracle.sql.BLOB[]) return java.lang.String';

  FUNCTION evaluate_arr(target BLOB, options VARCHAR2, rid VARCHAR2,
        idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2,
        tblSchema VARCHAR2, tblName VARCHAR2)
    RETURN BLOB_ARRAY PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
    'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.evaluateArr(oracle.sql.BLOB,
      java.lang.String, java.lang.String, java.lang.String,
      java.lang.String, java.lang.String, java.lang.String,
      java.lang.String) return oracle.sql.BLOB[]';


  FUNCTION hitColorAndAlign(tblSchema VARCHAR2, tblName VARCHAR2,
              colName VARCHAR2, query BLOB, rowids VARCHAR2,
              options VARCHAR2, hitColorAndAlignOptions VARCHAR2,
              scanId number)
              RETURN BLOB_ARRAY PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
    'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.hitColorAndAlignOptions(
      java.lang.String, java.lang.String, java.lang.String,
      oracle.sql.BLOB, java.lang.String, java.lang.String, java.lang.String, long)
      return oracle.sql.BLOB[]';

  FUNCTION index_scan(indexSchema VARCHAR2, idxName VARCHAR2,
                       indexPartition VARCHAR2, tblSchema VARCHAR2,
                       tblname VARCHAR2, colName VARCHAR2, colType VARCHAR2, optype VARCHAR2,
                       opflavor VARCHAR2, strt NUMBER, stop NUMBER, opFlag number,
                       query blob, options VARCHAR2, scanId number)
      RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME
	'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.indexScan(java.lang.String, 
			java.lang.String, java.lang.String, java.lang.String, java.lang.String,
			java.lang.String, java.lang.String, java.lang.String, 
                        java.lang.String, java.lang.Double, java.lang.Double,
                        int, oracle.sql.BLOB, java.lang.String, long)
          return java.lang.String'; 

  FUNCTION get_hit_count(tblSchema VARCHAR2, tblname VARCHAR2, colName VARCHAR2, 
                         query BLOB, options VARCHAR2)
        RETURN NUMBER PARALLEL_ENABLE AS LANGUAGE JAVA NAME
	'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.getHitCount(java.lang.String, 
			java.lang.String, java.lang.String, oracle.sql.BLOB,
                        java.lang.String) return int'; 

  FUNCTION similaritySearch(query BLOB, stype VARCHAR2,
                strt NUMBER, stop NUMBER, opFlag number, searchOptions
                VARCHAR2, indexSchema VARCHAR2, idxName VARCHAR2, idxPartition
                VARCHAR2, tblSchema VARCHAR2, tblname VARCHAR2, colName
                VARCHAR2, colType VARCHAR2, scanId NUMBER) RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
	'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.getSimilarity(
        oracle.sql.BLOB, java.lang.String, java.lang.Double, java.lang.Double, int,
        java.lang.String, java.lang.String, java.lang.String, java.lang.String,
        java.lang.String, java.lang.String, java.lang.String, java.lang.String, long) return
        java.lang.String';
  
  FUNCTION insert_mol_into_idxtable(str BLOB, 
		indexSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2,
        tblSchema VARCHAR2, tblname VARCHAR2, colName VARCHAR2, colType VARCHAR2, rid VARCHAR2)
        RETURN VARCHAR2
        PARALLEL_ENABLE AS LANGUAGE JAVA NAME
		'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.insertMolIntoIndexTable(
		oracle.sql.BLOB, java.lang.String, java.lang.String, java.lang.String,
		java.lang.String, java.lang.String, java.lang.String, java.lang.String, java.lang.String)
                return java.lang.String';

  FUNCTION update_mol_idxtable(oldval BLOB, str BLOB, 
		indexSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2,
                tblSchema VARCHAR2, tblname VARCHAR2, colName VARCHAR2, colType VARCHAR2,
                rid VARCHAR2)
        RETURN VARCHAR2
        PARALLEL_ENABLE AS LANGUAGE JAVA NAME
		'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.updateMolIndexTable(
                 oracle.sql.BLOB,
                 oracle.sql.BLOB, java.lang.String, java.lang.String, java.lang.String,
		 java.lang.String, java.lang.String, java.lang.String, java.lang.String, java.lang.String)
                return java.lang.String';

  FUNCTION calc_molProp(query BLOB, propName VARCHAR2, rid VARCHAR2, 
		idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2, 
        tblSchema VARCHAR2, tableName VARCHAR2) RETURN BLOB
        PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
		'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.calcMolProp(
        oracle.sql.BLOB, java.lang.String, java.lang.String, 
		java.lang.String, java.lang.String, java.lang.String,
		java.lang.String, java.lang.String) return oracle.sql.BLOB';

  FUNCTION calc_molPropNum(query BLOB, propName VARCHAR2, rid VARCHAR2, 
		idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2, 
        tblSchema VARCHAR2, tableName VARCHAR2) RETURN NUMBER
        PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
		'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.calcMolPropNum(
        oracle.sql.BLOB, java.lang.String, java.lang.String, 
		java.lang.String, java.lang.String, java.lang.String,
		java.lang.String, java.lang.String) return java.lang.Double';

  FUNCTION molconvertb(query VARCHAR2, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2, rid VARCHAR2,
        idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2,
        tblSchema VARCHAR2, tableName VARCHAR2, tmpBlob BLOB)
        RETURN BLOB PARALLEL_ENABLE AS LANGUAGE JAVA NAME
        'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.molConvertBV(java.lang.String,
        java.lang.String, java.lang.String, java.lang.String, java.lang.String,
        java.lang.String, java.lang.String, java.lang.String, 
        java.lang.String, java.lang.String, oracle.sql.BLOB) 
		return oracle.sql.BLOB';

  FUNCTION molconvertb(query CLOB, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2, rid VARCHAR2,
        idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2, 
        tblSchema VARCHAR2, tableName VARCHAR2, tmpBlob BLOB)
        RETURN BLOB PARALLEL_ENABLE AS LANGUAGE JAVA NAME
        'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.molConvertBC(
        oracle.sql.CLOB, java.lang.String,java.lang.String,java.lang.String,
        java.lang.String, java.lang.String, java.lang.String, java.lang.String,
        java.lang.String, java.lang.String, oracle.sql.BLOB)
        return oracle.sql.BLOB';

  FUNCTION molconvertb(query BLOB, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2, rid VARCHAR2,
        idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2, 
        tblSchema VARCHAR2, tableName VARCHAR2, tmpBlob BLOB)
        RETURN BLOB PARALLEL_ENABLE AS LANGUAGE JAVA NAME
        'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.molConvertBB(
        oracle.sql.BLOB, java.lang.String,java.lang.String,
        java.lang.String, java.lang.String,
        java.lang.String, java.lang.String, java.lang.String,
        java.lang.String, java.lang.String, oracle.sql.BLOB)
        return oracle.sql.BLOB';

  FUNCTION molconvertb_from_idx(rid VARCHAR2, type VARCHAR2,
        idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2, 
        tblSchema VARCHAR2, tableName VARCHAR2, tmpBlob BLOB)
        RETURN BLOB PARALLEL_ENABLE AS LANGUAGE JAVA NAME
        'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.molConvertBFromRowId(
        java.lang.String, java.lang.String, java.lang.String,
        java.lang.String, java.lang.String, java.lang.String, java.lang.String, oracle.sql.BLOB)
        return oracle.sql.BLOB';
  
  FUNCTION get_molweight(query BLOB, rid VARCHAR2, 
		idxSchema VARCHAR2, idxName VARCHAR2, idxPartition VARCHAR2,
        tblSchema VARCHAR2, tableName VARCHAR2) RETURN NUMBER
        PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
		'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.getMolweight(
        oracle.sql.BLOB, java.lang.String,java.lang.String,
        java.lang.String, java.lang.String, java.lang.String,
        java.lang.String) return java.lang.Double';

  FUNCTION get_molformula(query BLOB, rid VARCHAR2, 
		idxSchema VARCHAR2, idxName VARCHAR2,  idxPartition VARCHAR2,
        tblSchema VARCHAR2, tableName VARCHAR2) RETURN VARCHAR2
        PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
		'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.getMolformula(
        oracle.sql.BLOB, java.lang.String, java.lang.String,
        java.lang.String, java.lang.String, java.lang.String,
        java.lang.String) return java.lang.String';

  FUNCTION calc_molProp_from_idx(rid VARCHAR2, propName VARCHAR2, 
	idxSchema VARCHAR2, idxName VARCHAR2,  idxPartition VARCHAR2,
    tblSchema VARCHAR2, tableName VARCHAR2) RETURN BLOB
    PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
	'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.calcMolPropFromRowid(
    java.lang.String, java.lang.String, java.lang.String, java.lang.String,
    java.lang.String, java.lang.String, java.lang.String)
	return oracle.sql.BLOB';

  FUNCTION calc_molPropNum_from_idx(rid VARCHAR2, propName VARCHAR2, 
	idxSchema VARCHAR2, idxName VARCHAR2,  idxPartition VARCHAR2,
    tblSchema VARCHAR2, tableName VARCHAR2) RETURN NUMBER
    PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
	'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.calcMolPropFromRowidNum(
    java.lang.String, java.lang.String, java.lang.String, java.lang.String,
    java.lang.String, java.lang.String, java.lang.String)
	return java.lang.Double';

  FUNCTION send_user_func(name VARCHAR2, delim VARCHAR2, 
	params BLOB) RETURN BLOB AS LANGUAGE JAVA NAME 
		'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.sendUserFunc(
        java.lang.String, java.lang.String, oracle.sql.BLOB)
        return oracle.sql.BLOB';

  PROCEDURE send_user_func_batch(opName VARCHAR2, fieldName VARCHAR2, 
		operator_str VARCHAR2, params BLOB, indexSchema VARCHAR2,
        tableSchema VARCHAR2, tableName VARCHAR2,
        idxSchema VARCHAR2, scanId NUMBER) AS LANGUAGE JAVA NAME
	'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.sendUserFuncBatch(
        java.lang.String, java.lang.String,
		oracle.sql.BLOB, java.lang.String, java.lang.String,
		java.lang.String, java.lang.String, long)';

  FUNCTION getSqlForFormulaScan(idxSchema VARCHAR2, idxName VARCHAR2,
		 idxPartition VARCHAR2, query BLOB, pred VARCHAR2)
      RETURN VARCHAR2 PARALLEL_ENABLE AS LANGUAGE JAVA NAME
    'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.getSqlForFormulaScan( java.lang.String, java.lang.String, 
        java.lang.String, oracle.sql.BLOB, java.lang.String) return java.lang.String';

  FUNCTION react(reaction BLOB, reactant1 BLOB, reactant2 BLOB, reactant3 BLOB,
          reactant4 BLOB, options VARCHAR2, tempBlob BLOB, old_jc_react VARCHAR)
          RETURN BLOB PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
        'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.react( oracle.sql.BLOB,
        oracle.sql.BLOB, oracle.sql.BLOB, oracle.sql.BLOB, oracle.sql.BLOB,
        java.lang.String, oracle.sql.BLOB, java.lang.String)
        return oracle.sql.BLOB';

  FUNCTION react_arr(reaction BLOB, reactant1 BLOB, reactant2 BLOB, reactant3 BLOB,
          reactant4 BLOB, options VARCHAR2, tempBlob BLOB)
          RETURN BLOB_PRODUCT_ARRAY PARALLEL_ENABLE AS LANGUAGE JAVA NAME 
        'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.reactArr(oracle.sql.BLOB,
        oracle.sql.BLOB, oracle.sql.BLOB, oracle.sql.BLOB, oracle.sql.BLOB,
        java.lang.String, oracle.sql.BLOB)
        return oracle.sql.ARRAY';

  FUNCTION standardize(structure BLOB, param VARCHAR2, temp_blob BLOB)
      RETURN BLOB PARALLEL_ENABLE AS LANGUAGE JAVA NAME
      'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.standardize(
      oracle.sql.BLOB, java.lang.String, oracle.sql.BLOB)
      return oracle.sql.BLOB';

  FUNCTION fuse(first_structure IN BLOB, second_structure IN BLOB, inputformat IN VARCHAR2, tmpBlob BLOB)
      RETURN BLOB PARALLEL_ENABLE AS LANGUAGE JAVA NAME
              'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.fuse(
            oracle.sql.BLOB, oracle.sql.BLOB, java.lang.String, oracle.sql.BLOB) return oracle.sql.BLOB';	  

	FUNCTION equals(b1 BLOB, b2 BLOB) RETURN NUMBER PARALLEL_ENABLE
        AS LANGUAGE JAVA NAME 
		'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.equals(
        oracle.sql.BLOB, oracle.sql.BLOB) return int';

END jchem_blob_pkg;
/
show errors;

CREATE OR REPLACE PACKAGE jchem_misc_pkg AUTHID CURRENT_USER AS

  is_profiling boolean;

  PROCEDURE put_time as language java name 'chemaxon.jchem.cartridge.sharedorajcsrv.JCartLogFunctions.putTime()';

  procedure cache_stmts(yesno number) as language java name
      'chemaxon.jchem.cartridge.oraresident.dbsession.JavaStoredProcSession.cacheStatements(int)';

  function get_cached_stmt_count return number as language java name
      'chemaxon.jchem.cartridge.oraresident.dbsession.JavaStoredProcSession.getCachedStmtCount() return int';

  FUNCTION get_jchemstreams_url RETURN VARCHAR2;

  PROCEDURE memstat(dir VARCHAR2, name VARCHAR2) AS LANGUAGE JAVA NAME
      'chemaxon.jchem.cartridge.oraresident.JFunctions.memstat(java.lang.String,
                                                   java.lang.String)';
  PROCEDURE setDbCallback(host VARCHAR2, port number, dbname VARCHAR2) AS
        LANGUAGE JAVA NAME
      'chemaxon.jchem.cartridge.oraresident.JFunctions.setDbCallback(java.lang.String, int,
      java.lang.String)';

  FUNCTION start_search_profiling RETURN NUMBER AS LANGUAGE JAVA NAME
      'chemaxon.jchem.cartridge.oraresident.JFunctions.startSearchProfiling() return long';

  PROCEDURE stop_search_profiling AS LANGUAGE JAVA NAME
       'chemaxon.jchem.cartridge.oraresident.JFunctions.stopSearchProfiling()';

  PROCEDURE trace(msg VARCHAR2) AS LANGUAGE JAVA NAME
        'chemaxon.jchem.cartridge.sharedorajcsrv.JCartLogFunctions.trace(java.lang.String)';

  FUNCTION get_time RETURN VARCHAR2 as language java name 
	'chemaxon.jchem.cartridge.sharedorajcsrv.JCartLogFunctions.getTime() return java.lang.String';

  FUNCTION get_time_diff(time_p VARCHAR2) RETURN NUMBER as language java name
    'chemaxon.jchem.cartridge.oraresident.JFunctions.getTimeDiff(
    java.lang.String) return long';

  PROCEDURE put_any_line (pstrText IN VARCHAR2) ;

  --Functions for reading and writing row types
	-- Offers conversions raw and char formats
  FUNCTION chartoraw(v_char varchar2) return long raw;
  FUNCTION rawtochar(v_raw long raw) return varchar2;

	-- Offers conversions between decimal and hex format
  FUNCTION numtohex(v_hex number) return varchar2;
  FUNCTION hextonum(v_hex varchar2) return number;

END jchem_misc_pkg;
/
show errors;


CREATE OR REPLACE PACKAGE BODY jchem_core_pkg AS

  procedure set_use_arrays(use_arrays varchar2) is 
  begin
    if use_arrays='Y' then
      useArrays:=true;
	elsif use_arrays='N' then
      useArrays:=false;
	end if;
  end;

  function test(password varchar2) return varchar2 is
  begin
    jchem_core_pkg.use_password(null, password);
    return jchem_core_pkg.getEnvironment();
  end;

  procedure handle_java_error(errm varchar2) is
    em VARCHAR2(32767);
    pos integer;
    key1 VARCHAR2(4000) := 'Java call terminated by uncaught Java exception: ';
    key2 VARCHAR2(4000) := 'Exception: ';
    enostr varchar2(3);
    eno number := -20101;
  begin
    if errm is null then
      return;
    end if;

    em := errm;

    -- Starts with   [1][0-9][0-9]~   ?
    pos := instr(em, '~');
--    jchem_misc_pkg.trace('Processing `' || em || '`');
    if (pos = 4) then
      enostr := substr(em, 1, 3);
      begin
        eno := to_number(enostr);
        eno := -20000 - eno;
        em := substr(em, 5);
      exception
        when others then
        eno := eno;
      end;
    end if;

    raise_application_error(eno, em);
  end;

  FUNCTION get_cartowner_schema RETURN VARCHAR2 IS
  BEGIN
    return '&1';
  END;

  FUNCTION getJChemVersion RETURN VARCHAR2 IS
  BEGIN
    return '&2';
  END;

  FUNCTION getTableVersion RETURN NUMBER IS
  BEGIN
    return &3;
  END;

  procedure use_password(password varchar2) is
    errm varchar2(32767);
  begin
    errm := use_password0(null, password);
    if errm is not null then
      jchem_core_pkg.handle_java_error(errm);
    end if;
  end;

  procedure use_password(login varchar2, password varchar2) is
    errm varchar2(32767);
  begin
    errm := use_password0(login, password);
    if errm is not null then
      jchem_core_pkg.handle_java_error(errm);
    end if;
  end;

  function set_password(password varchar2) return number is
    errm varchar2(32767);
  begin
    errm := use_password0(null, password);
    if errm is null then
      return 0;
    else
      return -1;
    end if;
  end;

  procedure use_default_account is
    errm varchar2(32767);
  begin
    errm := use_default_account0();
    if errm is not null then
      jchem_core_pkg.handle_java_error(errm);
    end if;
  end;

  PROCEDURE checkTableVersionEx(func_name VARCHAR2,
                                indexctx IN SYS.ODCIIndexCtx) IS
    ia          SYS.ODCIIndexInfo;
    schema_name VARCHAR2(255);
    table_name  VARCHAR2(255);
    dobject_id  NUMBER;
    err_num     NUMBER;
  BEGIN

    IF indexctx.IndexInfo IS NOT NULL THEN
      -- The first argument has an index on it

      ia := indexctx.IndexInfo;
      IF ia.IndexSchema IS NOT NULL THEN
        -- Has index
        jchem_core_pkg.checkTableVersion(
                  ia.IndexSchema,
                  ia.IndexName,
                  ia.IndexPartition,
                  ia.IndexCols(1).TableSchema,
                  ia.IndexCols(1).TableName,
                  ia.IndexCols(1).ColName,
				  ia.IndexCols(1).ColTypeName);
      ELSE
        raise_application_error(-20102,
              'checkTableVersionEx: indexctx.IndexInfo NOT IS NULL, but indexctx.IndexInfo.IndexSchema IS NULL');
      END IF;
    ELSE 
        -- Does not have index
        IF indexctx.rid IS NOT NULL THEN
          -- The first argument is a column
          dobject_id := DBMS_ROWID.ROWID_OBJECT(indexctx.rid);

--            SELECT owner, object_name INTO schema_name, table_name
--                  FROM sys.dba_objects WHERE DATA_OBJECT_ID = dobject_id;

          raise_application_error(-20101,
                  'Please, create domain index on the column referenced in the operator '
                  || func_name || ' of the table with DATA_OBJECT_ID=' || dobject_id
                  || '. You can find out the name of the table by executing:
SELECT owner, object_name INTO schema_name, table_name
FROM sys.dba_objects WHERE DATA_OBJECT_ID = ' || dobject_id);

--          raise_application_error(-20101,
--                  'Please, create domain index on the column referenced in the operator '
--                  || func_name || ' of the table ' || schema_name
--                  || '.' || table_name);

        -- ELSE
        --   the first argument is a literal value instead of a column
        END IF;
    END IF;
  END; -- PROCEDURE checkTableVersionEx


  FUNCTION getEnvironment RETURN VARCHAR2 IS
    cursor cur is select banner from sys.v_$version;
    res varchar2(4000);
    remote_env varchar2(4000);
    errm varchar2(4000);
  BEGIN
    check_master_table();

    res := 'Oracle environment: ' || chr(10);
    for t in cur loop
        res := res || t.banner || chr(10);
        exit when cur%notfound;
    end loop;

    res := res || chr(10) || 'JChem owner: ' || get_cartowner_schema() || chr(10);

    res := res || chr(10) || 'JChem Server environment: ' || chr(10);

    errm := get_remote_environment(remote_env);
    if errm is null then
      return res || remote_env;
    else
      jchem_core_pkg.handle_java_error(errm);
    end if;
  END;

  FUNCTION jctf_autocalccts(index_schema VARCHAR2,
                            index_name VARCHAR2)
  RETURN COMP_CHAR_ARRAY PIPELINED AS
    carr CHAR_ARRAY;
  BEGIN
    IF index_schema IS NULL OR index_name IS NULL THEN
      RETURN;
    END IF;

    carr := jchem_core_pkg.autocalccts(index_schema, index_name, null);
    FOR i IN 1..carr.count LOOP
      PIPE ROW(COMPOSITE_CHAR(carr(i)));
    END LOOP;

    RETURN;
  END;

  FUNCTION jctf_autocalccts(index_schema VARCHAR2,
                            index_name VARCHAR2,
                            index_part VARCHAR2)
  RETURN COMP_CHAR_ARRAY PIPELINED AS
    carr CHAR_ARRAY;
  BEGIN
    IF index_schema IS NULL OR index_name IS NULL THEN
      RETURN;
    END IF;

    carr := jchem_core_pkg.autocalccts(index_schema, index_name, index_part);
    FOR i IN 1..carr.count LOOP
      PIPE ROW(COMPOSITE_CHAR(carr(i)));
    END LOOP;

    RETURN;
  END;

  FUNCTION jctf_autocalccts_bycol(table_schema VARCHAR2,
                                  table_name VARCHAR2,
                                  table_col VARCHAR2)
  RETURN COMP_CHAR_ARRAY PIPELINED AS
    carr CHAR_ARRAY;
  BEGIN
    IF table_schema IS NULL OR table_name IS NULL THEN
      RETURN;
    END IF;

    carr := jchem_core_pkg.autocalccts_bycol(table_schema, table_name, table_col);
    FOR i IN 1..carr.count LOOP
      PIPE ROW(COMPOSITE_CHAR(carr(i)));
    END LOOP;

    RETURN;
  END;

  PROCEDURE invert_result(res IN OUT NOCOPY CD_ID_ARRAY, tableName VARCHAR2) IS
    inv_array CD_ID_ARRAY := CD_ID_ARRAY(); 
    cdid INTEGER;
    nrows INTEGER;
    cnum INTEGER;
    done BOOLEAN := FALSE;
    equals BOOLEAN := FALSE;
  BEGIN 
    cnum := dbms_sql.open_cursor;
    dbms_sql.parse(cnum, 'select cd_id from ' || tableName,
		   dbms_sql.native);
    dbms_sql.define_column(cnum, 1, cdid);
    nrows := dbms_sql.execute(cnum);
    WHILE NOT done LOOP
        IF dbms_sql.fetch_rows(cnum) > 0 THEN
	    dbms_sql.column_value(cnum, 1, cdid);
	    FOR i in 1 .. res.count LOOP
		IF res(i) = cdid THEN
		    equals := TRUE;
		END IF;
   	    END LOOP;
	    IF equals = FALSE THEN
		inv_array.extend;
		inv_array(inv_array.count) := cdid;
	    END IF;
	    equals := FALSE;
	ELSE 
	    done := TRUE;
	END IF;
    END LOOP;
    dbms_sql.close_cursor(cnum);
    res := inv_array;
  END;

  PROCEDURE insert_into_idxtbl(idxtbl_name IN OUT NOCOPY VARCHAR2,
                               cdid_seq_name IN OUT NOCOPY VARCHAR2,
                               molprops IN OUT NOCOPY MOLPROPS_ARRAY_ARRAY) AS
    cn INTEGER;
  BEGIN
    cn := create_parse_insidxtbl_cur(idxtbl_name, cdid_seq_name, molprops(1));

    BEGIN
      FOR idx in 1 .. molprops.count LOOP
        insert_into_idxtbl_single(idxtbl_name, cdid_seq_name,
                                  molprops(idx), cn);
        molprops(idx).delete;
      END LOOP;

      molprops.delete;

      DBMS_SQL.CLOSE_CURSOR(cn);
    EXCEPTION
    WHEN OTHERS THEN
      IF DBMS_SQL.IS_OPEN(cn) THEN
        DBMS_SQL.CLOSE_CURSOR(cn);
      END IF;
      RAISE;
    END;
  END;

  PROCEDURE insert_into_idxtbl_single(idxtbl_name IN OUT NOCOPY VARCHAR2,
                                     cdid_seq_name IN OUT NOCOPY VARCHAR2,
                                     molprops IN OUT NOCOPY MOLPROPS_ARRAY,
                                     cn INTEGER)
  AS
    ignore NUMBER;
  BEGIN
    FOR idx in 1 .. molprops.count LOOP
      DBMS_SQL.BIND_VARIABLE(cn, ':var' || to_char(idx), molprops(idx));
    END LOOP;

    ignore := DBMS_SQL.EXECUTE(cn);
  END;


  FUNCTION create_parse_insidxtbl_cur(idxtbl_name IN OUT NOCOPY VARCHAR2,
                                      cdid_seq_name IN OUT NOCOPY VARCHAR2,
                                      molprops IN OUT NOCOPY MOLPROPS_ARRAY)
  RETURN INTEGER
  AS
    cn INTEGER;
    insert_stmt VARCHAR2(2000);
  BEGIN
    insert_stmt := 'INSERT INTO ' || idxtbl_name || ' VALUES(' ||
                      cdid_seq_name || '.nextval';

    FOR idx in 1 .. molprops.count LOOP
      insert_stmt := insert_stmt || ', :var' || to_char(idx);
    END LOOP;
    insert_stmt := insert_stmt || ')';

--    jchem_misc_pkg.trace('create_parse_insidxtbl_cur: stmt=' ||
--                         insert_stmt);

    cn := DBMS_SQL.OPEN_CURSOR;
    DBMS_SQL.PARSE(cn, insert_stmt, DBMS_SQL.NATIVE);
    RETURN cn;
  END;

END jchem_core_pkg;
/
SHOW ERRORS;

-- Defines procedures/functions with definer rights
create or replace package jchem_defright_pkg as
  function get_next_scanId return number;
end jchem_defright_pkg;
/
show errors;

create or replace package body jchem_defright_pkg as
  function get_next_scanId return number as
    sqname varchar2(50) := 'jchem_idxscan_no_sq';
    scanId number;
  begin
    execute immediate 'SELECT ' || sqname || '.nextval FROM dual' into scanId;
    return scanId;
  end;
end jchem_defright_pkg;
/
show errors


CREATE OR REPLACE PACKAGE BODY jchem_misc_pkg AS

  FUNCTION get_jchemstreams_url RETURN VARCHAR2 IS
  BEGIN
    return 'JCHEM_SERVLET_URL';
  END;

  PROCEDURE put_any_line (pstrText IN VARCHAR2) IS
    cnumLineSize  CONSTANT NUMBER := 255;   -- Maximum size of DBMS_OUTPUT 
    lstrLeft  VARCHAR2 (32767);
    lstrRight VARCHAR2 (32767);
    lnumPos   NUMBER;
  BEGIN
    IF LENGTH(pstrText) <= cnumLineSize THEN
         -- Line short enough to go to PUT_LINE
      DBMS_OUTPUT.PUT_LINE(pstrText);
    ELSE
         -- Line is too long, attempt to split it at a space
      lstrLeft := SUBSTR(pstrText,1,cnumLineSize);
      lnumPos := INSTR(lstrLeft,' ',-1);
      IF lnumPos = 0 THEN
            -- No spaces in the line, so just split it regardless
        lstrRight := SUBSTR(pstrText,cnumLineSize+1);
      ELSE
        lstrLeft := SUBSTR(lstrLeft,1,lnumPos-1);
        lstrRight := SUBSTR(pstrText,lnumPos+1);
      END IF;
        -- lstrLeft is small enough to send to put_line...
        DBMS_OUTPUT.PUT_LINE (lstrLeft);
        -- ...but lstrRight may not be. Make a recursive call to deal with it
        put_any_line (lstrRight);
      END IF;
   END;

  FUNCTION chartoraw(v_char varchar2) RETURN LONG RAW IS
    rawdata long raw;
    rawlen number;
    hex varchar2(32767);
    i number;
  BEGIN
    rawlen := length(v_char);
    i := 1;
    WHILE i <= rawlen LOOP
      hex := numtohex(ascii(substrb(v_char,i,1)));
      rawdata := rawdata || HEXTORAW(hex);
      i := i + 1;
    END LOOP;
    RETURN rawdata;
  END;


  FUNCTION rawtochar(v_raw long raw) RETURN VARCHAR2 IS
    rawlen number;
    hex varchar2(32767);
    rawparam varchar2(32767);
    i number;
  BEGIN
    hex := rawtohex(v_raw);
    rawlen := length(hex);
    i := 1;
    WHILE i <= rawlen LOOP
      rawparam := rawparam||CHR(HEXTONUM(substrb(hex,i,2)));
      i := i + 2;
    END LOOP;
    RETURN rawparam;
  END;


  FUNCTION numtohex(v_hex number) return varchar2 IS 
    hex varchar2(4);
    num1 number;
    num2 number;
  BEGIN
    num1 := trunc(v_hex/16);
    num2 := v_hex-(num1*16);

    IF ( num1 >= 0 and num1 <= 9 ) THEN
      hex := hex||to_char(num1);
    END IF; 
    if num1 = 10 then hex := hex||'A'; end if; 
    if num1 = 11 then hex := hex||'B'; end if; 
    if num1 = 12 then hex := hex||'C'; end if; 
    if num1 = 13 then hex := hex||'D'; end if; 
    if num1 = 14 then hex := hex||'E'; end if; 
    if num1 = 15 then hex := hex||'F'; end if; 

    if ( num2 >= 0 and num2 <= 9 ) then 
      hex := hex||to_char(num2);
    end if; 
    if num2 = 10 then hex := hex||'A'; end if; 
    if num2 = 11 then hex := hex||'B'; end if; 
    if num2 = 12 then hex := hex||'C'; end if; 
    if num2 = 13 then hex := hex||'D'; end if; 
    if num2 = 14 then hex := hex||'E'; end if; 
    if num2 = 15 then hex := hex||'F'; end if; 

    return hex;
  end;


  FUNCTION hextonum(v_hex varchar2) return number IS 
    hex varchar2(4);
    num number;
    num1 number;
    num2 number;
  BEGIN
    hex := substrb(v_hex,1,1);

    if ( hex >= '0' and hex <= '9' ) then 
      num1 := to_number(hex);
    end if; 
    if hex = 'A' then num1 := 10; end if; 
    if hex = 'B' then num1 := 11; end if; 
    if hex = 'C' then num1 := 12; end if; 
    if hex = 'D' then num1 := 13; end if; 
    if hex = 'E' then num1 := 14; end if; 
    if hex = 'F' then num1 := 15; end if; 

    hex := substrb(v_hex,2,1);

    if ( hex >= '0' and hex <= '9' ) then 
      num2 := to_number(hex);
    end if; 
    if hex = 'A' then num2 := 10; end if; 
    if hex = 'B' then num2 := 11; end if; 
    if hex = 'C' then num2 := 12; end if; 
    if hex = 'D' then num2 := 13; end if; 
    if hex = 'E' then num2 := 14; end if; 
    if hex = 'F' then num2 := 15; end if; 

   num := (num1*16)+num2;
   return num;
  END;

END jchem_misc_pkg;
/
show errors;


CREATE OR REPLACE PACKAGE jchem_func_pkg AUTHID CURRENT_USER AS

  function format_hits(target VARCHAR2, query VARCHAR2, options VARCHAR2,
                       tmpBlob BLOB)
        RETURN BLOB DETERMINISTIC PARALLEL_ENABLE  AS LANGUAGE JAVA NAME
        'chemaxon.jchem.cartridge.orastub.JcDisplayFunctions.formatHits(
                   java.lang.String, java.lang.String, java.lang.String,
                   oracle.sql.BLOB) return oracle.sql.BLOB';

  function format_hits(target BLOB, query BLOB, options VARCHAR2, tmpBlob BLOB)
        RETURN BLOB DETERMINISTIC PARALLEL_ENABLE AS LANGUAGE JAVA NAME
        'chemaxon.jchem.cartridge.orastub.JcDisplayFunctions.formatHits(
                   oracle.sql.BLOB, oracle.sql.BLOB, java.lang.String,
                   oracle.sql.BLOB) return oracle.sql.BLOB';
end jchem_func_pkg;
/
show errors;

create or replace package jchem_table_pkg authid current_user as

  procedure create_jctable(jchem_table_name varchar2,
                               jchem_property_table_name varchar2,
                               number_of_ints number,
                               number_of_ones number,
                               number_of_edges number,
                               coldefs varchar2,
                               standardizerConfig varchar2,
                               absoluteStereo number,
                               options varchar2) as language java name
  'chemaxon.jchem.cartridge.oraresident.JcTableFunctions.createJChemTable(java.lang.String,
  java.lang.String, int, int, int, java.lang.String, java.lang.String,
  int, java.lang.String)';
                               
  procedure drop_jctable(jchem_table_name varchar2,
                             jchem_property_table_name varchar2)
  as language java name
  'chemaxon.jchem.cartridge.oraresident.JcTableFunctions.dropJChemTable(java.lang.String,
  java.lang.String)';

  function list_arr(jcproptable varchar2)
    return char_array parallel_enable as language java name 
    'chemaxon.jchem.cartridge.oraresident.JcTableFunctions.listJChemTablesArr(java.lang.String)
    return java.lang.String[]';

  function list_jctables(jcproptable varchar2) return comp_char_array
  pipelined;

  FUNCTION insert_structure(str VARCHAR2, table_name VARCHAR2, jcprop_name VARCHAR2,
          dup_chk VARCHAR2, halt_on_err VARCHAR2, options VARCHAR2)
        RETURN CD_ID_ARRAY PARALLEL_ENABLE   AS LANGUAGE JAVA NAME
	    'chemaxon.jchem.cartridge.oraresident.JCartDml.insertMol(java.lang.String,
				  java.lang.String, java.lang.String,
				  java.lang.String, java.lang.String,
                  java.lang.String) return int[]';

  FUNCTION insert_structure(str CLOB, table_name VARCHAR2, jcprop_name VARCHAR2,
                    dup_chk VARCHAR2, halt_on_err VARCHAR2, options VARCHAR2)
        RETURN CD_ID_ARRAY PARALLEL_ENABLE AS LANGUAGE JAVA NAME
	    'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.insertMol(oracle.sql.CLOB,
            java.lang.String, java.lang.String, java.lang.String,
            java.lang.String, java.lang.String) return int[]';

  FUNCTION insert_structure(str BLOB, table_name VARCHAR2, jcprop_name VARCHAR2,
                    dup_chk VARCHAR2, halt_on_err VARCHAR2, options VARCHAR2)
        RETURN CD_ID_ARRAY PARALLEL_ENABLE AS LANGUAGE JAVA NAME
	    'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.insertMol(oracle.sql.BLOB,
        java.lang.String, java.lang.String, java.lang.String,
        java.lang.String, java.lang.String) return int[]';

  PROCEDURE update_structure(str VARCHAR2, table_name VARCHAR2, id NUMBER, 
	    jcprop_name VARCHAR2, options VARCHAR2)
          PARALLEL_ENABLE   AS LANGUAGE JAVA NAME
	    'chemaxon.jchem.cartridge.oraresident.JCartDml.updateMol(java.lang.String,
                  java.lang.String, int, java.lang.String, java.lang.String)';

  PROCEDURE update_structure(str CLOB, table_name VARCHAR2, id NUMBER, 
	    jcprop_name VARCHAR2, options VARCHAR2)
        PARALLEL_ENABLE AS LANGUAGE JAVA NAME
	    'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.updateMol(
        oracle.sql.CLOB, java.lang.String, int, java.lang.String,
        java.lang.String)';

  PROCEDURE update_structure(str BLOB, table_name VARCHAR2, id NUMBER, 
	    jcprop_name VARCHAR2, options VARCHAR2)
        PARALLEL_ENABLE AS LANGUAGE JAVA NAME
	    'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.updateMol(
        oracle.sql.BLOB, java.lang.String, int, java.lang.String,
        java.lang.String)';

  PROCEDURE delete_structure(table_name VARCHAR2, condition VARCHAR2,
                       jcprop_name VARCHAR2) AS LANGUAGE JAVA NAME
	    'chemaxon.jchem.cartridge.oraresident.JCartDml.deleteMol(java.lang.String,
				  java.lang.String, java.lang.String)';

  function jc_insert(str varchar2, table_name varchar2, 
                     jcprop_name varchar2 := '', duplicate_check varchar2 :=
                     '', halt_on_error varchar2 := '', options varchar2 := '')
                     return CD_ID_ARRAY;

  function jc_insert(str CLOB, table_name VARCHAR2,
                      jcprop_name VARCHAR2 := '', duplicate_check VARCHAR2 :=
                      '', halt_on_error VARCHAR2 := '', options VARCHAR2 := '')
        return cd_id_array;

  function jc_insert(str BLOB, table_name VARCHAR2,
                      jcprop_name VARCHAR2 := '', duplicate_check VARCHAR2 :=
                      '', halt_on_error VARCHAR2 := '', options VARCHAR2 := '')
        return cd_id_array;

  procedure jc_update(str varchar2, table_name varchar2,
                      id NUMBER, jcprop_name varchar2 := null,
                      options varchar2 := null);
  
  procedure jc_update(str CLOB, table_name VARCHAR2,
                           id NUMBER, jcprop_name VARCHAR2 := null,
                           options VARCHAR2 := null);

  procedure jc_update(str BLOB, table_name VARCHAR2,
                           id NUMBER, jcprop_name VARCHAR2 := null,
                           options VARCHAR2 := null);

  procedure jc_delete(table_name varchar2, condition varchar2,
                      jcprop_name varchar2 := null);

  function jc_insertb(str BLOB, table_name VARCHAR2,
                      jcprop_name VARCHAR2 := '', duplicate_check VARCHAR2 :=
                      '', halt_on_error VARCHAR2 := '', options VARCHAR2 := '')
        return cd_id_array;

  procedure jc_updateb(str BLOB, table_name VARCHAR2,
                           id NUMBER, jcprop_name VARCHAR2 := null,
                           options VARCHAR2 := null);
end jchem_table_pkg;
/
show errors;

create or replace package body jchem_table_pkg as
  function list_jctables(jcproptable varchar2) return comp_char_array pipelined
  as
    chararr char_array;
  begin
    chararr := list_arr(jcproptable);
    for i in 1..chararr.count loop
      pipe row(COMPOSITE_CHAR(chararr(i)));
    end loop;

    return;
  end;

  function jc_insert(str varchar2, table_name varchar2, 
                     jcprop_name varchar2 := '', duplicate_check varchar2 :=
                     '', halt_on_error varchar2 := '', options varchar2 := '')
                     return CD_ID_ARRAY as
    cdidarr cd_id_array;
  begin
    cdidarr := jchem_table_pkg.insert_structure(str, table_name,
        jcprop_name, duplicate_check, halt_on_error, options);
    if (cdidarr is null) then
      jchem_core_pkg.handle_java_error(jchem_core_pkg.lastError);
    else
      return cdidarr;
    end if;
  end;

  function jc_insert(str CLOB, table_name VARCHAR2,
                      jcprop_name VARCHAR2 := '', duplicate_check VARCHAR2 :=
                      '', halt_on_error VARCHAR2 := '', options VARCHAR2 := '')
        return cd_id_array as
    cdidarr cd_id_array;
  begin
    cdidarr := jchem_table_pkg.insert_structure(str, table_name,
        jcprop_name, duplicate_check, halt_on_error, options);
    if (cdidarr is null) then
      jchem_core_pkg.handle_java_error(jchem_core_pkg.lastError);
    else
      return cdidarr;
    end if;
  end;

  function jc_insert(str BLOB, table_name VARCHAR2,
                      jcprop_name VARCHAR2 := '', duplicate_check VARCHAR2 :=
                      '', halt_on_error VARCHAR2 := '', options VARCHAR2 := '')
        return cd_id_array as
    cdidarr cd_id_array;
  begin
    cdidarr := jchem_table_pkg.insert_structure(str, table_name,
        jcprop_name, duplicate_check, halt_on_error, options);
    if (cdidarr is null) then
      jchem_core_pkg.handle_java_error(jchem_core_pkg.lastError);
    else
      return cdidarr;
    end if;
  end;

  procedure jc_update(str varchar2, table_name varchar2,
                      id NUMBER, jcprop_name varchar2 := null,
                      options varchar2 := null) as
  begin
    jchem_table_pkg.update_structure(str, table_name, id, jcprop_name, options);
  end;
  
  procedure jc_update(str CLOB, table_name VARCHAR2,
                           id NUMBER, jcprop_name VARCHAR2 := null,
                           options VARCHAR2 := null) AS
  begin
    jchem_table_pkg.update_structure(str, table_name, id, jcprop_name, options);
  end;

  procedure jc_update(str BLOB, table_name VARCHAR2,
                           id NUMBER, jcprop_name VARCHAR2 := null,
                           options VARCHAR2 := null) AS
  begin
    jchem_table_pkg.update_structure(str, table_name, id, jcprop_name, options);
  end;

  procedure jc_delete(table_name varchar2, condition varchar2,
                      jcprop_name varchar2 := null) as
  begin
    jchem_table_pkg.delete_structure(table_name, condition, jcprop_name);
  end;

  function jc_insertb(str BLOB, table_name VARCHAR2,
                      jcprop_name VARCHAR2 := '', duplicate_check VARCHAR2 :=
                      '', halt_on_error VARCHAR2 := '', options VARCHAR2 := '')
        return cd_id_array as
    cdidarr cd_id_array;
  begin
      cdidarr := jchem_table_pkg.insert_structure(str, table_name,
          jcprop_name, duplicate_check, halt_on_error, options);
      if (cdidarr is null) then
        jchem_core_pkg.handle_java_error(jchem_core_pkg.lastError);
      else
        return cdidarr;
      end if;
  end;

  procedure jc_updateb(str BLOB, table_name VARCHAR2,
                           id NUMBER, jcprop_name VARCHAR2 := null,
                           options VARCHAR2 := null) AS
  begin
    jchem_table_pkg.update_structure(str, table_name, id, jcprop_name, options);
  end;

end jchem_table_pkg;
/
show errors;

---------------------------------
--  CREATE IMPLEMENTATION UNIT --
---------------------------------

-- CREATE TYPE BODY
CREATE OR REPLACE TYPE BODY jc_idxtype_im IS

   STATIC FUNCTION ODCIGetInterfaces(ifclist OUT sys.ODCIObjectList) 
							RETURN NUMBER IS
   BEGIN
       ifclist := sys.ODCIObjectList(sys.ODCIObject('SYS','ODCIINDEX2'));
       RETURN ODCIConst.Success;
   END ODCIGetInterfaces;

  STATIC FUNCTION ODCIIndexCreate (ia sys.odciindexinfo, 
                                   parms VARCHAR2,
                                   env sys.ODCIEnv) RETURN NUMBER IS
    errmsg varchar2(32767);
	endstatus number;
	rindex binary_integer;
	slno binary_integer;
  BEGIN
      IF ( (env.CallProperty IS NULL) OR
           (env.CallProperty = SYS.ODCIConst.IntermediateCall) ) THEN
        errmsg := jchem_core_pkg.indexCreate(ia.IndexSchema, 
            ia.IndexName, ia.IndexPartition, ia.IndexCols(1).TableSchema, 
            ia.IndexCols(1).TableName, ia.IndexCols(1).TablePartition,
            ia.IndexCols(1).colName, ia.IndexCols(1).colTypeName, parms);
      END IF;

      if errmsg is null then
        RETURN ODCICONST.SUCCESS;
      else
        jchem_core_pkg.handle_java_error(errmsg);
      end if;

      RETURN ODCIConst.SUCCESS;
  END;

  STATIC FUNCTION ODCIIndexAlter(ia SYS.ODCIIndexInfo, parms IN OUT VARCHAR2, 
                                 alter_option NUMBER, env sys.ODCIEnv)
  RETURN NUMBER IS
    ret NUMBER;
  BEGIN
    /*
    begin
    */
      ret := jchem_core_pkg.indexAlter(ia.IndexSchema, 
          ia.IndexName, ia.IndexPartition, ia.IndexCols(1).TableSchema,
          ia.IndexCols(1).TableName, ia.IndexCols(1).TablePartition,
          ia.IndexCols(1).colName, ia.IndexCols(1).colTypeName, parms, alter_option);
      RETURN ret;
    /*
    exception
    when others then
      if sqlcode = -29532 then
        jchem_core_pkg.handle_java_error(sqlerrm);
      else
        raise;
      end if;
    end;
    */
  END;

  STATIC FUNCTION ODCIIndexDrop(ia sys.odciindexinfo, env sys.ODCIEnv)
    RETURN NUMBER IS
    junk PLS_INTEGER;
  BEGIN
    /*
    begin
    */

      /*
      DEBUG 
      jchem_misc_pkg.trace('ODCIIndexDrop: BEGIN');
      if env.CallProperty IS NULL then
        jchem_misc_pkg.trace('ODCIIndexDrop: env.CallProperty is null');
      else
        jchem_misc_pkg.trace('ODCIIndexDrop: env.CallProperty is ' || to_char(env.CallProperty));
      end if;
      */

      IF ( (env.CallProperty IS NULL) OR
           (env.CallProperty = SYS.ODCIConst.IntermediateCall) ) THEN
        RETURN jchem_core_pkg.indexDrop(ia.IndexSchema, ia.IndexName,
                  ia.IndexPartition, ia.IndexCols(1).TableSchema,
                  ia.IndexCols(1).TableName);
      END IF;
      RETURN ODCIConst.Success;
    /*
    exception
    when others then
      if sqlcode = -29532 then
        jchem_core_pkg.handle_java_error(sqlerrm);
      else
        raise;
      end if;
    end;
    */
  END;

  STATIC FUNCTION ODCIIndexTruncate(ia sys.odciindexinfo, env sys.ODCIEnv)
    RETURN NUMBER IS
    junk PLS_INTEGER;
  BEGIN
    /*
    begin
    */
      IF ( (env.CallProperty IS NULL) OR
           (env.CallProperty = SYS.ODCIConst.IntermediateCall) ) THEN
        RETURN jchem_core_pkg.indexTruncate(ia.IndexSchema, ia.IndexName,
                  ia.IndexPartition, ia.IndexCols(1).TableSchema,
                  ia.IndexCols(1).TableName);
      END IF;
      RETURN ODCIConst.Success;
    /*
    exception
    when others then
      if sqlcode = -29532 then
        jchem_core_pkg.handle_java_error(sqlerrm);
      else
        raise;
      end if;
    end;
    */
  END;

  STATIC FUNCTION ODCIIndexStart(sctx IN OUT jc_idxtype_im, 
                        ia sys.odciindexinfo,
                        op sys.odciPredInfo, 
                        qi sys.ODCIQueryInfo, 
                        strt NUMBER, 
                        stop NUMBER,
                        env sys.ODCIEnv) RETURN NUMBER IS
    scanId number;
    cnum INTEGER;
    stmt VARCHAR2(32000);
    rid ROWID;
    nrows INTEGER;

    errmsg varchar2(32767);
  BEGIN
      jchem_core_pkg.checkTableVersion(
                ia.IndexSchema,
                ia.IndexName,
                ia.IndexPartition,
                ia.IndexCols(1).TableSchema,
                ia.IndexCols(1).TableName,
                ia.IndexCols(1).ColName,
				ia.IndexCols(1).ColTypeName);

      scanId := jchem_defright_pkg.get_next_scanId();

      IF op.ObjectName = 'JC_MOLWEIGHT' OR op.ObjectName = 'JC_MOLWEIGHTB' THEN
        stmt := jchem_core_pkg.getSqlForMolweightScan(ia.IndexSchema,
                      ia.IndexName, ia.IndexPartition, TO_CHAR(strt), TO_CHAR(stop), op.Flags);
        cnum := dbms_sql.open_cursor;
        dbms_sql.parse(cnum, stmt, dbms_sql.native);
        dbms_sql.define_column_rowid(cnum, 1, rid);
        nrows := dbms_sql.execute(cnum);
        sctx := jc_idxtype_im(ia, scanId, cnum);
        return ODCICONST.SUCCESS;
      ELSE
          errmsg := jchem_clob_pkg.index_scan(ia.IndexSchema, ia.IndexName,
                          ia.IndexPartition, ia.IndexCols(1).TableSchema,
                          ia.IndexCols(1).TableName, ia.IndexCols(1).colName, ia.IndexCols(1).colTypeName,
                          op.ObjectName, 'userdef', strt, stop, op.Flags, null,
                          op.ObjectName, scanId);
          -- dbms_output.put_line('user def batch');
          --   raise_application_error(-20101, 'Unsupported operator: '
          --              || op.ObjectName);
      END IF;

      -- jchem_misc_pkg.new_profile_point();
      sctx := jc_idxtype_im(ia, scanId, null);

      if errmsg is null then
        RETURN ODCICONST.SUCCESS;
      else
        jchem_core_pkg.handle_java_error(errmsg);
      end if;

      RETURN ODCICONST.SUCCESS;
  END;  

  STATIC FUNCTION ODCIIndexStart(sctx IN OUT NOCOPY jc_idxtype_im, 
          ia sys.odciindexinfo,
          op sys.odciPredInfo, 
          qi sys.ODCIQueryInfo, 
          strt NUMBER, 
          stop NUMBER,
          query CLOB,
          env sys.ODCIEnv) RETURN NUMBER IS
    scanId number;
    cnum INTEGER;
    stmt varchar2(32000);
    rid ROWID;
    nrows INTEGER;
    regCode VARCHAR2(20);
    errmsg varchar2(32767);
  BEGIN

--    jchem_core_pkg.profile_array := TIMESTAMP_ARRAY();

      jchem_core_pkg.checkTableVersion(
                ia.IndexSchema,
                ia.IndexName,ia.IndexPartition,
                ia.IndexCols(1).TableSchema,
                ia.IndexCols(1).TableName,
                ia.IndexCols(1).ColName,
				ia.IndexCols(1).ColTypeName);

      scanId := jchem_defright_pkg.get_next_scanId();

      IF op.ObjectName = 'JC_CONTAINS' or op.ObjectName = 'JC_EQUALS' or op.ObjectName = 'JC_MATCHCOUNT'  or
	     op.ObjectName= 'JC_EVALUATE' THEN 

          errmsg := jchem_clob_pkg.index_scan(ia.IndexSchema, ia.IndexName,
                          ia.IndexPartition, ia.IndexCols(1).TableSchema,
                          ia.IndexCols(1).TableName, ia.IndexCols(1).colName, ia.IndexCols(1).colTypeName,
                          op.ObjectName, null, strt, stop, op.Flags, query,
                          null, scanId);

      ELSIF op.ObjectName = 'JC_TANIMOTO' OR op.ObjectName = 'JC_DISSIMILARITY' THEN

        errmsg := jchem_clob_pkg.similaritySearch(query, op.ObjectName,
                          strt, stop, op.Flags, 'dissimilarityMetric:tanimoto',
                          ia.IndexSchema, ia.IndexName,ia.IndexPartition,
                          ia.IndexCols(1).TableSchema,
                          ia.IndexCols(1).TableName, ia.IndexCols(1).colName, ia.IndexCols(1).colTypeName,
                          scanId);

      ELSIF op.ObjectName = 'JC_FORMULA_EQB' or op.ObjectName = 'JC_FORMULA_EQ' THEN

        stmt := jchem_clob_pkg.getSqlForFormulaScan(ia.IndexSchema,
                          ia.IndexName, ia.IndexPartition, query, strt);

        cnum := dbms_sql.open_cursor; 
        dbms_sql.parse(cnum, stmt, dbms_sql.native);
        dbms_sql.define_column_rowid(cnum, 1, rid);    
        nrows := dbms_sql.execute(cnum);

        sctx := jc_idxtype_im(ia, scanId, cnum);
        return ODCICONST.SUCCESS;
      ELSE
          errmsg := jchem_clob_pkg.index_scan(ia.IndexSchema, ia.IndexName,
                          ia.IndexPartition, ia.IndexCols(1).TableSchema,
                          ia.IndexCols(1).TableName, ia.IndexCols(1).colName, ia.IndexCols(1).colTypeName,
                          op.ObjectName, 'userdef', strt, stop, op.Flags,
                          query, op.ObjectName, scanId);
          --raise_application_error(-20101, 
          --		'Unsupported operator: '|| op.ObjectName);
      END IF;

      sctx := jc_idxtype_im(ia, scanId, null);

  --    jchem_misc_pkg.trace(op.ObjectName || ': BATCH END');

      -- jchem_misc_pkg.new_profile_point();

      if errmsg is null then
        RETURN ODCICONST.SUCCESS;
      else
        jchem_core_pkg.handle_java_error(errmsg);
      end if;

      RETURN ODCICONST.SUCCESS;
  END;
 
  STATIC FUNCTION ODCIIndexStart(sctx IN OUT NOCOPY jc_idxtype_im, 
          ia sys.odciindexinfo,
          op sys.odciPredInfo, 
          qi sys.ODCIQueryInfo, 
          strt NUMBER, 
          stop NUMBER,
          query BLOB,
          env sys.ODCIEnv) RETURN NUMBER IS
    scanId number;
    cnum INTEGER;
    stmt varchar2(32000);
    rid ROWID;
    nrows INTEGER;
    regCode VARCHAR2(20);
    errmsg varchar2(32767);
  BEGIN

--    jchem_core_pkg.profile_array := TIMESTAMP_ARRAY();

      jchem_core_pkg.checkTableVersion(
                ia.IndexSchema,
                ia.IndexName,ia.IndexPartition,
                ia.IndexCols(1).TableSchema,
                ia.IndexCols(1).TableName,
                ia.IndexCols(1).ColName,
				ia.IndexCols(1).ColTypeName);

      scanId := jchem_defright_pkg.get_next_scanId();

      IF op.ObjectName = 'JC_CONTAINSB' or op.ObjectName = 'JC_CONTAINS' or op.ObjectName = 'JC_EQUALSB' or op.ObjectName = 'JC_EQUALS' or
	     op.ObjectName = 'JC_MATCHCOUNTB' or op.ObjectName = 'JC_MATCHCOUNT' or op.ObjectName = 'JC_EVALUATEB' OR op.ObjectName = 'JC_EVALUATE' THEN

          errmsg := jchem_blob_pkg.index_scan(ia.IndexSchema, ia.IndexName,
                          ia.IndexPartition, ia.IndexCols(1).TableSchema,
                          ia.IndexCols(1).TableName, ia.IndexCols(1).colName, ia.IndexCols(1).colTypeName,
                          op.ObjectName, null, strt, stop, op.Flags, query,
                          null, scanId);

      ELSIF op.ObjectName = 'JC_TANIMOTOB' or op.ObjectName = 'JC_TANIMOTO' or op.ObjectName = 'JC_DISSIMILARITYB' or op.ObjectName = 'JC_DISSIMILARITY' THEN

        errmsg := jchem_blob_pkg.similaritySearch(query, op.ObjectName,
                          strt, stop, op.Flags, 'dissimilarityMetric:tanimoto',
                          ia.IndexSchema, ia.IndexName,ia.IndexPartition,
                          ia.IndexCols(1).TableSchema,
                          ia.IndexCols(1).TableName, ia.IndexCols(1).colName, ia.IndexCols(1).colTypeName,
                          scanId);

      ELSIF op.ObjectName = 'JC_FORMULA_EQB' or op.ObjectName = 'JC_FORMULA_EQ' THEN

        stmt := jchem_blob_pkg.getSqlForFormulaScan(
                          ia.IndexSchema, ia.IndexName, ia.IndexPartition,
                          query, strt);

        cnum := dbms_sql.open_cursor; 
        dbms_sql.parse(cnum, stmt, dbms_sql.native);
        dbms_sql.define_column_rowid(cnum, 1, rid);    
        nrows := dbms_sql.execute(cnum);

        sctx := jc_idxtype_im(ia, scanId, cnum);
        return ODCICONST.SUCCESS;
      ELSE
          raise_application_error(-20101, 'Unsupported operator: '
              || op.ObjectName);
      END IF;

      sctx := jc_idxtype_im(ia, scanId, null);

  --    jchem_misc_pkg.trace(op.ObjectName || ': BATCH END');

      if errmsg is null then
        RETURN ODCICONST.SUCCESS;
      else
        jchem_core_pkg.handle_java_error(errmsg);
      end if;

      -- jchem_misc_pkg.new_profile_point();
      RETURN ODCICONST.SUCCESS;
  END;
 
  STATIC FUNCTION ODCIIndexStart(sctx IN OUT jc_idxtype_im, 
                                 ia sys.odciindexinfo,
                                 op sys.odciPredInfo,
                                 qi sys.ODCIQueryInfo,
                                 strt NUMBER,
                                 stop NUMBER,
                                 query CLOB,
                                 param VARCHAR2,
                                 env sys.ODCIEnv) RETURN NUMBER IS
    scanId number;
    errmsg varchar2(32767);
  BEGIN

--    jchem_core_pkg.profile_array := TIMESTAMP_ARRAY();

      jchem_core_pkg.checkTableVersion(
                ia.IndexSchema,
                ia.IndexName, ia.IndexPartition,
                ia.IndexCols(1).TableSchema,
                ia.IndexCols(1).TableName,
                ia.IndexCols(1).ColName,
				ia.IndexCols(1).ColTypeName);

      scanId := jchem_defright_pkg.get_next_scanId();

      IF op.ObjectName = 'JC_COMPARE' or op.ObjectName = 'JC_EVALUATE' THEN

          errmsg := jchem_clob_pkg.index_scan(ia.IndexSchema, ia.IndexName,
                          ia.IndexPartition, ia.IndexCols(1).TableSchema,
                          ia.IndexCols(1).TableName, ia.IndexCols(1).colName, ia.IndexCols(1).colTypeName,
                          op.ObjectName, null, strt, stop, op.Flags, query,
                          param, scanId);

      ELSIF op.ObjectName = 'JC_DISSIMILARITY' THEN

        errmsg := jchem_clob_pkg.similaritySearch(query, op.ObjectName,
                          strt, stop, op.Flags, param, ia.IndexSchema,
                          ia.IndexName,ia.IndexPartition,
                          ia.IndexCols(1).TableSchema,
                          ia.IndexCols(1).TableName, ia.IndexCols(1).colName, ia.IndexCols(1).colTypeName,
                          scanId);
        
      ELSE
          raise_application_error(-20101, 'Unsupported operator: '
              || op.ObjectName);
      END IF;

      sctx := jc_idxtype_im(ia, scanId, null); 
    
      if errmsg is null then
        RETURN ODCICONST.SUCCESS;
      else
        jchem_core_pkg.handle_java_error(errmsg);
      end if;

      -- jchem_misc_pkg.new_profile_point();
      RETURN ODCICONST.SUCCESS;
  END;


  STATIC FUNCTION ODCIIndexStart(sctx IN OUT jc_idxtype_im, 
                                 ia sys.odciindexinfo,
                                 op sys.odciPredInfo,
                                 qi sys.ODCIQueryInfo,
                                 strt NUMBER,
                                 stop NUMBER,
                                 query BLOB,
                                 param VARCHAR2,
                                 env sys.ODCIEnv) RETURN NUMBER IS
    scanId number;
    errmsg varchar2(32767);
  BEGIN

--    jchem_core_pkg.profile_array := TIMESTAMP_ARRAY();

      jchem_core_pkg.checkTableVersion(
                ia.IndexSchema,
                ia.IndexName, ia.IndexPartition,
                ia.IndexCols(1).TableSchema,
                ia.IndexCols(1).TableName,
                ia.IndexCols(1).ColName,
				ia.IndexCols(1).ColTypeName);

      scanId := jchem_defright_pkg.get_next_scanId();

      IF op.ObjectName = 'JC_COMPAREB' or op.ObjectName = 'JC_COMPARE_VB' or op.ObjectName = 'JC_COMPARE' or op.ObjectName = 'JC_EVALUATE' THEN

          errmsg := jchem_blob_pkg.index_scan(ia.IndexSchema, ia.IndexName,
                          ia.IndexPartition, ia.IndexCols(1).TableSchema,
                          ia.IndexCols(1).TableName, ia.IndexCols(1).colName, ia.IndexCols(1).colTypeName,
                          op.ObjectName, null, strt, stop, op.Flags, query,
                          param, scanId);

      ELSIF op.ObjectName = 'JC_DISSIMILARITY' THEN

        errmsg := jchem_blob_pkg.similaritySearch(query, op.ObjectName,
                          strt, stop, op.Flags, param, ia.IndexSchema,
                          ia.IndexName,ia.IndexPartition,
                          ia.IndexCols(1).TableSchema,
                          ia.IndexCols(1).TableName, ia.IndexCols(1).colName, ia.IndexCols(1).colTypeName,
                          scanId);
        
      END IF;

      sctx := jc_idxtype_im(ia, scanId, null); 
    
      if errmsg is null then
        RETURN ODCICONST.SUCCESS;
      else
        jchem_core_pkg.handle_java_error(errmsg);
      end if;

      -- jchem_misc_pkg.new_profile_point();
      RETURN ODCICONST.SUCCESS;
  END;

  STATIC FUNCTION ODCIIndexStart(sctx IN OUT jc_idxtype_im, 
                                 ia sys.odciindexinfo,
                                 op sys.odciPredInfo,
                                 qi sys.ODCIQueryInfo,
                                 strt NUMBER,
                                 stop NUMBER,
                                 query CLOB,
                                 param1 number, param2 number,
                                 env sys.ODCIEnv) RETURN NUMBER IS
    scanId number;
    errmsg varchar2(32767);
  BEGIN

--    jchem_core_pkg.profile_array := TIMESTAMP_ARRAY();

      jchem_core_pkg.checkTableVersion(
                ia.IndexSchema,
                ia.IndexName, ia.IndexPartition,
                ia.IndexCols(1).TableSchema,
                ia.IndexCols(1).TableName,
                ia.IndexCols(1).ColName,
				ia.IndexCols(1).ColTypeName);

      scanId := jchem_defright_pkg.get_next_scanId();

      IF op.ObjectName = 'JC_TVERSKY' THEN
          errmsg := jchem_clob_pkg.similaritySearch(query,
                          op.ObjectName,strt, stop, op.Flags,
                          'dissimilarityMetric:tversky' ||
                          jchem_core_pkg.simCalcSeparator || param2 ||
                          jchem_core_pkg.simCalcSeparator || param1,
                          ia.IndexSchema, ia.IndexName,ia.IndexPartition,
                          ia.IndexCols(1).TableSchema,
                          ia.IndexCols(1).TableName, ia.IndexCols(1).colName, ia.IndexCols(1).colTypeName,
                          scanId);
      ELSE
          raise_application_error(-20101, 'Unsupported operator: '
              || op.ObjectName);
      END IF;

      sctx := jc_idxtype_im(ia, scanId, null); 
    
      if errmsg is null then
        RETURN ODCICONST.SUCCESS;
      else
        jchem_core_pkg.handle_java_error(errmsg);
      end if;

      -- jchem_misc_pkg.new_profile_point();
      RETURN ODCICONST.SUCCESS;

  END;


  STATIC FUNCTION ODCIIndexStart(sctx IN OUT jc_idxtype_im, 
                                 ia sys.odciindexinfo,
                                 op sys.odciPredInfo,
                                 qi sys.ODCIQueryInfo,
                                 strt NUMBER,
                                 stop NUMBER,
                                 query BLOB,
                                 param1 number, param2 number,
                                 env sys.ODCIEnv) RETURN NUMBER IS
    scanId number;
    errmsg varchar2(32767);
  BEGIN

--    jchem_core_pkg.profile_array := TIMESTAMP_ARRAY();

      jchem_core_pkg.checkTableVersion(
                ia.IndexSchema,
                ia.IndexName, ia.IndexPartition,
                ia.IndexCols(1).TableSchema,
                ia.IndexCols(1).TableName,
                ia.IndexCols(1).ColName,
				ia.IndexCols(1).ColTypeName);

      scanId := jchem_defright_pkg.get_next_scanId();

      IF op.ObjectName = 'JC_TVERSKY' THEN
          errmsg := jchem_blob_pkg.similaritySearch(query, op.ObjectName,
                          strt, stop, op.Flags, 'dissimilarityMetric=tversky'
                          || jchem_core_pkg.simCalcSeparator || param2 ||
                          jchem_core_pkg.simCalcSeparator || param1,
                          ia.IndexSchema, ia.IndexName,ia.IndexPartition,
                          ia.IndexCols(1).TableSchema,
                          ia.IndexCols(1).TableName, ia.IndexCols(1).colName, ia.IndexCols(1).colTypeName,
                          scanId);
      ELSE
          raise_application_error(-20101, 'Unsupported operator: '
              || op.ObjectName);
      END IF;

      sctx := jc_idxtype_im(ia, scanId, null); 
    
      if errmsg is null then
        RETURN ODCICONST.SUCCESS;
      else
        jchem_core_pkg.handle_java_error(errmsg);
      end if;

      -- jchem_misc_pkg.new_profile_point();
      RETURN ODCICONST.SUCCESS;

  END;



  MEMBER FUNCTION ODCIIndexFetch(self IN OUT NOCOPY jc_idxtype_im, nrows NUMBER,
                                 rids OUT SYS.ODCIRidList, env sys.ODCIEnv)
  RETURN NUMBER IS
    l_idx_out    INTEGER := 1;
	l_rowids	 RIDARRAY;
    l_rlist      sys.ODCIRidList := sys.ODCIRidList();
    l_refcur     jchem_refcur_pkg.refcur_t;
	type t_rid IS TABLE OF ROWID;
    l_rid        t_rid;
    l_s          VARCHAR2(4000);

    l_count      integer := 1;
    l_done       boolean := false;

    l_nrRemainingRowids integer;

    invalid_rowid exception;
    pragma exception_init(invalid_rowid, -1410);

    errmsg VARCHAR2(32767);
  BEGIN
--    jchem_misc_pkg.trace('scanId=' || to_char(self.scan_num));

      if cnum is not null then
        WHILE not l_done LOOP
          IF l_count > nrows THEN
            l_done := TRUE;
          ELSE
            l_rlist.extend;
            IF dbms_sql.fetch_rows(cnum) > 0 THEN
              dbms_sql.column_value_rowid(cnum, 1, l_rlist(l_count));
              l_count := l_count + 1;
            ELSE
              l_rlist(l_count) := null;
              l_done := TRUE;
              dbms_sql.close_cursor(cnum);
              cnum := null;
            END IF;
          END IF;
        END LOOP;
      else
        errmsg := jchem_core_pkg.nr_remaining_rowids(self.scan_num, l_nrRemainingRowids);
        jchem_core_pkg.handle_java_error(errmsg);
        if l_nrRemainingRowids = 0 then
          l_rlist.extend;
          l_rlist(l_rlist.count) := null;
        else
		  if jchem_core_pkg.usearrays then 
            errmsg := jchem_core_pkg.get_rowids_array(self.scan_num, nrows, l_rowids);
		  else 
            errmsg := jchem_core_pkg.get_rowids_resultset(self.scan_num, nrows, l_refcur);
		  end if;
          jchem_core_pkg.handle_java_error(errmsg);
		  if jchem_core_pkg.usearrays then
		    SELECT CAST(l_rowids as sys.ODCIRidList) INTO l_rlist FROM dual;
		  else
            begin
              loop
                fetch l_refcur bulk collect into l_rid;
                if l_rid is not null and l_rid.count > 0 then
    --              jchem_misc_pkg.trace(l_rid);
	  			  l_rlist.extend(l_rid.count);
				  for curr_rid in l_rid.first .. l_rid.last loop
					  l_rlist(l_idx_out) := l_rid(curr_rid);
					  l_idx_out := l_idx_out + 1;
				  end loop;
                else
                  exit;
                end if;
                exit when l_refcur%notfound;
              end loop;
            exception
              when invalid_rowid then
                fetch l_refcur into l_s;
                raise_application_error(-20101, 
                    'Fake ROWID: '''|| l_s || ''', l_idx_out=' ||
                    to_char(l_idx_out));
			end;
          end if;		  		  
          errmsg := jchem_core_pkg.nr_remaining_rowids(self.scan_num, l_nrRemainingRowids);
          jchem_core_pkg.handle_java_error(errmsg);
          if l_nrRemainingRowids = 0 then
            l_rlist.extend;
            l_rlist(l_rlist.count) := null;
          end if;
		  if NOT jchem_core_pkg.usearrays then 
            close l_refcur;
		  end if;
        end if;
      end if;

      rids := l_rlist;
      RETURN ODCICONST.SUCCESS;
  END;

  MEMBER FUNCTION ODCIIndexClose(self IN OUT jc_idxtype_im, 
					env sys.ODCIEnv) RETURN NUMBER IS
    errmsg varchar2(32767);
  BEGIN
--    jchem_misc_pkg.trace('ODCIIndexClose BEGIN');
    if cnum is not null and dbms_sql.is_open(cnum) then
      dbms_sql.close_cursor(cnum);
      cnum := null;
    end if;
    errmsg := jchem_core_pkg.close_scan_resultset(self.scan_num);
    jchem_core_pkg.handle_java_error(errmsg);
    -- jchem_misc_pkg.new_profile_point();
    RETURN ODCICONST.SUCCESS;
  END;

  STATIC FUNCTION ODCIIndexInsert(ia sys.odciindexinfo, rid VARCHAR2, 
                          newval VARCHAR2, env sys.ODCIEnv) RETURN NUMBER IS 
    errmsg varchar2(32767);
  BEGIN
--    jchem_core_pkg.checkTableVersion(
--              ia.IndexSchema,
--              ia.IndexName,
--              ia.IndexCols(1).TableSchema,
--              ia.IndexCols(1).TableName,
--              ia.IndexCols(1).ColName,
--				ia.IndexCols(1).ColTypeName);

      errmsg := jchem_core_pkg.insert_mol_into_idxtable(newval,
            ia.IndexSchema, ia.IndexName, ia.IndexPartition,
            ia.IndexCols(1).TableSchema, ia.IndexCols(1).TableName, 
            ia.IndexCols(1).ColName, ia.IndexCols(1).ColTypeName, 
            rid);

      if errmsg is null then
        RETURN ODCICONST.SUCCESS;
      else
        jchem_core_pkg.handle_java_error(errmsg);
      end if;

      RETURN ODCICONST.SUCCESS;
  END;  

  STATIC FUNCTION ODCIIndexDelete(ia sys.odciindexinfo, rid VARCHAR2, 
			  oldval VARCHAR2, env sys.ODCIEnv) RETURN NUMBER IS 
   BEGIN
    /*
    begin
    */
      jchem_core_pkg.delete_mol_from_idxtable(
          ia.IndexSchema, ia.IndexName, ia.IndexPartition,
                  ia.IndexCols(1).TableSchema, ia.IndexCols(1).TableName, rid);
      RETURN ODCICONST.SUCCESS;
    /*
    exception
    when others then
      if sqlcode = -29532 then
        jchem_core_pkg.handle_java_error(sqlerrm);
      else
        raise;
      end if;
    end;
    */
  END;
  
  STATIC FUNCTION ODCIIndexUpdate(ia sys.odciindexinfo, rid VARCHAR2, 
         oldval VARCHAR2, newval VARCHAR2, env sys.ODCIEnv) RETURN NUMBER IS 
    errmsg varchar2(32767);
  BEGIN
--    jchem_core_pkg.checkTableVersion(
--              ia.IndexSchema,
--              ia.IndexName,
--              ia.IndexCols(1).TableSchema,
--              ia.IndexCols(1).TableName,
--              ia.IndexCols(1).ColName,
--				ia.IndexCols(1).ColTypeName);

      errmsg := jchem_core_pkg.update_mol_idxtable(oldval, newval, 
                      ia.IndexSchema, ia.IndexName, ia.IndexPartition,
                      ia.IndexCols(1).TableSchema, ia.IndexCols(1).TableName,
                      ia.IndexCols(1).ColName, ia.IndexCols(1).ColTypeName, rid);    

      if errmsg is null then
        RETURN ODCICONST.SUCCESS;
      else
        jchem_core_pkg.handle_java_error(errmsg);
      end if;

      RETURN ODCICONST.SUCCESS;
  END;

  STATIC FUNCTION ODCIIndexInsert(ia sys.odciindexinfo, rid VARCHAR2, 
                                  newval CLOB, env sys.ODCIEnv) RETURN NUMBER IS 
    errmsg varchar2(32767);
  BEGIN
  --    jchem_core_pkg.checkTableVersion(
  --              ia.IndexSchema,
  --              ia.IndexName,
  --              ia.IndexCols(1).TableSchema,
  --              ia.IndexCols(1).TableName,
  --              ia.IndexCols(1).ColName,
  --     		  ia.IndexCols(1).ColTypeName);

      errmsg := jchem_clob_pkg.insert_mol_into_idxtable(newval, 
                    ia.IndexSchema, ia.IndexName, ia.IndexPartition,
                    ia.IndexCols(1).TableSchema, ia.IndexCols(1).TableName, 
                    ia.IndexCols(1).ColName, ia.IndexCols(1).ColTypeName, 
                    rid);

      if errmsg is null then
        RETURN ODCICONST.SUCCESS;
      else
        jchem_core_pkg.handle_java_error(errmsg);
      end if;

      RETURN ODCICONST.SUCCESS;
  END;  

  STATIC FUNCTION ODCIIndexDelete(ia sys.odciindexinfo, rid VARCHAR2, 
                                  oldval CLOB, env sys.ODCIEnv) RETURN NUMBER IS 
  BEGIN
    /*
    begin
    */
      jchem_core_pkg.delete_mol_from_idxtable(
              ia.IndexSchema, ia.IndexName, ia.IndexPartition,
              ia.IndexCols(1).TableSchema, ia.IndexCols(1).TableName, rid);
      RETURN ODCICONST.SUCCESS;
    /*
    exception
    when others then
      if sqlcode = -29532 then
        jchem_core_pkg.handle_java_error(sqlerrm);
      else
        raise;
      end if;
    end;
    */
  END;

  STATIC FUNCTION ODCIIndexUpdate(ia sys.odciindexinfo, rid VARCHAR2, 
         oldval CLOB, newval CLOB, env sys.ODCIEnv) RETURN NUMBER IS 
    errmsg varchar2(32767);
  BEGIN
--    jchem_core_pkg.checkTableVersion(
--              ia.IndexSchema,
--              ia.IndexName,
--              ia.IndexCols(1).TableSchema,
--              ia.IndexCols(1).TableName,
--              ia.IndexCols(1).ColName,
--				ia.IndexCols(1).ColTypeName);

      errmsg := jchem_clob_pkg.update_mol_idxtable(oldval, newval, 
                      ia.IndexSchema, ia.IndexName, ia.IndexPartition,
                      ia.IndexCols(1).TableSchema, ia.IndexCols(1).TableName,
                      ia.IndexCols(1).ColName, ia.IndexCols(1).ColTypeName, rid);

      if errmsg is null then
        RETURN ODCICONST.SUCCESS;
      else
        jchem_core_pkg.handle_java_error(errmsg);
      end if;

      RETURN ODCICONST.SUCCESS;
  END;

  STATIC FUNCTION ODCIIndexInsert(ia sys.odciindexinfo, rid VARCHAR2, 
                          newval BLOB, env sys.ODCIEnv) RETURN NUMBER IS 
    errmsg varchar2(32767);
  BEGIN
--    jchem_core_pkg.checkTableVersion(
--              ia.IndexSchema,
--              ia.IndexName,
--              ia.IndexCols(1).TableSchema,
--              ia.IndexCols(1).TableName,
--              ia.IndexCols(1).ColName,
--				ia.IndexCols(1).ColTypeName);

      errmsg := jchem_blob_pkg.insert_mol_into_idxtable(newval, 
                      ia.IndexSchema, ia.IndexName, ia.IndexPartition,
                      ia.IndexCols(1).TableSchema, ia.IndexCols(1).TableName,
                      ia.IndexCols(1).ColName, ia.IndexCols(1).ColTypeName, 
                      rid);

      if errmsg is null then
        RETURN ODCICONST.SUCCESS;
      else
        jchem_core_pkg.handle_java_error(errmsg);
      end if;

      RETURN ODCICONST.SUCCESS;
  END;  

  STATIC FUNCTION ODCIIndexDelete(ia sys.odciindexinfo, rid VARCHAR2, 
			  oldval BLOB, env sys.ODCIEnv) RETURN NUMBER IS 
  BEGIN
    /*
    begin
    */
      jchem_core_pkg.delete_mol_from_idxtable(
          ia.IndexSchema, ia.IndexName, ia.IndexPartition,
                  ia.IndexCols(1).TableSchema, ia.IndexCols(1).TableName, rid);
      RETURN ODCICONST.SUCCESS;
    /*
    exception
    when others then
      if sqlcode = -29532 then
        jchem_core_pkg.handle_java_error(sqlerrm);
      else
        raise;
      end if;
    end;
    */
  END;
  
  STATIC FUNCTION ODCIIndexUpdate(ia sys.odciindexinfo, rid VARCHAR2, 
         oldval BLOB, newval BLOB, env sys.ODCIEnv) RETURN NUMBER IS 
    errmsg varchar2(32767);
  BEGIN
--    jchem_core_pkg.checkTableVersion(
--              ia.IndexSchema,
--              ia.IndexName,
--              ia.IndexCols(1).TableSchema,
--              ia.IndexCols(1).TableName,
--              ia.IndexCols(1).ColName,
--				ia.IndexCols(1).ColTypeName);

      errmsg := jchem_blob_pkg.update_mol_idxtable(oldval, newval,
                      ia.IndexSchema, ia.IndexName, ia.IndexPartition,
                      ia.IndexCols(1).TableSchema, ia.IndexCols(1).TableName,
                      ia.IndexCols(1).ColName, ia.IndexCols(1).ColTypeName, rid);    

      if errmsg is null then
        RETURN ODCICONST.SUCCESS;
      else
        jchem_core_pkg.handle_java_error(errmsg);
      end if;

      RETURN ODCICONST.SUCCESS;
  END;

  STATIC FUNCTION ODCIIndexExchangePartition(ia sys.ODCIIndexInfo, 
        ia1 sys.ODCIIndexInfo, env sys.ODCIEnv) RETURN NUMBER IS
  BEGIN
    /*
    begin
    */
      RETURN jchem_core_pkg.exchangePartitions(ia.IndexSchema,
          ia.IndexName, ia.IndexPartition, ia1.IndexSchema, ia1.IndexName);
    /*
    exception
    when others then
      if sqlcode = -29532 then
        jchem_core_pkg.handle_java_error(sqlerrm);
      else
        raise;
      end if;
    end;
    */
  END;

  STATIC FUNCTION ODCIIndexMergePartition(ia sys.ODCIIndexInfo, 
        part_name1 sys.ODCIPartInfo, part_name2 sys.ODCIPartInfo, 
        parms VARCHAR2, env sys.ODCIEnv) RETURN NUMBER IS
    ret number;
    errmsg varchar2(32767);
  BEGIN
      IF (ia.IndexPartition IS NOT NULL) THEN
        ret := jchem_core_pkg.indexDrop(ia.IndexSchema, ia.IndexName,
                  ia.IndexPartition, null, null);
      END IF;

      IF (part_name1.IndexPartition IS NOT NULL) THEN
        ret := jchem_core_pkg.indexDrop(ia.IndexSchema,
            ia.IndexName, part_name1.IndexPartition,
            null, null);
      END IF;

      IF (part_name2.IndexPartition IS NOT NULL) THEN
        errmsg := jchem_core_pkg.indexCreate(ia.IndexSchema, 
            ia.IndexName, part_name2.IndexPartition, ia.IndexCols(1).TableSchema, 
            ia.IndexCols(1).TableName, part_name2.TablePartition,
            ia.IndexCols(1).colName, ia.IndexCols(1).colTypeName, parms);

        if errmsg is null then
          ret := ODCICONST.SUCCESS;
        else
          jchem_core_pkg.handle_java_error(errmsg);
        end if;

      END IF;

      RETURN ret;
  END;

  STATIC FUNCTION ODCIIndexSplitPartition(ia sys.ODCIIndexInfo, 
        part_name1 sys.ODCIPartInfo, part_name2 sys.ODCIPartInfo, 
        parms VARCHAR2, env sys.ODCIEnv) RETURN NUMBER IS
    ret number;
    errmsg varchar2(32767);
  BEGIN
      IF (ia.IndexPartition IS NOT NULL) THEN
        ret := jchem_core_pkg.indexDrop(ia.IndexSchema, ia.IndexName,
                  ia.IndexPartition, null, null);
      END IF;

      IF (part_name1.IndexPartition IS NOT NULL) THEN
        errmsg := jchem_core_pkg.indexCreate(ia.IndexSchema, 
            ia.IndexName, part_name1.IndexPartition, ia.IndexCols(1).TableSchema, 
            ia.IndexCols(1).TableName, part_name1.TablePartition,
            ia.IndexCols(1).colName, ia.IndexCols(1).colTypeName, parms);

        if errmsg is null then
          ret := ODCICONST.SUCCESS;
        else
          jchem_core_pkg.handle_java_error(errmsg);
        end if;
      END IF;

      IF (part_name2.IndexPartition IS NOT NULL) THEN
        errmsg := jchem_core_pkg.indexCreate(ia.IndexSchema, 
            ia.IndexName, part_name2.IndexPartition, ia.IndexCols(1).TableSchema, 
            ia.IndexCols(1).TableName, part_name2.TablePartition,
            ia.IndexCols(1).colName, ia.IndexCols(1).colTypeName, parms);

        if errmsg is null then
          ret := ODCICONST.SUCCESS;
        else
          jchem_core_pkg.handle_java_error(errmsg);
        end if;
      END IF;

      RETURN ret;
  END;
END;
/
show errors;

-- Calls to this package incur a significant overhead -- always!!!
-- Use these procedures only where the call overhead is not relevant.
-- Turning log levels off will not reduce call overhead (unlike wth the logger
-- in Java code), it will reduce only the amount of logging.
create or replace package jcart_logger authid current_user as
    procedure set_log_level(logger_name varchar2, level number)
          parallel_enable as language java name
      'chemaxon.jchem.cartridge.oraresident.util.JCartPlsqlLogger.setLogLevel(java.lang.String, int)';

    procedure error(logger_name varchar2, message varchar2)
          parallel_enable as language java name
      'chemaxon.jchem.cartridge.oraresident.util.JCartPlsqlLogger.error(java.lang.String,
                                                            java.lang.String)';

    procedure warning(logger_name varchar2, message varchar2)
          parallel_enable as language java name
      'chemaxon.jchem.cartridge.oraresident.util.JCartPlsqlLogger.warning(java.lang.String,
                                                              java.lang.String)';
    procedure info(logger_name varchar2, message varchar2)
          parallel_enable as language java name
      'chemaxon.jchem.cartridge.oraresident.util.JCartPlsqlLogger.info(java.lang.String,
                                                           java.lang.String)';
    procedure debug(logger_name varchar2, message varchar2)
          parallel_enable as language java name
      'chemaxon.jchem.cartridge.oraresident.util.JCartPlsqlLogger.debug(java.lang.String,
                                                           java.lang.String)';

    procedure deprecated_method(method_name varchar2, replace_method_name varchar2)
          parallel_enable;
														   
end jcart_logger;
/
show errors;

create or replace package body jcart_logger as

    procedure deprecated_method(method_name varchar2, replace_method_name varchar2)
          parallel_enable is
	  warningstr VARCHAR2(4000);
	begin
	  warningstr:='Method '||method_name||' is deprecated';
	  if replace_method_name is not null then
	    warningstr:=warningstr||', use '||replace_method_name||' method instead.';
	  else
	    warningstr:=warningstr||'.';
	  end if;
	  warning('Deprecated method', warningstr);
	end;
	
end jcart_logger;
/
show errors;

------------------------------------------------------
-- CREATE OR REPLACE FUNCTIONAL IMPLEMENTATIONS for operators
------------------------------------------------------

CREATE OR REPLACE FUNCTION Exec_Func(sqlOperator VARCHAR2, target VARCHAR2,
                    query VARCHAR2, options VARCHAR2,
                    indexctx IN SYS.ODCIIndexCtx,
                    scanctx IN OUT jc_idxtype_im, scanflg IN NUMBER,
                    queryRequired IN NUMBER := 1)
      RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
  error varchar2(32767);
  result varchar2(32767);
  scanId number;
BEGIN
  jchem_core_pkg.checkTableVersionEx('JC_' || sqlOperator, indexctx);

  -- jchem_misc_pkg.trace('scanflg=' || to_char(scanflg));
  if scanctx is null then
    scanId := jchem_defright_pkg.get_next_scanId();
    scanctx := jc_idxtype_im(indexctx.IndexInfo, scanId, null);
  else
    scanId := scanctx.scan_num;
  end if;

  -- jchem_misc_pkg.trace('scanId=' || to_char(scanId));
  if scanflg = 1 then
    -- jchem_misc_pkg.trace('calling jchem_core_pkg.exec_function...');
    error := jchem_core_pkg.exec_function(
              null, null, null, null,
              null, null, null, null, null, null, null, null, scanId, scanflg,
              result);
    if error is null then
      -- jchem_misc_pkg.trace('error is null');
      return null;
    else
      -- jchem_misc_pkg.trace('error is not null');
      jchem_core_pkg.handle_java_error(error);
    end if;
  end if;

  IF query is null and queryRequired = 1 THEN 
    return null;
  END IF;

  IF target is null THEN
    IF indexctx.IndexInfo IS NOT NULL THEN
      error := jchem_core_pkg.exec_function(
                    sqlOperator, null, query, options,
                    indexctx.Rid,
                    indexctx.IndexInfo.IndexSchema,
                    indexctx.IndexInfo.IndexName,
                    indexctx.IndexInfo.IndexPartition,
                    indexctx.IndexInfo.IndexCols(1).TableSchema,
                    indexctx.IndexInfo.IndexCols(1).TableName,
                    indexctx.IndexInfo.IndexCols(1).ColName,
					indexctx.IndexInfo.IndexCols(1).ColTypeName,
                    scanId,
                    scanFlg,
                    result);
    ELSE 
      result := null;
    END IF;
  ELSE
    IF indexctx.IndexInfo IS NOT NULL THEN
      error := jchem_core_pkg.exec_function(
                    sqlOperator, target, query, options,
                    indexctx.Rid,
                    indexctx.IndexInfo.IndexSchema,
                    indexctx.IndexInfo.IndexName,
                    indexctx.IndexInfo.IndexPartition,
                    indexctx.IndexInfo.IndexCols(1).TableSchema,
                    indexctx.IndexInfo.IndexCols(1).TableName,
                    indexctx.IndexInfo.IndexCols(1).ColName,
					indexctx.IndexInfo.IndexCols(1).ColTypeName,
                    scanId,
                    scanFlg,
                    result);
    ELSE
      error := jchem_core_pkg.exec_function(
                sqlOperator, target, query, options,
                null, null, null, null, null, null, null, null, null, null,
                result);
    END IF;
  END IF;

--  raise_application_error(-20102,
--    'error=' || error || ', result=' || result);

  if error is null then
    return result;
  else
    jchem_core_pkg.handle_java_error(error);
  end if;
END;
/
show errors;

CREATE OR REPLACE FUNCTION Exec_FuncV(sqlOperator VARCHAR2,
                              target CLOB, query CLOB, options VARCHAR2,
                              indexctx IN SYS.ODCIIndexCtx, 
                              scanctx IN OUT jc_idxtype_im,
                              scanflg IN NUMBER,
                              queryRequired IN NUMBER := 1)
      RETURN VARCHAR2 AUTHID CURRENT_USER PARALLEL_ENABLE AS
  error varchar2(32767);
  result varchar2(32767);
  scanId number;
BEGIN
  jchem_core_pkg.checkTableVersionEx('JC_' || sqlOperator, indexctx);

  if indexctx.IndexInfo is not null then
    if scanctx is null then
      scanId := jchem_defright_pkg.get_next_scanId();
      scanctx := jc_idxtype_im(indexctx.IndexInfo, scanId, null);
    else
      scanId := scanctx.scan_num;
    end if;
  end if;

  IF query is null and queryRequired = 1 THEN 
    return null; 
  END IF;
  IF target is null THEN
        return NULL;
  ELSE
    IF indexctx.IndexInfo IS NOT NULL THEN
      error := jchem_clob_pkg.exec_function(
                    sqlOperator, target, query, options,
                    indexctx.Rid,
                    indexctx.IndexInfo.IndexSchema,
                    indexctx.IndexInfo.IndexName,
                    indexctx.IndexInfo.IndexPartition,
                    indexctx.IndexInfo.IndexCols(1).TableSchema,
                    indexctx.IndexInfo.IndexCols(1).TableName,
                    indexctx.IndexInfo.IndexCols(1).ColName, 
					indexctx.IndexInfo.IndexCols(1).ColTypeName,
					scanId,
                    scanflg, result);

    ELSE
      error := jchem_clob_pkg.exec_function(
                sqlOperator, target, query, options,
                null, null, null, null, null, null, null, null, null, null, result);
    END IF;
  END IF;

  if error is null then
    return result;
  else
    jchem_core_pkg.handle_java_error(error);
  end if;
END;
/
show errors;

CREATE OR REPLACE FUNCTION Exec_FuncC(sqlOperator VARCHAR2,
                              target CLOB, query CLOB, options VARCHAR2,
                              indexctx IN SYS.ODCIIndexCtx, 
                              scanctx IN OUT jc_idxtype_im,
                              scanflg IN NUMBER,
                              queryRequired IN NUMBER := 1)
      RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
  error varchar2(32767);
  result varchar2(32767);
  scanId number;
BEGIN
  jchem_core_pkg.checkTableVersionEx('JC_' || sqlOperator, indexctx);

  if indexctx.IndexInfo is not null then
    if scanctx is null then
      scanId := jchem_defright_pkg.get_next_scanId();
      scanctx := jc_idxtype_im(indexctx.IndexInfo, scanId, null);
    else
      scanId := scanctx.scan_num;
    end if;
  end if;

  IF query is null and queryRequired = 1 THEN 
    return null; 
  END IF;
  IF target is null THEN
        return NULL;
  ELSE
    IF indexctx.IndexInfo IS NOT NULL THEN
      error := jchem_clob_pkg.exec_function(
                    sqlOperator, target, query, options,
                    indexctx.Rid,
                    indexctx.IndexInfo.IndexSchema,
                    indexctx.IndexInfo.IndexName,
                    indexctx.IndexInfo.IndexPartition,
                    indexctx.IndexInfo.IndexCols(1).TableSchema,
                    indexctx.IndexInfo.IndexCols(1).TableName,
                    indexctx.IndexInfo.IndexCols(1).ColName,
					indexctx.IndexInfo.IndexCols(1).ColTypeName,
                    scanId, scanflg, result);

    ELSE
      error := jchem_clob_pkg.exec_function(
                sqlOperator, target, query, options,
                null, null, null, null, null, null, null, null,
                null, null, result);
    END IF;
  END IF;

  if error is null then
    return result;
  else
    jchem_core_pkg.handle_java_error(error);
  end if;
END;
/
show errors;

CREATE OR REPLACE FUNCTION Exec_FuncB(sqlOperator VARCHAR2,
                              target BLOB, query BLOB, options VARCHAR2,
                              indexctx IN SYS.ODCIIndexCtx, 
                              scanctx IN OUT jc_idxtype_im,
                              scanflg IN NUMBER,
                              queryRequired IN NUMBER := 1)
      RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
  error varchar2(32767);
  result varchar2(32767);
  scanId number;
BEGIN
  jchem_core_pkg.checkTableVersionEx('JC_' || sqlOperator, indexctx);

  if indexctx.IndexInfo is not null then
    if scanctx is null then
      scanId := jchem_defright_pkg.get_next_scanId();
      scanctx := jc_idxtype_im(indexctx.IndexInfo, scanId, null);
    else
      scanId := scanctx.scan_num;
    end if;
  end if;

  IF query is null and queryRequired = 1 THEN 
    return null; 
  END IF;
  IF target is null THEN
        return NULL;
  ELSE
    IF indexctx.IndexInfo IS NOT NULL THEN
      error := jchem_blob_pkg.exec_function(
                    sqlOperator, target, query, options,
                    indexctx.Rid,
                    indexctx.IndexInfo.IndexSchema,
                    indexctx.IndexInfo.IndexName,
                    indexctx.IndexInfo.IndexPartition,
                    indexctx.IndexInfo.IndexCols(1).TableSchema,
                    indexctx.IndexInfo.IndexCols(1).TableName,
                    indexctx.IndexInfo.IndexCols(1).ColName,
					indexctx.IndexInfo.IndexCols(1).ColTypeName,
                    scanId, scanflg, result);

    ELSE
      error := jchem_blob_pkg.exec_function(
                sqlOperator, target, query, options,
                null, null, null, null, null, null, null, null,
                null, null, result);
    END IF;
  END IF;

  if error is null then
    return result;
  else
    jchem_core_pkg.handle_java_error(error);
  end if;
END;
/
show errors;

CREATE OR REPLACE FUNCTION Exec_FuncCC(sqlOperator VARCHAR2,
                              target CLOB, query CLOB, options VARCHAR2,
                              indexctx IN SYS.ODCIIndexCtx, 
                              scanctx IN OUT jc_idxtype_im,
                              scanflg IN NUMBER,
                              queryRequired IN NUMBER := 1)
      RETURN CLOB AUTHID CURRENT_USER PARALLEL_ENABLE AS
  error varchar2(32767);
  result CLOB;
  scanId number;
BEGIN
  jchem_core_pkg.checkTableVersionEx('JC_' || sqlOperator, indexctx);

  if indexctx.IndexInfo is not null then
    if scanctx is null then
      scanId := jchem_defright_pkg.get_next_scanId();
      scanctx := jc_idxtype_im(indexctx.IndexInfo, scanId, null);
    else
      scanId := scanctx.scan_num;
    end if;
  end if;

  IF query is null and queryRequired = 1 THEN 
    return null; 
  END IF;
  IF target is null THEN
        return NULL;
  ELSE
    IF indexctx.IndexInfo IS NOT NULL THEN
      error := jchem_clob_pkg.exec_function__c(
                    sqlOperator, target, query, options,
                    indexctx.Rid,
                    indexctx.IndexInfo.IndexSchema,
                    indexctx.IndexInfo.IndexName,
                    indexctx.IndexInfo.IndexPartition,
                    indexctx.IndexInfo.IndexCols(1).TableSchema,
                    indexctx.IndexInfo.IndexCols(1).TableName,
                    indexctx.IndexInfo.IndexCols(1).ColName,
					indexctx.IndexInfo.IndexCols(1).ColTypeName,
                    scanId, scanflg, result);

    ELSE
      error := jchem_clob_pkg.exec_function__c(
                sqlOperator, target, query, options,
                null, null, null, null, null, null, null, null,
                null, null, result);
    END IF;
  END IF;

  if error is null then
    return result;
  else
    jchem_core_pkg.handle_java_error(error);
  end if;
END;
/
show errors;

CREATE OR REPLACE FUNCTION Exec_FuncBB(sqlOperator VARCHAR2,
                              target BLOB, query BLOB, options VARCHAR2,
                              indexctx IN SYS.ODCIIndexCtx, 
                              scanctx IN OUT jc_idxtype_im,
                              scanflg IN NUMBER,
                              queryRequired IN NUMBER := 1)
      RETURN BLOB AUTHID CURRENT_USER PARALLEL_ENABLE AS
  error varchar2(32767);
  result BLOB;
  scanId number;
BEGIN
  jchem_core_pkg.checkTableVersionEx('JC_' || sqlOperator, indexctx);

  if indexctx.IndexInfo is not null then
    if scanctx is null then
      scanId := jchem_defright_pkg.get_next_scanId();
      scanctx := jc_idxtype_im(indexctx.IndexInfo, scanId, null);
    else
      scanId := scanctx.scan_num;
    end if;
  end if;

  IF query is null and queryRequired = 1 THEN 
    return null; 
  END IF;
  IF target is null THEN
        return NULL;
  ELSE
    IF indexctx.IndexInfo IS NOT NULL THEN
      error := jchem_blob_pkg.exec_function__b(
                    sqlOperator, target, query, options,
                    indexctx.Rid,
                    indexctx.IndexInfo.IndexSchema,
                    indexctx.IndexInfo.IndexName,
                    indexctx.IndexInfo.IndexPartition,
                    indexctx.IndexInfo.IndexCols(1).TableSchema,
                    indexctx.IndexInfo.IndexCols(1).TableName,
                    indexctx.IndexInfo.IndexCols(1).ColName,
					indexctx.IndexInfo.IndexCols(1).ColTypeName,
                    scanId, scanflg, result);

    ELSE
      error := jchem_blob_pkg.exec_function__b(
                sqlOperator, target, query, options,
                null, null, null, null, null, null, null, null,
                null, null, result);
    END IF;
  END IF;

  if error is null then
    return result;
  else
    jchem_core_pkg.handle_java_error(error);
  end if;
END;
/
show errors;

CREATE OR REPLACE FUNCTION Exec_FuncBV(sqlOperator VARCHAR2,
                              target BLOB, param VARCHAR2, options VARCHAR2,
                              indexctx IN SYS.ODCIIndexCtx, 
                              scanctx IN OUT jc_idxtype_im,
                              scanflg IN NUMBER)
      RETURN VARCHAR2 AUTHID CURRENT_USER PARALLEL_ENABLE AS
  error varchar2(32767);
  result VARCHAR2(32767);
  scanId number;
BEGIN
  jchem_core_pkg.checkTableVersionEx('JC_' || sqlOperator, indexctx);

  if indexctx.IndexInfo is not null then
    if scanctx is null then
      scanId := jchem_defright_pkg.get_next_scanId();
      scanctx := jc_idxtype_im(indexctx.IndexInfo, scanId, null);
    else
      scanId := scanctx.scan_num;
    end if;
  end if;

  IF target is null THEN
        return NULL;
  ELSE
    IF indexctx.IndexInfo IS NOT NULL THEN
      error := jchem_blob_pkg.exec_function_bv(
                    sqlOperator, target, param, options,
                    indexctx.Rid,
                    indexctx.IndexInfo.IndexSchema,
                    indexctx.IndexInfo.IndexName,
                    indexctx.IndexInfo.IndexPartition,
                    indexctx.IndexInfo.IndexCols(1).TableSchema,
                    indexctx.IndexInfo.IndexCols(1).TableName,
                    indexctx.IndexInfo.IndexCols(1).ColName,
					indexctx.IndexInfo.IndexCols(1).ColTypeName,
                    scanId, scanflg, result);                
        
    ELSE
      error := jchem_blob_pkg.exec_function_bv(
                sqlOperator, target, param, options,
                null, null, null, null, null, null, null, null,
                null, null, result);
    END IF;
  END IF;

  if error is null then
    return result;
  else
    jchem_core_pkg.handle_java_error(error);
  end if;
END;
/
show errors;

CREATE OR REPLACE FUNCTION Contains_Func_deprecated(target VARCHAR2, query VARCHAR2,
                              indexctx IN SYS.ODCIIndexCtx, 
                              scanctx IN OUT jc_idxtype_im,
                              scanflg IN NUMBER)
RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
BEGIN
  jcart_logger.deprecated_method('JC_CONTAINS','JC_COMPARE');
  return Exec_Func('CONTAINS', target, query, null, indexctx, scanctx, scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION Contains_FuncC_deprecated(target CLOB, query CLOB,
                              indexctx IN SYS.ODCIIndexCtx, 
                              scanctx IN OUT jc_idxtype_im,
                              scanflg IN NUMBER)
RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
BEGIN
  jcart_logger.deprecated_method('JC_CONTAINS','JC_COMPARE');
  return Exec_FuncC('CONTAINS', target, query, null, indexctx, scanctx, scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION Contains_FuncB_deprecated(target BLOB, query BLOB,
                              indexctx IN SYS.ODCIIndexCtx, 
                              scanctx IN OUT jc_idxtype_im,
                              scanflg IN NUMBER)
RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
BEGIN
  jcart_logger.deprecated_method('JC_CONTAINS','JC_COMPARE');
  return Exec_FuncB('CONTAINS', target, query, null, indexctx, scanctx, scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION ContainsB_FuncB_deprecated(target BLOB, query BLOB,
                              indexctx IN SYS.ODCIIndexCtx, 
                              scanctx IN OUT jc_idxtype_im,
                              scanflg IN NUMBER)
RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
BEGIN
  jcart_logger.deprecated_method('JC_CONTAINSB','JC_COMPARE');
  return Contains_FuncB_deprecated(target,query,indexctx,scanctx,scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION Equals_Func_deprecated(target VARCHAR2, query VARCHAR2,
                            indexctx IN SYS.ODCIIndexCtx, 
                            scanctx IN OUT jc_idxtype_im,
                            scanflg IN NUMBER)
    RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
BEGIN
  jcart_logger.deprecated_method('JC_EQUALS','JC_COMPARE');
  return Exec_Func('EQUALS', target, query, null, indexctx, scanctx, scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION Equals_FuncC_deprecated(target CLOB, query CLOB,
                            indexctx IN SYS.ODCIIndexCtx, 
                            scanctx IN OUT jc_idxtype_im,
                            scanflg IN NUMBER)
    RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
BEGIN
  jcart_logger.deprecated_method('JC_EQUALS','JC_COMPARE');
  return Exec_FuncC('EQUALS', target, query, null, indexctx, scanctx, scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION Equals_FuncB_deprecated(target BLOB, query BLOB,
                            indexctx IN SYS.ODCIIndexCtx, 
                            scanctx IN OUT jc_idxtype_im,
                            scanflg IN NUMBER)
    RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
BEGIN
  jcart_logger.deprecated_method('JC_EQUALS','JC_COMPARE');
  return Exec_FuncB('EQUALS', target, query, null, indexctx, scanctx, scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION EqualsB_FuncB_deprecated(target BLOB, query BLOB,
                            indexctx IN SYS.ODCIIndexCtx, 
                            scanctx IN OUT jc_idxtype_im,
                            scanflg IN NUMBER)
RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
BEGIN
  jcart_logger.deprecated_method('JC_EQUALSB','JC_COMPARE');
  return Equals_FuncB_deprecated(target,query,indexctx,scanctx,scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION Evaluate_Func(target VARCHAR2, ct VARCHAR2,
                                         indexctx IN SYS.ODCIIndexCtx, 
                                         scanctx IN OUT jc_idxtype_im,
                                         scanflg IN NUMBER)
      RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
BEGIN
  return Exec_Func('EVALUATE', target, null, ct,
                   indexctx, scanctx, scanflg, 0);
END;
/
show errors;

CREATE OR REPLACE FUNCTION Evaluate_FuncC(target CLOB, ct VARCHAR2,
                                         indexctx IN SYS.ODCIIndexCtx, 
                                         scanctx IN OUT jc_idxtype_im,
                                         scanflg IN NUMBER)
      RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
BEGIN
  return Exec_FuncC('EVALUATE', target, null, ct,
                    indexctx, scanctx, scanflg, 0);
END;
/
show errors;

CREATE OR REPLACE FUNCTION Evaluate_FuncB(target BLOB, ct VARCHAR2,
                                         indexctx IN SYS.ODCIIndexCtx, 
                                         scanctx IN OUT jc_idxtype_im,
                                         scanflg IN NUMBER)
      RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
BEGIN
  return Exec_FuncB('EVALUATE', target, null, ct,
                    indexctx, scanctx, scanflg, 0);
END;
/
show errors;

CREATE OR REPLACE FUNCTION Evaluate3_Func(target VARCHAR2, ct VARCHAR2,
                                          options VARCHAR2,
                                          indexctx IN SYS.ODCIIndexCtx, 
                                          scanctx IN OUT jc_idxtype_im,
                                          scanflg IN NUMBER)
      RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
BEGIN
  return Exec_Func('EVALUATE', target, options, ct,
                   indexctx, scanctx, scanflg, 0);
END;
/
show errors;

CREATE OR REPLACE FUNCTION Evaluate3_FuncC(target CLOB, ct VARCHAR2,
                                          options VARCHAR2,
                                          indexctx IN SYS.ODCIIndexCtx, 
                                          scanctx IN OUT jc_idxtype_im,
                                          scanflg IN NUMBER)
      RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
BEGIN
  return Exec_FuncC('EVALUATE', target, options, ct,
                   indexctx, scanctx, scanflg, 0);
END;
/
show errors;

CREATE OR REPLACE FUNCTION Evaluate3_FuncB(target BLOB, ct VARCHAR2,
                                          options VARCHAR2,
                                          indexctx IN SYS.ODCIIndexCtx, 
                                          scanctx IN OUT jc_idxtype_im,
                                          scanflg IN NUMBER)
      RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
BEGIN
  return Exec_FuncBV('EVALUATE', target, options, ct,
                   indexctx, scanctx, scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION EvaluateX_Func(target VARCHAR2, options VARCHAR2,
                                         indexctx IN SYS.ODCIIndexCtx, 
                                         scanctx IN OUT jc_idxtype_im,
                                         scanflg IN NUMBER)
      RETURN VARCHAR2 AUTHID CURRENT_USER PARALLEL_ENABLE AS
BEGIN
  return Exec_FuncV('EVALUATE_X', target, null, options, indexctx, scanctx, scanflg, 0);
END;
/
show errors;

CREATE OR REPLACE FUNCTION EvaluateX_FuncC(target CLOB, options VARCHAR2,
                           indexctx IN SYS.ODCIIndexCtx, 
                           scanctx IN OUT jc_idxtype_im,
                           scanflg IN NUMBER)
      RETURN CLOB AUTHID CURRENT_USER PARALLEL_ENABLE AS
BEGIN
  return Exec_FuncCC('EVALUATE_X', target, null, options, indexctx, scanctx, scanflg, 0);
END;
/
show errors;

CREATE OR REPLACE FUNCTION EvaluateX_FuncB(target BLOB, options VARCHAR2,
                           indexctx IN SYS.ODCIIndexCtx, 
                           scanctx IN OUT jc_idxtype_im,
                           scanflg IN NUMBER)
      RETURN BLOB AUTHID CURRENT_USER PARALLEL_ENABLE AS
BEGIN
  return Exec_FuncBB('EVALUATE_X', target, null, options, indexctx, scanctx, scanflg, 0);
END;
/
show errors;

CREATE OR REPLACE FUNCTION Compare_Func(target VARCHAR2, query VARCHAR2,
                             options VARCHAR2, indexctx IN SYS.ODCIIndexCtx, 
                             scanctx IN OUT jc_idxtype_im,
                             scanflg IN NUMBER)
      RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
BEGIN
  return Exec_Func('COMPARE', target, query, options, indexctx, scanctx, scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION Compare_FuncC(target CLOB, query CLOB,
                             options VARCHAR2, indexctx IN SYS.ODCIIndexCtx, 
                             scanctx IN OUT jc_idxtype_im,
                             scanflg IN NUMBER)
      RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
BEGIN
  return Exec_FuncC('COMPARE', target, query, options, indexctx, scanctx, scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION Compare_FuncCV(target CLOB, query VARCHAR2,
                             options VARCHAR2, indexctx IN SYS.ODCIIndexCtx, 
                             scanctx IN OUT jc_idxtype_im,
                             scanflg IN NUMBER)
      RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
  error varchar2(32767);
  result varchar2(32767);
  scanId NUMBER;
  BEGIN
    jchem_core_pkg.checkTableVersionEx('JC_COMPARE', indexctx);

  if indexctx.IndexInfo is not null then
    if scanctx is null then
      scanId := jchem_defright_pkg.get_next_scanId();
      scanctx := jc_idxtype_im(indexctx.IndexInfo, scanId, null);
    else
      scanId := scanctx.scan_num;
    end if;
  end if;

    IF query is null THEN 
      return null; 
    END IF;
    -- TODO: Should not we check 'indexctx IS NULL'?
    -- TODO: Should not we check 'scanflg = SYS.CleanUpCall'?
    IF target is null THEN
	  return NULL;
    ELSE
      IF indexctx.IndexInfo IS NOT NULL THEN
        error := jchem_clob_pkg.exec_function_cv('COMPARE',
				target, query, options,
                                indexctx.Rid,
				indexctx.IndexInfo.IndexSchema,
				indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
				indexctx.IndexInfo.IndexCols(1).TableSchema, 
				indexctx.IndexInfo.IndexCols(1).TableName,
                indexctx.IndexInfo.IndexCols(1).ColName,
				indexctx.IndexInfo.IndexCols(1).ColTypeName,
                scanId, scanflg, result);
      ELSE
        error := jchem_clob_pkg.exec_function_cv('COMPARE',
                    target, query, options,
                    null, null, null, null, null, null, null, null,
                    null, null, result);
      END IF;
    END IF;

  if error is null then
    return result;
  else
    jchem_core_pkg.handle_java_error(error);
  end if;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION Compare_FuncCB(target CLOB, query BLOB,
                             options VARCHAR2, indexctx IN SYS.ODCIIndexCtx, 
                             scanctx IN OUT jc_idxtype_im,
                             scanflg IN NUMBER)
RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
  error varchar2(32767);
  result varchar2(32767);
  scanId number;
  BEGIN
    jchem_core_pkg.checkTableVersionEx('JC_COMPARE', indexctx);

  if indexctx.IndexInfo is not null then
    if scanctx is null then
      scanId := jchem_defright_pkg.get_next_scanId();
      scanctx := jc_idxtype_im(indexctx.IndexInfo, scanId, null);
    else
      scanId := scanctx.scan_num;
    end if;
  end if;

    IF query is null THEN 
      return null; 
    END IF;
    -- TODO: Should not we check 'indexctx IS NULL'?
    -- TODO: Should not we check 'scanflg = SYS.CleanUpCall'?
    IF target is null THEN
	  return NULL;
    ELSE
      IF indexctx.IndexInfo IS NOT NULL THEN
        error := jchem_clob_pkg.exec_function_cb('COMPARE',
				target, query, options,
                                indexctx.Rid,
				indexctx.IndexInfo.IndexSchema,
				indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
				indexctx.IndexInfo.IndexCols(1).TableSchema, 
				indexctx.IndexInfo.IndexCols(1).TableName,
                indexctx.IndexInfo.IndexCols(1).ColName,
				indexctx.IndexInfo.IndexCols(1).ColTypeName,
                scanId, scanflg, result);
      ELSE
        error := jchem_clob_pkg.exec_function_cb('COMPARE',
                    target, query, options,
                    null, null, null, null, null, null, null, null,
                    null, null, result);
      END IF;
    END IF;

  if error is null then
    return result;
  else
    jchem_core_pkg.handle_java_error(error);
  end if;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION Compare_FuncB(target BLOB, query BLOB,
                             options VARCHAR2, indexctx IN SYS.ODCIIndexCtx, 
                             scanctx IN OUT jc_idxtype_im,
                             scanflg IN NUMBER)
      RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
BEGIN
  return Exec_FuncB('COMPARE', target, query, options, indexctx, scanctx, scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION Compare_FuncBV(target BLOB, query VARCHAR2,
                             options VARCHAR2, indexctx IN SYS.ODCIIndexCtx, 
                             scanctx IN OUT jc_idxtype_im,
                             scanflg IN NUMBER)
RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
  error varchar2(32767);
  result varchar2(32767);
  scanId number;
  BEGIN
    jchem_core_pkg.checkTableVersionEx('JC_COMPARE', indexctx);

  if indexctx.IndexInfo is not null then
    if scanctx is null then
      scanId := jchem_defright_pkg.get_next_scanId();
      scanctx := jc_idxtype_im(indexctx.IndexInfo, scanId, null);
    else
      scanId := scanctx.scan_num;
    end if;
  end if;

    IF query is null THEN 
      return null; 
    END IF;
    -- TODO: Should not we check 'indexctx IS NULL'?
    -- TODO: Should not we check 'scanflg = SYS.CleanUpCall'?
    IF target is null THEN
	  return NULL;
    ELSE
      IF indexctx.IndexInfo IS NOT NULL THEN
       error := jchem_blob_pkg.exec_function_bv('COMPARE',
				target, query, options,
                                indexctx.Rid,
				indexctx.IndexInfo.IndexSchema,
				indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
				indexctx.IndexInfo.IndexCols(1).TableSchema, 
				indexctx.IndexInfo.IndexCols(1).TableName,
                indexctx.IndexInfo.IndexCols(1).ColName,
				indexctx.IndexInfo.IndexCols(1).ColTypeName,
                scanId, scanflg, result);
      ELSE
        error := jchem_blob_pkg.exec_function_bv('COMPARE', target, query, options,
                                           null, null, null, null, null, null, null, null,
                                           null, null, result);
      END IF;
    END IF;

  if error is null then
    return result;
  else
    jchem_core_pkg.handle_java_error(error);
  end if;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION compare_func_bquery(target VARCHAR2, query BLOB,
                             options VARCHAR2, indexctx IN SYS.ODCIIndexCtx, 
                             scanctx IN OUT jc_idxtype_im,
                             scanflg IN NUMBER)
      RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
  error varchar2(32767);
  result varchar2(32767);
  scanId number;
BEGIN
    jchem_core_pkg.checkTableVersionEx('JC_COMPARE', indexctx);

  if indexctx.IndexInfo is not null then
    if scanctx is null then
      scanId := jchem_defright_pkg.get_next_scanId();
      scanctx := jc_idxtype_im(indexctx.IndexInfo, scanId, null);
    else
      scanId := scanctx.scan_num;
    end if;
  end if;

    IF query is null THEN 
      return null; 
    END IF;
    -- TODO: Should not we check 'indexctx IS NULL'?
    -- TODO: Should not we check 'scanflg = SYS.CleanUpCall'?
    IF target is null THEN
	  return NULL;
    ELSE
      IF indexctx.IndexInfo IS NOT NULL THEN
       error := jchem_blob_pkg.exec_function_vb('COMPARE',
				target, query, options,
                                indexctx.Rid,
				indexctx.IndexInfo.IndexSchema,
				indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
				indexctx.IndexInfo.IndexCols(1).TableSchema, 
				indexctx.IndexInfo.IndexCols(1).TableName,
                indexctx.IndexInfo.IndexCols(1).ColName,
				indexctx.IndexInfo.IndexCols(1).ColTypeName,
                scanId, scanflg, result);
      ELSE
        error := jchem_blob_pkg.exec_function_vb('COMPARE', target, query, options,
                                               null, null, null, null, null, null, null, null,
                                               null, null, result);
      END IF;
    END IF;
  if error is null then
    return result;
  else
    jchem_core_pkg.handle_java_error(error);
  end if;
END;
/
show errors;

CREATE OR REPLACE FUNCTION MatchCount_Func(target VARCHAR2, query VARCHAR2,
                                indexctx IN SYS.ODCIIndexCtx,
                                scanctx IN OUT jc_idxtype_im,
                                scanflg IN NUMBER)
      RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
BEGIN
  return Exec_Func('MATCHCOUNT', target, query, null, indexctx, scanctx, scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION MatchCount_FuncC(target CLOB, query CLOB,
				indexctx IN SYS.ODCIIndexCtx, 
                                scanctx IN OUT jc_idxtype_im,
                                scanflg IN NUMBER)
      RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
BEGIN
  return Exec_FuncC('MATCHCOUNT', target, query, null, indexctx, scanctx, scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION MatchCount_FuncB(target BLOB, query BLOB,
                                      indexctx IN SYS.ODCIIndexCtx, 
                                      scanctx IN OUT jc_idxtype_im,
                                      scanflg IN NUMBER)
      RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
BEGIN
  return Exec_FuncB('MATCHCOUNT', target, query, null, indexctx, scanctx, scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION MatchCountB_FuncB_deprecated(target BLOB, query BLOB,
                                      indexctx IN SYS.ODCIIndexCtx, 
                                      scanctx IN OUT jc_idxtype_im,
                                      scanflg IN NUMBER)
      RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
BEGIN
  jcart_logger.deprecated_method('JC_MATCHCOUNTB','JC_MATCHCOUNT');
  return MatchCount_FuncB(target,query,indexctx,scanctx,scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION Tanimoto_Func(target VARCHAR2, query VARCHAR2, 
          indexctx IN SYS.ODCIIndexCtx, scanctx IN OUT jc_idxtype_im,
          scanflg IN NUMBER) RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
  res VARCHAR2(32767);
BEGIN
  return Exec_Func('TANIMOTO', target, query, 'dissimilarityMetric:tanimoto',
                   indexctx, scanctx, scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION Tanimoto_FuncC(target CLOB, query CLOB, 
          indexctx IN SYS.ODCIIndexCtx, scanctx IN OUT jc_idxtype_im,
          scanflg IN NUMBER)
RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
  res VARCHAR2(32767);
BEGIN
  return Exec_FuncC('TANIMOTO', target, query, 'dissimilarityMetric:tanimoto',
                   indexctx, scanctx, scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION Tanimoto_FuncB(target BLOB, query BLOB, 
          indexctx IN SYS.ODCIIndexCtx, scanctx IN OUT jc_idxtype_im,
          scanflg IN NUMBER)
RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
  res VARCHAR2(32767);
BEGIN
  return Exec_FuncB('TANIMOTO', target, query, 'dissimilarityMetric:tanimoto',
                   indexctx, scanctx, scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION TanimotoB_FuncB_deprecated(target BLOB, query BLOB, 
          indexctx IN SYS.ODCIIndexCtx, scanctx IN OUT jc_idxtype_im,
          scanflg IN NUMBER)
RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
  res VARCHAR2(32767);
BEGIN
  jcart_logger.deprecated_method('JC_TANIMOTOB','JC_TANIMOTO');
  return Tanimoto_FuncB(target,query,indexctx,scanctx,scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION Tversky_Func(target in VARCHAR2, query in VARCHAR2, 
          targetWeight in number, queryWeight in number,
          indexctx IN SYS.ODCIIndexCtx, scanctx IN OUT jc_idxtype_im,
          scanflg IN NUMBER) RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
    res VARCHAR2(32767);
BEGIN
  return Exec_Func('TVERSKY', target, query, 'dissimilarityMetric:tversky'
                    || jchem_core_pkg.simCalcSeparator 
                    || queryWeight || jchem_core_pkg.simCalcSeparator || targetWeight,
                   indexctx, scanctx, scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION Tversky_FuncC(target CLOB, query CLOB, 
          targetWeight in number, queryWeight in number,
          indexctx IN SYS.ODCIIndexCtx, 
          scanctx IN OUT jc_idxtype_im,
          scanflg IN NUMBER)
RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
    res VARCHAR2(32767);
BEGIN
  return Exec_FuncC('TVERSKY', target, query, 'dissimilarityMetric:tversky'
                    || jchem_core_pkg.simCalcSeparator 
                    || queryWeight || jchem_core_pkg.simCalcSeparator || targetWeight,
                   indexctx, scanctx, scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION Tversky_FuncB(target BLOB, query BLOB, 
                  targetWeight in number, queryWeight in number,
                  indexctx IN SYS.ODCIIndexCtx, 
                  scanctx IN OUT jc_idxtype_im,
                  scanflg IN NUMBER)
RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
    res VARCHAR2(32767);
BEGIN
  return Exec_FuncB('TVERSKY', target, query, 'dissimilarityMetric:tversky'
                    || jchem_core_pkg.simCalcSeparator 
                    || queryWeight || jchem_core_pkg.simCalcSeparator || targetWeight,
                   indexctx, scanctx, scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION Dissimilarity_Func(target VARCHAR2, query VARCHAR2, 
          indexctx IN SYS.ODCIIndexCtx, scanctx IN OUT jc_idxtype_im,
          scanflg IN NUMBER)
RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS 
    res VARCHAR2(32767);
BEGIN
  return Exec_Func('DISSIMILARITY', target, query, 'dissimilarityMetric:tanimoto',
                   indexctx, scanctx, scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION Dissimilarity_FuncC(target CLOB, query CLOB, 
          indexctx IN SYS.ODCIIndexCtx, scanctx IN OUT jc_idxtype_im,
          scanflg IN NUMBER)
RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS 
  res VARCHAR2(32767);
BEGIN
  return Exec_FuncC('DISSIMILARITY', target, query, 'dissimilarityMetric:tanimoto',
                   indexctx, scanctx, scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION Dissimilarity_FuncB(target BLOB, query BLOB, 
          indexctx IN SYS.ODCIIndexCtx, scanctx IN OUT jc_idxtype_im,
          scanflg IN NUMBER)
RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS 
    res VARCHAR2(32767);
BEGIN
  return Exec_FuncB('DISSIMILARITY', target, query, 'dissimilarityMetric:tanimoto',
                   indexctx, scanctx, scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION DissimB_FuncB_deprecated(target BLOB, query BLOB, 
          indexctx IN SYS.ODCIIndexCtx, scanctx IN OUT jc_idxtype_im,
          scanflg IN NUMBER)
RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS 
    res VARCHAR2(32767);
BEGIN
  jcart_logger.deprecated_method('JC_DISSIMILARITYB','JC_DISSIMILARITY');
  return Dissimilarity_FuncB(target,query,indexctx,scanctx,scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION DissimilarityX_Func(target VARCHAR2, query VARCHAR2, 
          options VARCHAR2,
          indexctx IN SYS.ODCIIndexCtx, scanctx IN OUT jc_idxtype_im,
          scanflg IN NUMBER)
RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS 
    res VARCHAR2(32767);
BEGIN
  return Exec_Func('DISSIMILARITY', target, query, options,
                   indexctx, scanctx, scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION DissimilarityX_FuncC(target CLOB, query CLOB, 
          options VARCHAR2,
          indexctx IN SYS.ODCIIndexCtx, scanctx IN OUT jc_idxtype_im,
          scanflg IN NUMBER)
RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS 
  res VARCHAR2(32767);
BEGIN
  return Exec_FuncC('DISSIMILARITY', target, query, options,
                   indexctx, scanctx, scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION DissimilarityX_FuncB(target BLOB, query BLOB, 
          options VARCHAR2,
          indexctx IN SYS.ODCIIndexCtx, scanctx IN OUT jc_idxtype_im,
          scanflg IN NUMBER)
RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS 
    res VARCHAR2(32767);
BEGIN
  return Exec_FuncB('DISSIMILARITY', target, query, options,
                   indexctx, scanctx, scanflg);
END;
/
show errors;

CREATE OR REPLACE FUNCTION Molweight_Func(
          query VARCHAR2, indexctx IN SYS.ODCIIndexCtx, 
          scanctx IN OUT jc_idxtype_im,
          scanflg IN NUMBER) RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    jchem_core_pkg.checkTableVersionEx('JC_MOLWEIGHT', indexctx);

    IF query is null THEN
      IF indexctx.IndexInfo IS NOT NULL THEN
	return jchem_core_pkg.calc_molPropNum_from_idx(indexctx.Rid, 'molweight',
				indexctx.IndexInfo.IndexSchema,
				indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
				indexctx.IndexInfo.IndexCols(1).TableSchema,
			     	indexctx.IndexInfo.IndexCols(1).TableName);
      ELSE 
	return NULL;
      END IF;
    ELSE
      IF indexctx.IndexInfo IS NOT NULL THEN
        return jchem_core_pkg.calc_molPropNum( query, 'molweight', indexctx.Rid, 
				indexctx.IndexInfo.IndexSchema, indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
				indexctx.IndexInfo.IndexCols(1).TableSchema,
			     	indexctx.IndexInfo.IndexCols(1).TableName);
      ELSE
        return jchem_core_pkg.calc_molPropNum(query, 'molweight', null, 
					   null, null, null, null, null);
      END IF;
    END IF;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION Molweight_FuncC(query CLOB, indexctx IN SYS.ODCIIndexCtx, 
                                scanctx IN OUT jc_idxtype_im,
                                scanflg IN NUMBER)
RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    jchem_core_pkg.checkTableVersionEx('JC_MOLWEIGHT', indexctx);

    IF query is null THEN
      IF indexctx.IndexInfo IS NOT NULL THEN
	    return jchem_core_pkg.calc_molPropNum_from_idx(indexctx.Rid,
                    'molweight',
                    indexctx.IndexInfo.IndexSchema,
                    indexctx.IndexInfo.IndexName,
				    indexctx.IndexInfo.IndexPartition,
                    indexctx.IndexInfo.IndexCols(1).TableSchema,
			     	indexctx.IndexInfo.IndexCols(1).TableName);
      ELSE 
        return NULL;
      END IF;
    ELSE
      IF indexctx.IndexInfo IS NOT NULL THEN
        return jchem_clob_pkg.get_molweight(
                    query, indexctx.Rid, 
                    indexctx.IndexInfo.IndexSchema,
                    indexctx.IndexInfo.IndexName,
				    indexctx.IndexInfo.IndexPartition,
                    indexctx.IndexInfo.IndexCols(1).TableSchema,
                    indexctx.IndexInfo.IndexCols(1).TableName);
      ELSE
        return jchem_clob_pkg.get_molweight(
                     query, null, null, null, null, null, null);
      END IF;
    END IF;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION Molweight_FuncB(query BLOB, indexctx IN SYS.ODCIIndexCtx, 
                                scanctx IN OUT jc_idxtype_im,
                                scanflg IN NUMBER)
RETURN NUMBER AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    jchem_core_pkg.checkTableVersionEx('JC_MOLWEIGHT', indexctx);

    IF query is null THEN
      IF indexctx.IndexInfo IS NOT NULL THEN
	    return jchem_core_pkg.calc_molPropNum_from_idx(indexctx.Rid,
                    'molweight',
                    indexctx.IndexInfo.IndexSchema,
                    indexctx.IndexInfo.IndexName,
				    indexctx.IndexInfo.IndexPartition,
                    indexctx.IndexInfo.IndexCols(1).TableSchema,
			     	indexctx.IndexInfo.IndexCols(1).TableName);
      ELSE 
        return NULL;
      END IF;
    ELSE
      IF indexctx.IndexInfo IS NOT NULL THEN
        return jchem_blob_pkg.get_molweight(
                    query, indexctx.Rid, 
                    indexctx.IndexInfo.IndexSchema,
                    indexctx.IndexInfo.IndexName,
				    indexctx.IndexInfo.IndexPartition,
                    indexctx.IndexInfo.IndexCols(1).TableSchema,
                    indexctx.IndexInfo.IndexCols(1).TableName);
      ELSE
        return jchem_blob_pkg.get_molweight(
                     query, null, null, null, null, null, null);
      END IF;
    END IF;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION Formula_Func(query VARCHAR2, indexctx IN SYS.ODCIIndexCtx, 
                                scanctx IN OUT jc_idxtype_im,
                                scanflg IN NUMBER)
RETURN VARCHAR2 AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    jchem_core_pkg.checkTableVersionEx('JC_FORMULA', indexctx);
    
    IF query is null THEN
      IF indexctx.IndexInfo IS NOT NULL THEN
	return jchem_core_pkg.calc_molProp_from_idx(indexctx.Rid, 'molformula', 
				indexctx.IndexInfo.IndexSchema,
				indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
				indexctx.IndexInfo.IndexCols(1).TableSchema,
				indexctx.IndexInfo.IndexCols(1).TableName);
      ELSE 
	return NULL;
      END IF;
    ELSE
      IF indexctx.IndexInfo IS NOT NULL THEN
        return jchem_core_pkg.calc_molProp(
                                query, 'molformula', indexctx.Rid,
				indexctx.IndexInfo.IndexSchema,
				indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
				indexctx.IndexInfo.IndexCols(1).TableSchema,
				indexctx.IndexInfo.IndexCols(1).TableName);
      ELSE
        return jchem_core_pkg.calc_molProp(
                       query, 'molformula', null, 
					   null, null, null, null, null);
      END IF;
    END IF;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION Formula_FuncC(query CLOB, indexctx IN SYS.ODCIIndexCtx, 
                                scanctx IN OUT jc_idxtype_im,
                                scanflg IN NUMBER)
RETURN VARCHAR2 AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    jchem_core_pkg.checkTableVersionEx('JC_FORMULAB', indexctx);
    
    IF query is null THEN
      IF indexctx.IndexInfo IS NOT NULL THEN
	return jchem_core_pkg.calc_molProp_from_idx(indexctx.Rid,
                'molformula', 
				indexctx.IndexInfo.IndexSchema,
				indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
				indexctx.IndexInfo.IndexCols(1).TableSchema,
				indexctx.IndexInfo.IndexCols(1).TableName);
      ELSE 
	return NULL;
      END IF;
    ELSE
      IF indexctx.IndexInfo IS NOT NULL THEN
        return jchem_clob_pkg.get_molformula(
                    query, indexctx.Rid, 
                    indexctx.IndexInfo.IndexSchema,
                    indexctx.IndexInfo.IndexName,
				    indexctx.IndexInfo.IndexPartition,
                    indexctx.IndexInfo.IndexCols(1).TableSchema,
                    indexctx.IndexInfo.IndexCols(1).TableName);
      ELSE
        return jchem_clob_pkg.get_molformula(
                     query, null, null, null, null, null, null);
      END IF;
    END IF;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION Formula_FuncB(query BLOB, indexctx IN SYS.ODCIIndexCtx, 
                                scanctx IN OUT jc_idxtype_im,
                                scanflg IN NUMBER)
RETURN VARCHAR2 AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    jchem_core_pkg.checkTableVersionEx('JC_FORMULAB', indexctx);
    
    IF query is null THEN
      IF indexctx.IndexInfo IS NOT NULL THEN
	return jchem_core_pkg.calc_molProp_from_idx(indexctx.Rid,
                'molformula', 
				indexctx.IndexInfo.IndexSchema,
				indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
				indexctx.IndexInfo.IndexCols(1).TableSchema,
				indexctx.IndexInfo.IndexCols(1).TableName);
      ELSE 
	return NULL;
      END IF;
    ELSE
      IF indexctx.IndexInfo IS NOT NULL THEN
        return jchem_blob_pkg.get_molformula(
                    query, indexctx.Rid, 
                    indexctx.IndexInfo.IndexSchema,
                    indexctx.IndexInfo.IndexName,
				    indexctx.IndexInfo.IndexPartition,
                    indexctx.IndexInfo.IndexCols(1).TableSchema,
                    indexctx.IndexInfo.IndexCols(1).TableName);
      ELSE
        return jchem_blob_pkg.get_molformula(
                     query, null, null, null, null, null, null);
      END IF;
    END IF;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION MolconvertV_Func(query VARCHAR2, type VARCHAR2,
                                indexctx IN SYS.ODCIIndexCtx,
                                scanctx IN OUT jc_idxtype_im,
                                scanflg IN NUMBER)
RETURN VARCHAR2 AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    jchem_core_pkg.checkTableVersionEx('JC_MOLCONVERTV', indexctx);

    IF query is null THEN
      IF indexctx.IndexInfo IS NOT NULL THEN
    return jchem_core_pkg.molconvertv_from_idx(
                indexctx.Rid, type,
                indexctx.IndexInfo.IndexSchema,
                indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
                indexctx.IndexInfo.IndexCols(1).TableSchema,
                    indexctx.IndexInfo.IndexCols(1).TableName);
      ELSE
    return NULL;
      END IF;
    ELSE
      IF indexctx.IndexInfo IS NOT NULL THEN
        return jchem_core_pkg.molconvertv(
                query, null, type, null, indexctx.Rid,
                indexctx.IndexInfo.IndexSchema,
                indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
                indexctx.IndexInfo.IndexCols(1).TableSchema,
                    indexctx.IndexInfo.IndexCols(1).TableName);
      ELSE
        return jchem_core_pkg.molconvertv(
                    query, null, type, null, null, null, null, null, null, null);
      END IF;
    END IF;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION MolconvertV_FuncC(query CLOB, type VARCHAR2,
                                indexctx IN SYS.ODCIIndexCtx,
                                scanctx IN OUT jc_idxtype_im,
                                scanflg IN NUMBER)
RETURN VARCHAR2 AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    jchem_core_pkg.checkTableVersionEx('JC_MOLCONVERTV', indexctx);

    IF query is null THEN
      IF indexctx.IndexInfo IS NOT NULL THEN
    return jchem_core_pkg.molconvertv_from_idx(
                indexctx.Rid, type,
                indexctx.IndexInfo.IndexSchema,
                indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
                indexctx.IndexInfo.IndexCols(1).TableSchema,
                    indexctx.IndexInfo.IndexCols(1).TableName);
      ELSE
    return NULL;
      END IF;
    ELSE
      IF indexctx.IndexInfo IS NOT NULL THEN
        return jchem_core_pkg.molconvertv(
                query, null, type, null, indexctx.Rid,
                indexctx.IndexInfo.IndexSchema,
                indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
                indexctx.IndexInfo.IndexCols(1).TableSchema,
                    indexctx.IndexInfo.IndexCols(1).TableName);
      ELSE
        return jchem_core_pkg.molconvertv(
                    query, null, type, null, null, null, null, null, null, null);
      END IF;
    END IF;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION MolconvertV_FuncB(query BLOB, type VARCHAR2,
                                indexctx IN SYS.ODCIIndexCtx,
                                scanctx IN OUT jc_idxtype_im,
                                scanflg IN NUMBER)
RETURN VARCHAR2 AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    jchem_core_pkg.checkTableVersionEx('JC_MOLCONVERTV', indexctx);

    IF query is null THEN
      IF indexctx.IndexInfo IS NOT NULL THEN
    return jchem_core_pkg.molconvertv_from_idx(
                indexctx.Rid, type,
                indexctx.IndexInfo.IndexSchema,
                indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
                indexctx.IndexInfo.IndexCols(1).TableSchema,
                    indexctx.IndexInfo.IndexCols(1).TableName);
      ELSE
    return NULL;
      END IF;
    ELSE
      IF indexctx.IndexInfo IS NOT NULL THEN
        return jchem_core_pkg.molconvertv(
                query, null, type, null, indexctx.Rid,
                indexctx.IndexInfo.IndexSchema,
                indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
                indexctx.IndexInfo.IndexCols(1).TableSchema,
                    indexctx.IndexInfo.IndexCols(1).TableName);
      ELSE
        return jchem_core_pkg.molconvertv(
                    query, null, type, null, null, null, null, null, null, null);
      END IF;
    END IF;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION MolconvertC_Func(query VARCHAR2, type VARCHAR2,
                                indexctx IN SYS.ODCIIndexCtx,
                                scanctx IN OUT jc_idxtype_im,
                                scanflg IN NUMBER)
RETURN CLOB AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    jchem_core_pkg.checkTableVersionEx('JC_MOLCONVERTC', indexctx);

    IF query is null THEN
      IF indexctx.IndexInfo IS NOT NULL THEN
        return jchem_clob_pkg.molconvertc_from_idx(
                indexctx.Rid, type,
                indexctx.IndexInfo.IndexSchema,
                indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
                indexctx.IndexInfo.IndexCols(1).TableSchema,
                indexctx.IndexInfo.IndexCols(1).TableName, null);
      ELSE
        return NULL;
      END IF;
    ELSE
      IF indexctx.IndexInfo IS NOT NULL THEN
        return jchem_clob_pkg.molconvertc(
                query, null, type, null, indexctx.Rid,
                indexctx.IndexInfo.IndexSchema,
                indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
                indexctx.IndexInfo.IndexCols(1).TableSchema,
                    indexctx.IndexInfo.IndexCols(1).TableName, null);
      ELSE
        return jchem_clob_pkg.molconvertc(
                    query, null, type, null, null, null, null, null, null, null, null);
      END IF;
    END IF;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION MolconvertC_FuncC(query CLOB, type VARCHAR2,
                                indexctx IN SYS.ODCIIndexCtx,
                                scanctx IN OUT jc_idxtype_im,
                                scanflg IN NUMBER)
RETURN CLOB AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    jchem_core_pkg.checkTableVersionEx('JC_MOLCONVERTC', indexctx);

    IF query is null THEN
      IF indexctx.IndexInfo IS NOT NULL THEN
        return jchem_clob_pkg.molconvertc_from_idx(
                indexctx.Rid, type,
                indexctx.IndexInfo.IndexSchema,
                indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
                indexctx.IndexInfo.IndexCols(1).TableSchema,
                indexctx.IndexInfo.IndexCols(1).TableName, null);
      ELSE
        return NULL;
      END IF;
    ELSE
      IF indexctx.IndexInfo IS NOT NULL THEN
        return jchem_clob_pkg.molconvertc(
                query, null, type, null, indexctx.Rid,
                indexctx.IndexInfo.IndexSchema,
                indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
                indexctx.IndexInfo.IndexCols(1).TableSchema,
                    indexctx.IndexInfo.IndexCols(1).TableName, null);
      ELSE
        return jchem_clob_pkg.molconvertc(
                    query, null, type, null, null, null, null, null, null, null, null);
      END IF;
    END IF;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION MolconvertC_FuncB(query BLOB, type VARCHAR2,
                                indexctx IN SYS.ODCIIndexCtx,
                                scanctx IN OUT jc_idxtype_im,
                                scanflg IN NUMBER)
RETURN CLOB AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    jchem_core_pkg.checkTableVersionEx('JC_MOLCONVERTC', indexctx);

    IF query is null THEN
      IF indexctx.IndexInfo IS NOT NULL THEN
        return jchem_clob_pkg.molconvertc_from_idx(
                indexctx.Rid, type,
                indexctx.IndexInfo.IndexSchema,
                indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
                indexctx.IndexInfo.IndexCols(1).TableSchema,
                indexctx.IndexInfo.IndexCols(1).TableName, null);
      ELSE
        return NULL;
      END IF;
    ELSE
      IF indexctx.IndexInfo IS NOT NULL THEN
        return jchem_clob_pkg.molconvertc(
                query, null, type, null, indexctx.Rid,
                indexctx.IndexInfo.IndexSchema,
                indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
                indexctx.IndexInfo.IndexCols(1).TableSchema,
                    indexctx.IndexInfo.IndexCols(1).TableName, null);
      ELSE
        return jchem_clob_pkg.molconvertc(
                    query, null, type, null, null, null, null, null, null, null, null);
      END IF;
    END IF;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION MolconvertB_Func(query VARCHAR2, type VARCHAR2,
                                indexctx IN SYS.ODCIIndexCtx,
                                scanctx IN OUT jc_idxtype_im,
                                scanflg IN NUMBER)
RETURN BLOB AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    jchem_core_pkg.checkTableVersionEx('JC_MOLCONVERTB', indexctx);

    IF query is null THEN
      IF indexctx.IndexInfo IS NOT NULL THEN
    return jchem_blob_pkg.molconvertb_from_idx(
                indexctx.Rid, type,
                indexctx.IndexInfo.IndexSchema,
                indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
                indexctx.IndexInfo.IndexCols(1).TableSchema,
                    indexctx.IndexInfo.IndexCols(1).TableName, null);
      ELSE
    return NULL;
      END IF;
    ELSE
      IF indexctx.IndexInfo IS NOT NULL THEN
        return jchem_blob_pkg.molconvertb(
                query, null, type, null, indexctx.Rid,
                indexctx.IndexInfo.IndexSchema,
                indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
                indexctx.IndexInfo.IndexCols(1).TableSchema,
                    indexctx.IndexInfo.IndexCols(1).TableName, null);
      ELSE
        return jchem_blob_pkg.molconvertb(
                    query, null, type, null, null, null, null, null, null, null, null);
      END IF;
    END IF;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION MolconvertB_FuncC(query CLOB, type VARCHAR2,
                                indexctx IN SYS.ODCIIndexCtx,
                                scanctx IN OUT jc_idxtype_im,
                                scanflg IN NUMBER)
RETURN BLOB AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    jchem_core_pkg.checkTableVersionEx('JC_MOLCONVERTB', indexctx);

    IF query is null THEN
      IF indexctx.IndexInfo IS NOT NULL THEN
        return jchem_blob_pkg.molconvertb_from_idx(
                indexctx.Rid, type,
                indexctx.IndexInfo.IndexSchema,
                indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
                indexctx.IndexInfo.IndexCols(1).TableSchema,
                indexctx.IndexInfo.IndexCols(1).TableName, null);
      ELSE
        return NULL;
      END IF;
    ELSE
      IF indexctx.IndexInfo IS NOT NULL THEN
        return jchem_blob_pkg.molconvertb(
                query, null, type, null, indexctx.Rid,
                indexctx.IndexInfo.IndexSchema,
                indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
                indexctx.IndexInfo.IndexCols(1).TableSchema,
                    indexctx.IndexInfo.IndexCols(1).TableName, null);
      ELSE
        return jchem_blob_pkg.molconvertb(
                    query, null, type, null, null, null, null, null, null, null, null);
      END IF;
    END IF;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION MolconvertB_FuncB(query BLOB, type VARCHAR2,
                                indexctx IN SYS.ODCIIndexCtx,
                                scanctx IN OUT jc_idxtype_im,
                                scanflg IN NUMBER)
RETURN BLOB AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    jchem_core_pkg.checkTableVersionEx('JC_MOLCONVERTB', indexctx);

    IF query is null THEN
      IF indexctx.IndexInfo IS NOT NULL THEN
        return jchem_blob_pkg.molconvertb_from_idx(
                indexctx.Rid, type,
                indexctx.IndexInfo.IndexSchema,
                indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
                indexctx.IndexInfo.IndexCols(1).TableSchema,
                indexctx.IndexInfo.IndexCols(1).TableName, null);
      ELSE
        return NULL;
      END IF;
    ELSE
      IF indexctx.IndexInfo IS NOT NULL THEN
        return jchem_blob_pkg.molconvertb(
                query, null, type, null, indexctx.Rid,
                indexctx.IndexInfo.IndexSchema,
                indexctx.IndexInfo.IndexName,
				indexctx.IndexInfo.IndexPartition,
                indexctx.IndexInfo.IndexCols(1).TableSchema,
                    indexctx.IndexInfo.IndexCols(1).TableName, null);
      ELSE
        return jchem_blob_pkg.molconvertb(
                    query, null, type, null, null, null, null, null, null, null, null);
      END IF;
    END IF;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION MolconvertB_FuncB_deprecated(query BLOB, type VARCHAR2,
                                indexctx IN SYS.ODCIIndexCtx,
                                scanctx IN OUT jc_idxtype_im,
                                scanflg IN NUMBER)
RETURN BLOB AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    jcart_logger.deprecated_method('JC_MOLCONVERTBB','JC_MOLCONVERTB');
    return MolconvertB_FuncB(query,type,indexctx,scanctx,scanflg);
  END;
/
show errors;

CREATE OR REPLACE FUNCTION Transform_Func(reaction VARCHAR2, reactants VARCHAR2,
                                indexctx IN SYS.ODCIIndexCtx,
                                scanctx IN OUT jc_idxtype_im,
                                scanflg IN NUMBER)
RETURN VARCHAR2 AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
	jcart_logger.deprecated_method('JC_TRANSFORM','JC_REACT');
    IF reaction IS NULL OR reactants IS NULL THEN
	  return NULL;
    ELSE
        return jchem_core_pkg.react(reaction, reactants, null, null, null,
                      'method:n mappingstyle:d permuteReactants:y', 't');
    END IF;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION React_Func(reaction VARCHAR2, reactants VARCHAR2,
                           options VARCHAR2,
                           indexctx IN SYS.ODCIIndexCtx,
                           scanctx IN OUT jc_idxtype_im,
                           scanflg IN NUMBER)
RETURN VARCHAR2 AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    IF reaction IS NULL OR reactants IS NULL THEN
	  return NULL;
    ELSE
        return jchem_core_pkg.react(reaction, reactants, null, null, null,
                                    options, 't');
    END IF;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION React_FuncB(reaction BLOB, reactants BLOB,
                           options VARCHAR2, indexctx IN SYS.ODCIIndexCtx,
                           scanctx IN OUT jc_idxtype_im, scanflg IN NUMBER)
RETURN BLOB AUTHID CURRENT_USER AS
BEGIN
  IF reaction IS NULL OR reactants IS NULL THEN
    return NULL;
  ELSE
      return jchem_blob_pkg.react(reaction, reactants, null, null, null,
              options, null, null);
  END IF;
END;
/
show errors;

CREATE OR REPLACE FUNCTION React4_Func(reaction VARCHAR2, reactant1 VARCHAR2,
                           reactant2 VARCHAR2, reactant3 VARCHAR2,reactant4
                           VARCHAR2, options VARCHAR2, indexctx IN
                           SYS.ODCIIndexCtx, scanctx IN OUT jc_idxtype_im,
                           scanflg IN NUMBER)
RETURN VARCHAR2 AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    IF reaction IS NULL OR reactant1 IS NULL THEN
	  return NULL;
    ELSE
        return jchem_core_pkg.react(reaction, reactant1, reactant2, reactant3,
                                    reactant4, options, null);
    END IF;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION React4_FuncC(reaction CLOB, reactant1 CLOB,
                           reactant2 CLOB, reactant3 CLOB, reactant4 CLOB,
                           options VARCHAR2, indexctx IN SYS.ODCIIndexCtx,
                           scanctx IN OUT jc_idxtype_im, scanflg IN NUMBER)
RETURN CLOB AUTHID CURRENT_USER AS
BEGIN
  IF reaction IS NULL OR reactant1 IS NULL THEN
    return NULL;
  ELSE
      return jchem_clob_pkg.react(reaction, reactant1, reactant2, reactant3,
              reactant4, options, null,  null);
  END IF;
END;
/
show errors;

CREATE OR REPLACE FUNCTION React4_FuncB(reaction BLOB, reactant1 BLOB,
                           reactant2 BLOB, reactant3 BLOB, reactant4 BLOB,
                           options VARCHAR2, indexctx IN SYS.ODCIIndexCtx,
                           scanctx IN OUT jc_idxtype_im, scanflg IN NUMBER)
RETURN BLOB AUTHID CURRENT_USER AS
BEGIN
  IF reaction IS NULL OR reactant1 IS NULL THEN
    return NULL;
  ELSE
      return jchem_blob_pkg.react(reaction, reactant1, reactant2, reactant3,
              reactant4, options, null,  null);
  END IF;
END;
/
show errors;

CREATE OR REPLACE FUNCTION Standardize_Func(structure VARCHAR2, param VARCHAR2,
                           indexctx IN SYS.ODCIIndexCtx,
                           scanctx IN OUT jc_idxtype_im,
                           scanflg IN NUMBER)
RETURN VARCHAR2 AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    jchem_core_pkg.checkTableVersionEx('JC_STANDARDIZE', indexctx);

    IF structure IS NULL THEN
        return NULL;
    ELSE
        return jchem_core_pkg.standardize(structure, param);
    END IF;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION Standardize_FuncC(structure CLOB, param VARCHAR2,
                           indexctx IN SYS.ODCIIndexCtx,
                           scanctx IN OUT jc_idxtype_im,
                           scanflg IN NUMBER) RETURN CLOB AUTHID CURRENT_USER AS
  BEGIN
    jchem_core_pkg.checkTableVersionEx('JC_STANDARDIZEB', indexctx);

    IF structure IS NULL THEN
        return NULL;
    ELSE
        return jchem_clob_pkg.standardize(structure, param, null);
    END IF;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION Standardize_FuncB(structure BLOB, param VARCHAR2,
                           indexctx IN SYS.ODCIIndexCtx,
                           scanctx IN OUT jc_idxtype_im,
                           scanflg IN NUMBER) RETURN BLOB AUTHID CURRENT_USER AS
  BEGIN
    jchem_core_pkg.checkTableVersionEx('JC_STANDARDIZEB', indexctx);

    IF structure IS NULL THEN
        return NULL;
    ELSE
        return jchem_blob_pkg.standardize(structure, param, null);
    END IF;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION Formula_Func_eq(mol VARCHAR2, query VARCHAR2) RETURN NUMBER
    AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    IF mol is null THEN return null; END IF;
    IF query is null THEN return null; END IF;
    IF jchem_clob_pkg.get_molformula(mol, null, null, null, null, null, null) = query THEN 
      RETURN 1;
    ELSE
      RETURN 0;
    END IF;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION Formula_Func_eqC(mol CLOB, query VARCHAR2) RETURN NUMBER
    AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    IF mol is null THEN return null; END IF;
    IF query is null THEN return null; END IF;
    IF jchem_clob_pkg.get_molformula(mol, null, null, null, null, null, null)=query THEN 
      RETURN 1;
    ELSE
      RETURN 0;
    END IF;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION Formula_Func_eqB(mol BLOB, query VARCHAR2) RETURN NUMBER
    AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    IF mol is null THEN return null; END IF;
    IF query is null THEN return null; END IF;
    IF jchem_blob_pkg.get_molformula(mol, null, null, null, null, null, null)=query THEN 
      RETURN 1;
    ELSE
      RETURN 0;
    END IF;
  END;
/
show errors;

CREATE OR REPLACE FUNCTION Formula_Func_eqB_deprecated(mol BLOB, query VARCHAR2) RETURN NUMBER
    AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    jcart_logger.deprecated_method('JC_FORMULA_EQB','JC_FORMULA_EQ');
    return Formula_Func_eqB(mol,query);
  END;
/
show errors;

CREATE OR REPLACE FUNCTION Fuse_Func(first_structure VARCHAR2, second_structure VARCHAR2)
RETURN VARCHAR2 AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    IF first_structure is null THEN
	  return second_structure;
	END IF;
	IF second_structure is null THEN
	  return first_structure;
	END IF;
    return jchem_core_pkg.fuse(first_structure, second_structure, NULL);
  END;
/
show errors;

CREATE OR REPLACE FUNCTION Fuse_FuncC(first_structure CLOB, second_structure CLOB)
RETURN CLOB AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    IF first_structure is null THEN
	  return second_structure;
	END IF;
	IF second_structure is null THEN
	  return first_structure;
	END IF;
    return jchem_clob_pkg.fuse(first_structure, second_structure, NULL, NULL);
  END;
/
show errors;

CREATE OR REPLACE FUNCTION Fuse_FuncB(first_structure BLOB, second_structure BLOB)
RETURN BLOB AUTHID CURRENT_USER PARALLEL_ENABLE AS
  BEGIN
    IF first_structure is null THEN
	  return second_structure;
	END IF;
	IF second_structure is null THEN
	  return first_structure;
	END IF;
    return jchem_blob_pkg.fuse(first_structure, second_structure, NULL, NULL);
  END;
/
show errors;

CREATE OR REPLACE FUNCTION User_Def_Func(name VARCHAR2, delim VARCHAR2, 
				params VARCHAR2) RETURN VARCHAR2 authid current_user AS
  BEGIN
    jcart_logger.deprecated_method('SEND_USER_FUNC','chemical terms');
    return jchem_core_pkg.send_user_func(name, delim, params);
  END;
/
show errors;

------------------------
-- CREATE OR REPLACE OPERATORS
------------------------

CREATE OR REPLACE OPERATOR jc_contains BINDING
    (VARCHAR2, VARCHAR2) RETURN NUMBER
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Contains_Func_deprecated,
    (CLOB, CLOB) RETURN NUMBER
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Contains_FuncC_deprecated,
    (BLOB, BLOB) RETURN NUMBER
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Contains_FuncB_deprecated
/

CREATE OR REPLACE OPERATOR jc_equals BINDING
    (VARCHAR2, VARCHAR2) RETURN NUMBER 
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Equals_Func_deprecated,
    (CLOB, CLOB) RETURN NUMBER 
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Equals_FuncC_deprecated,
    (BLOB, BLOB) RETURN NUMBER 
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Equals_FuncB_deprecated
/

CREATE OR REPLACE OPERATOR jc_compare BINDING
    (VARCHAR2, VARCHAR2, VARCHAR2) RETURN NUMBER 
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Compare_Func,
    (CLOB, CLOB, VARCHAR2) RETURN NUMBER 
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Compare_FuncC,
    (CLOB, VARCHAR2, VARCHAR2) RETURN NUMBER 
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Compare_FuncCV,
    (CLOB, BLOB, VARCHAR2) RETURN NUMBER 
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Compare_FuncCB,
    (BLOB, BLOB, VARCHAR2) RETURN NUMBER 
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Compare_FuncB,
    (BLOB, VARCHAR2, VARCHAR2) RETURN NUMBER 
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Compare_FuncBV
/

CREATE OR REPLACE OPERATOR jc_compare_vb BINDING(VARCHAR2, BLOB, VARCHAR2) RETURN NUMBER 
  WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
  USING compare_func_bquery
/

CREATE OR REPLACE OPERATOR jc_evaluate BINDING
    (VARCHAR2, VARCHAR2) RETURN NUMBER 
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Evaluate_Func,
    (CLOB, VARCHAR2) RETURN NUMBER 
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Evaluate_FuncC,
    (BLOB, VARCHAR2) RETURN NUMBER 
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Evaluate_FuncB,
    (VARCHAR2, VARCHAR2, VARCHAR2) RETURN NUMBER 
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Evaluate3_Func,
    (CLOB, VARCHAR2, VARCHAR2) RETURN NUMBER 
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Evaluate3_FuncC,
    (BLOB, VARCHAR2, VARCHAR2) RETURN NUMBER 
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Evaluate3_FuncB
/

CREATE OR REPLACE OPERATOR jc_evaluate_x BINDING
    (VARCHAR2, VARCHAR2) RETURN VARCHAR2
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING EvaluateX_Func,
    (CLOB, VARCHAR2) RETURN CLOB
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING EvaluateX_FuncC,
    (BLOB, VARCHAR2) RETURN BLOB
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING EvaluateX_FuncB
/

CREATE OR REPLACE OPERATOR jc_matchcount BINDING
    (VARCHAR2, VARCHAR2) RETURN NUMBER 
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING MatchCount_Func,
    (CLOB, CLOB) RETURN NUMBER 
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING MatchCount_FuncC,
    (BLOB, BLOB) RETURN NUMBER 
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING MatchCount_FuncB
/

CREATE OR REPLACE OPERATOR jc_tanimoto BINDING
    (VARCHAR2, VARCHAR2)
      RETURN NUMBER WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im 
      USING Tanimoto_Func,
    (CLOB, CLOB)
      RETURN NUMBER WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im 
      USING Tanimoto_FuncC,
    (BLOB, BLOB)
      RETURN NUMBER WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im 
      USING Tanimoto_FuncB
/

CREATE OR REPLACE OPERATOR jc_tversky BINDING
    (VARCHAR2, VARCHAR2, NUMBER, NUMBER)
      RETURN NUMBER WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im 
      USING Tversky_Func,
    (CLOB, CLOB, NUMBER, NUMBER)
      RETURN NUMBER WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im 
      USING Tversky_FuncC,
    (BLOB, BLOB, NUMBER, NUMBER)
      RETURN NUMBER WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im 
      USING Tversky_FuncB
/

CREATE OR REPLACE OPERATOR jc_dissimilarity BINDING
    (VARCHAR2, VARCHAR2)
      RETURN NUMBER WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im 
      USING Dissimilarity_Func,
    (CLOB, CLOB)
      RETURN NUMBER WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im 
      USING Dissimilarity_FuncC,
    (BLOB, BLOB)
      RETURN NUMBER WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im 
      USING Dissimilarity_FuncB,
    (VARCHAR2, VARCHAR2, VARCHAR2)
      RETURN NUMBER WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im 
      USING DissimilarityX_Func,
    (CLOB, CLOB, VARCHAR2)
      RETURN NUMBER WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im 
      USING DissimilarityX_FuncC,
    (BLOB, BLOB, VARCHAR2)
      RETURN NUMBER WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im 
      USING DissimilarityX_FuncB

/

CREATE OR REPLACE OPERATOR jc_molweight BINDING
    (VARCHAR2) RETURN NUMBER
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Molweight_Func,
    (CLOB) RETURN NUMBER
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Molweight_FuncC,
    (BLOB) RETURN NUMBER
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Molweight_FuncB
/

CREATE OR REPLACE OPERATOR jc_formula BINDING
    (VARCHAR2) RETURN VARCHAR2
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Formula_Func,
    (CLOB) RETURN VARCHAR2
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Formula_FuncC,
    (BLOB) RETURN VARCHAR2
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Formula_FuncB
/

CREATE OR REPLACE OPERATOR jc_molconvert BINDING
    (VARCHAR2, VARCHAR2) RETURN VARCHAR2
  WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
  USING MolconvertV_Func,
    (CLOB, VARCHAR2) RETURN CLOB
  WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
  USING MolconvertC_FuncC,
    (BLOB, VARCHAR2) RETURN BLOB
  WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
  USING MolconvertB_FuncB
/

CREATE OR REPLACE OPERATOR jc_molconvertv BINDING
    (VARCHAR2, VARCHAR2) RETURN VARCHAR2
  WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
  USING MolconvertV_Func,
    (CLOB, VARCHAR2) RETURN VARCHAR2
  WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
  USING MolconvertV_FuncC,
    (BLOB, VARCHAR2) RETURN VARCHAR2
  WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
  USING MolconvertV_FuncB
/

CREATE OR REPLACE OPERATOR jc_molconvertc BINDING
    (VARCHAR2, VARCHAR2) RETURN CLOB
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING MolconvertC_Func,
    (CLOB, VARCHAR2) RETURN CLOB
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING MolconvertC_FuncC,
    (BLOB, VARCHAR2) RETURN CLOB
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING MolconvertC_FuncB
/

CREATE OR REPLACE OPERATOR jc_molconvertb BINDING
    (VARCHAR2, VARCHAR2) RETURN BLOB
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING MolconvertB_Func,
    (CLOB, VARCHAR2) RETURN BLOB
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING MolconvertB_FuncC,
    (BLOB, VARCHAR2) RETURN BLOB
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING MolconvertB_FuncB
/

CREATE OR REPLACE OPERATOR jc_formula_eq BINDING
    (VARCHAR2, VARCHAR2) RETURN NUMBER 
        USING Formula_Func_eq,
    (CLOB, VARCHAR2) RETURN NUMBER 
        USING Formula_Func_eqC,
    (BLOB, VARCHAR2) RETURN NUMBER 
        USING Formula_Func_eqb
/

CREATE OR REPLACE OPERATOR jc_react BINDING(VARCHAR2, VARCHAR2, VARCHAR2) RETURN VARCHAR2
  WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
  USING React_Func
/

CREATE OR REPLACE OPERATOR jc_react4 BINDING
    (VARCHAR2, VARCHAR2, VARCHAR2, VARCHAR2, VARCHAR2, VARCHAR2) RETURN VARCHAR2
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING React4_Func,
    (CLOB, CLOB, CLOB, CLOB, CLOB, VARCHAR2) RETURN CLOB
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING React4_FuncC,
    (BLOB, BLOB, BLOB, BLOB, BLOB, VARCHAR2) RETURN BLOB
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING React4_FuncB
/

CREATE OR REPLACE OPERATOR jc_transform BINDING(VARCHAR2, VARCHAR2) RETURN VARCHAR2
  WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
  USING Transform_Func
/

CREATE OR REPLACE OPERATOR jc_standardize BINDING
    (VARCHAR2, VARCHAR2) RETURN VARCHAR2
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Standardize_Func,
    (CLOB, VARCHAR2) RETURN CLOB
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Standardize_FuncC,
    (BLOB, VARCHAR2) RETURN BLOB
      WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
      USING Standardize_FuncB
/

CREATE OR REPLACE OPERATOR jc_fuse BINDING
    (VARCHAR2, VARCHAR2) RETURN VARCHAR2
  USING Fuse_Func,
    (CLOB, CLOB) RETURN CLOB
  USING Fuse_FuncC,
    (BLOB, BLOB) RETURN BLOB
  USING Fuse_FuncB
/

CREATE OR REPLACE FUNCTION User_Def_FuncB(name VARCHAR2, delim VARCHAR2, 
				params BLOB) RETURN BLOB AUTHID CURRENT_USER AS
  BEGIN
    jcart_logger.deprecated_method('SEND_USER_FUNC','chemical terms');
    return jchem_blob_pkg.send_user_func(name, delim, params);
  END;
/
show errors;

------------------------
-- CREATE OR REPLACE OPERATORS
------------------------

CREATE OR REPLACE OPERATOR jc_containsb BINDING(BLOB, BLOB) RETURN NUMBER
  WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
  USING ContainsB_FuncB_deprecated
/

CREATE OR REPLACE OPERATOR jc_equalsb BINDING(BLOB, BLOB) RETURN NUMBER 
  WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
  USING EqualsB_FuncB_deprecated
/

CREATE OR REPLACE OPERATOR jc_compareb BINDING(BLOB, BLOB, VARCHAR2) RETURN NUMBER 
  WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
  USING Compare_FuncB
/

CREATE OR REPLACE OPERATOR jc_evaluateb BINDING(BLOB, VARCHAR2) RETURN NUMBER 
  WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
  USING Evaluate_FuncB
/

CREATE OR REPLACE OPERATOR jc_evaluateb_x BINDING(BLOB, VARCHAR2) RETURN BLOB
  WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
  USING EvaluateX_FuncB
/

CREATE OR REPLACE OPERATOR jc_matchcountb BINDING(BLOB, BLOB) RETURN NUMBER 
  WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
  USING MatchCountB_Funcb_deprecated
/

CREATE OR REPLACE OPERATOR jc_tanimotob BINDING(BLOB, BLOB)
  RETURN NUMBER WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im 
  USING TanimotoB_Funcb_deprecated
/

CREATE OR REPLACE OPERATOR jc_dissimilarityb BINDING(BLOB, BLOB)
  RETURN NUMBER WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im 
  USING DissimB_Funcb_deprecated
/

CREATE OR REPLACE OPERATOR jc_molweightb BINDING(BLOB) RETURN NUMBER
  WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
  USING Molweight_Funcb
/

CREATE OR REPLACE OPERATOR jc_formulab BINDING(BLOB) RETURN VARCHAR2
  WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
  USING Formula_Funcb
/

CREATE OR REPLACE OPERATOR jc_molconvertbb BINDING(BLOB, VARCHAR2) RETURN BLOB
  WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
  USING MolconvertB_Funcb_deprecated
/

CREATE OR REPLACE OPERATOR jc_formula_eqb BINDING(BLOB, VARCHAR2) RETURN NUMBER 
  USING Formula_Func_eqb
/

CREATE OR REPLACE OPERATOR jc_reactb BINDING(BLOB, BLOB, VARCHAR2) RETURN BLOB
  WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
  USING React_FuncB
/

CREATE OR REPLACE OPERATOR jc_reactb4 BINDING(BLOB, BLOB, BLOB, BLOB, BLOB, VARCHAR2) RETURN BLOB
  WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
  USING React4_FuncB
/

CREATE OR REPLACE OPERATOR jc_standardizeb BINDING(BLOB, VARCHAR2) RETURN BLOB
  WITH INDEX CONTEXT, SCAN CONTEXT jc_idxtype_im
  USING Standardize_FuncB
/

-------------------------------------------------------------------
-- CREATE INDEXTYPE
-------------------------------------------------------------------

-- set serveroutput on;
CALL dbms_java.set_output(2000)
/

CREATE OR REPLACE INDEXTYPE jc_idxtype
FOR
   jc_contains(VARCHAR2, VARCHAR2),
   jc_contains(CLOB, CLOB),
   jc_contains(BLOB, BLOB),
   jc_equals(VARCHAR2, VARCHAR2),
   jc_equals(CLOB, CLOB),
   jc_equals(BLOB, BLOB),
   jc_matchcount(VARCHAR2, VARCHAR2),
   jc_matchcount(CLOB, CLOB),
   jc_matchcount(BLOB, BLOB),
   jc_dissimilarity(VARCHAR2, VARCHAR2),
   jc_dissimilarity(CLOB, CLOB),
   jc_dissimilarity(BLOB, BLOB),
   jc_dissimilarity(VARCHAR2, VARCHAR2, VARCHAR2),
   jc_dissimilarity(CLOB, CLOB, VARCHAR2),
   jc_dissimilarity(BLOB, BLOB, VARCHAR2),
   jc_tanimoto(VARCHAR2, VARCHAR2),
   jc_tanimoto(CLOB, CLOB),
   jc_tanimoto(BLOB, BLOB),
   jc_tversky(VARCHAR2, VARCHAR2, NUMBER, NUMBER),
   jc_tversky(CLOB, CLOB, NUMBER, NUMBER),
   jc_tversky(BLOB, BLOB, NUMBER, NUMBER),
   jc_compare(VARCHAR2, VARCHAR2, VARCHAR2),
   jc_compare(CLOB, VARCHAR2, VARCHAR2),
   jc_compare(CLOB, CLOB, VARCHAR2),
   jc_compare(CLOB, BLOB, VARCHAR2),
   jc_compare(BLOB, BLOB, VARCHAR2),
   jc_compare(BLOB, VARCHAR2, VARCHAR2),
   jc_compare_vb(VARCHAR2, BLOB, VARCHAR2),
   jc_evaluate(VARCHAR2, VARCHAR2),
   jc_evaluate(CLOB, VARCHAR2),
   jc_evaluate(BLOB, VARCHAR2),
   jc_evaluate(VARCHAR2, VARCHAR2, VARCHAR2),
   jc_evaluate(CLOB, VARCHAR2, VARCHAR2),
   jc_evaluate(BLOB, VARCHAR2, VARCHAR2),
   jc_evaluate_x(VARCHAR2, VARCHAR2),
   jc_evaluate_x(CLOB, VARCHAR2),
   jc_evaluate_x(BLOB, VARCHAR2),
   jc_molweight(VARCHAR2),
   jc_molweight(CLOB),
   jc_molweight(BLOB),
   jc_formula(VARCHAR2),
   jc_formula(CLOB),
   jc_formula(BLOB),
   jc_formula_eq(VARCHAR2, VARCHAR2),
   jc_formula_eq(CLOB, VARCHAR2),
   jc_formula_eq(BLOB, VARCHAR2),
   jc_molconvert(VARCHAR2, VARCHAR2),
   jc_molconvert(CLOB, VARCHAR2),
   jc_molconvert(BLOB, VARCHAR2),
   jc_molconvertv(VARCHAR2, VARCHAR2),
   jc_molconvertv(CLOB, VARCHAR2),
   jc_molconvertv(BLOB, VARCHAR2),
   jc_molconvertc(VARCHAR2, VARCHAR2),
   jc_molconvertc(CLOB, VARCHAR2),
   jc_molconvertc(BLOB, VARCHAR2),
   jc_molconvertb(VARCHAR2, VARCHAR2),
   jc_molconvertb(CLOB, VARCHAR2),
   jc_molconvertb(BLOB, VARCHAR2),
   jc_react(VARCHAR2, VARCHAR2, VARCHAR2),
   jc_react4 (VARCHAR2, VARCHAR2, VARCHAR2, VARCHAR2, VARCHAR2, VARCHAR2),
   jc_react4 (CLOB, CLOB, CLOB, CLOB, CLOB, VARCHAR2),
   jc_react4 (BLOB, BLOB, BLOB, BLOB, BLOB, VARCHAR2),
   jc_standardize(VARCHAR2, VARCHAR2),
   jc_standardize(CLOB, VARCHAR2),
   jc_standardize(BLOB, VARCHAR2),
   jc_transform(VARCHAR2, VARCHAR2),
   jc_containsb(BLOB, BLOB),
   jc_equalsb(BLOB, BLOB),
   jc_matchcountb(BLOB, BLOB),
   jc_dissimilarityb(BLOB, BLOB),
   jc_tanimotob(BLOB, BLOB),
   jc_compareb(BLOB, BLOB, VARCHAR2),
   jc_evaluateb(BLOB, VARCHAR2),
   jc_evaluateb_x(BLOB, VARCHAR2),
   jc_molweightb(BLOB),
   jc_formulab(BLOB),
   jc_formula_eqb(BLOB, VARCHAR2),
   jc_molconvertbb(BLOB, VARCHAR2),
   jc_reactb (BLOB, BLOB, VARCHAR2),
   jc_reactb4 (BLOB, BLOB, BLOB, BLOB, BLOB, VARCHAR2),
   jc_standardizeb(BLOB, VARCHAR2),
   jc_fuse(VARCHAR2, VARCHAR2),
   jc_fuse(CLOB, CLOB),
   jc_fuse(BLOB, BLOB)
USING jc_idxtype_im
WITH LOCAL RANGE PARTITION
/


CREATE OR REPLACE FUNCTION jc_insert(str VARCHAR2, table_name VARCHAR2, 
                          jcprop_name VARCHAR2 := '',
                          duplicate_check VARCHAR2 := '',
                          halt_on_error VARCHAR2 := '',
                          options VARCHAR2 := '') RETURN CD_ID_ARRAY
        AUTHID CURRENT_USER AS
  BEGIN
    begin
      return jchem_table_pkg.insert_structure(str, table_name,
          jcprop_name, duplicate_check, halt_on_error, options);
    exception
    when others then
      if sqlcode = -29532 then
        jchem_core_pkg.handle_java_error(sqlerrm);
      else
        raise;
      end if;
    end;
  END;
/
show errors;

CREATE OR REPLACE PROCEDURE jc_update(str VARCHAR2, table_name VARCHAR2,
                           id NUMBER, jcprop_name VARCHAR2 := null,
                           options VARCHAR2 := null)
        AUTHID CURRENT_USER AS
  BEGIN
    jchem_table_pkg.update_structure(str, table_name, id, jcprop_name, options);
  END;
/
show errors;

--
-- Condition is a WHERE clause in the form of "WHERE ..."
--
CREATE OR REPLACE PROCEDURE jc_delete(table_name VARCHAR2, condition VARCHAR2,
                           jcprop_name VARCHAR2 := null) AUTHID CURRENT_USER AS
  BEGIN
    jchem_table_pkg.delete_structure(table_name, condition, jcprop_name);
  END;
/
show errors;

CREATE OR REPLACE FUNCTION jc_insertb(str BLOB, table_name VARCHAR2, 
                          jcprop_name VARCHAR2 := '', duplicate_check VARCHAR2 := '',
                          halt_on_error VARCHAR2 := '', options VARCHAR2 := '')
        RETURN CD_ID_ARRAY AUTHID CURRENT_USER AS
  BEGIN
    begin
      return jchem_table_pkg.insert_structure(str, table_name,
          jcprop_name, duplicate_check, halt_on_error, options);
    exception
    when others then
      if sqlcode = -29532 then
        jchem_core_pkg.handle_java_error(sqlerrm);
      else
        raise;
      end if;
    end;
  END;
/
show errors;

CREATE OR REPLACE PROCEDURE jc_updateB(str BLOB, table_name VARCHAR2,
                           id NUMBER, jcprop_name VARCHAR2 := null,
                           options VARCHAR2 := null) AUTHID CURRENT_USER AS
  BEGIN
    jchem_table_pkg.update_structure(str, table_name, id, jcprop_name, options);
  END;
/
show errors;


create sequence jchem_idxscan_no_sq
/

create sequence jchem_sessionid_sq
/

create or replace procedure testFS2621(dbhost VARCHAR2, dbport NUMBER,
                        dbname VARCHAR2, password VARCHAR2, webappUrl VARCHAR2)
        authid current_user AS
begin
    jchem_core_pkg.set_master_property(null, 'jchem.service.endPoint.url.1', webappUrl);
    jchem_misc_pkg.setDbCallback(dbhost, dbport, dbname);
    jchem_core_pkg.use_password(password);
    
    execute immediate 'create table nci_10k (id number, structure varchar2(4000))';
end;
/
show errors;

create or replace package jcf authid current_user as
    FUNCTION contains(target VARCHAR2, query VARCHAR2)
        RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION contains(target CLOB, query CLOB)
        RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION contains(target BLOB, query BLOB)
        RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;

    FUNCTION equals(target VARCHAR2, query VARCHAR2)
        RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION equals(target CLOB, query CLOB)
        RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION equals(target BLOB, query BLOB)
        RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;

    FUNCTION evaluate(target VARCHAR2, filter VARCHAR2)
                            RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION evaluate(target CLOB, filter VARCHAR2)
                            RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION evaluate(target BLOB, filter VARCHAR2)
                            RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;

    FUNCTION evaluate_x(target VARCHAR2, filter VARCHAR2)
          RETURN VARCHAR2 DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION evaluate_x(target CLOB, filter VARCHAR2, tempClob CLOB := null)
          RETURN CLOB DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION evaluate_x(target BLOB, filter VARCHAR2, tempBlob BLOB := null)
          RETURN BLOB DETERMINISTIC PARALLEL_ENABLE;

    FUNCTION t_evaluate(target VARCHAR2, filter VARCHAR2)
        RETURN COMP_CHAR_ARRAY DETERMINISTIC PARALLEL_ENABLE PIPELINED;
    FUNCTION t_evaluate(target CLOB, filter VARCHAR2)
        RETURN COMP_CLOB_ARRAY DETERMINISTIC PARALLEL_ENABLE PIPELINED;
    FUNCTION t_evaluate(target BLOB, filter VARCHAR2)
        RETURN COMP_BLOB_ARRAY DETERMINISTIC PARALLEL_ENABLE PIPELINED;

    FUNCTION compare(target VARCHAR2, query VARCHAR2, options VARCHAR2)
                RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION compare(target CLOB, query CLOB, options VARCHAR2)
                RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION compare(target CLOB, query VARCHAR2, options VARCHAR2)
                RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION compare(target BLOB, query BLOB, options VARCHAR2)
                RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;

    FUNCTION compare(target VARCHAR2, query BLOB, options VARCHAR2)
                RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;

    FUNCTION hitColorAndAlign(tblSchema VARCHAR2, tblName VARCHAR2,
                colName VARCHAR2, query CLOB, rowids VARCHAR2,
                options VARCHAR2, hitColorAndAlignOptions VARCHAR2)
                RETURN COMP_CLOB_ARRAY DETERMINISTIC PARALLEL_ENABLE PIPELINED;
    FUNCTION hitColorAndAlign(tblSchema VARCHAR2, tblName VARCHAR2,
                colName VARCHAR2, query BLOB, rowids VARCHAR2,
                options VARCHAR2, hitColorAndAlignOptions VARCHAR2)
                RETURN COMP_BLOB_ARRAY DETERMINISTIC PARALLEL_ENABLE PIPELINED;

    FUNCTION matchcount(target VARCHAR2, query VARCHAR2)
              RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION matchcount(target CLOB, query CLOB)
              RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION matchcount(target BLOB, query BLOB)
              RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;

    FUNCTION tanimoto(target VARCHAR2, query VARCHAR2)
              RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION tanimoto(target CLOB, query CLOB)
              RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION tanimoto(target BLOB, query BLOB)
              RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;

    FUNCTION tversky(target VARCHAR2, query VARCHAR2,
                     targetWeight number, queryWeight number)
              RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION tversky(target CLOB, query CLOB,
                     targetWeight number, queryWeight number)
              RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION tversky(target BLOB, query BLOB,
                     targetWeight number, queryWeight number)
              RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;

    FUNCTION dissimilarity(target VARCHAR2, query VARCHAR2)
                          RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION dissimilarity(target CLOB, query CLOB)
                          RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION dissimilarity(target BLOB, query BLOB)
                          RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;

    FUNCTION molweight(query VARCHAR2)
                        RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION molweight(query CLOB)
                        RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION molweight(query BLOB)
                        RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;

    FUNCTION formula(query VARCHAR2)
                RETURN VARCHAR2 DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION formula(query CLOB)
                RETURN VARCHAR2 DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION formula(query BLOB)
                RETURN VARCHAR2 DETERMINISTIC PARALLEL_ENABLE;

    FUNCTION formula_search(target VARCHAR2, query VARCHAR2, searchType VARCHAR2)
                RETURN VARCHAR2 DETERMINISTIC PARALLEL_ENABLE;

    FUNCTION molconvert(query VARCHAR2, type VARCHAR2)
                RETURN VARCHAR2 DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION molconvert(query CLOB, type VARCHAR2, tmpClob CLOB := null)
                RETURN CLOB DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION molconvert(query BLOB, type VARCHAR2, tmpBlob BLOB := null)
                RETURN BLOB DETERMINISTIC PARALLEL_ENABLE;

    FUNCTION molconvert(query VARCHAR2, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2 := null)
                RETURN VARCHAR2 DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION molconvert(query CLOB, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2 := null, tmpClob CLOB := null)
                RETURN CLOB DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION molconvert(query BLOB, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2 := null, tmpBlob BLOB := null)
                RETURN BLOB DETERMINISTIC PARALLEL_ENABLE;

    FUNCTION molconvertv(query VARCHAR2, type VARCHAR2)
                RETURN VARCHAR2 DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION molconvertv(query CLOB, type VARCHAR2)
                RETURN VARCHAR2 DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION molconvertv(query BLOB, type VARCHAR2)
                RETURN VARCHAR2 DETERMINISTIC PARALLEL_ENABLE;

    FUNCTION molconvertv(query VARCHAR2, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2 := null)
                RETURN VARCHAR2 DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION molconvertv(query CLOB, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2 := null)
                RETURN VARCHAR2 DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION molconvertv(query BLOB, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2 := null)
                RETURN VARCHAR2 DETERMINISTIC PARALLEL_ENABLE;

	FUNCTION molconvertc(query VARCHAR2, type VARCHAR2, tmpClob CLOB := null)
                RETURN CLOB DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION molconvertc(query CLOB, type VARCHAR2, tmpClob CLOB := null)
                RETURN CLOB DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION molconvertc(query BLOB, type VARCHAR2, tmpClob CLOB := null)
                RETURN CLOB DETERMINISTIC PARALLEL_ENABLE;

    FUNCTION molconvertc(query VARCHAR2, inputFormat VARCHAR2, type VARCHAR2,
                         options VARCHAR2 := null, tmpClob CLOB := null) 
				RETURN CLOB DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION molconvertc(query CLOB, inputFormat VARCHAR2, type VARCHAR2,
                         options VARCHAR2 := null, tmpClob CLOB := null) 
				RETURN CLOB DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION molconvertc(query BLOB, inputFormat VARCHAR2, type VARCHAR2,
                         options VARCHAR2 := null, tmpClob CLOB := null) 
				RETURN CLOB DETERMINISTIC PARALLEL_ENABLE;

	FUNCTION molconvertb(query VARCHAR2, type VARCHAR2, tmpBlob BLOB := null)
                RETURN BLOB DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION molconvertb(query CLOB, type VARCHAR2, tmpBlob BLOB := null)
                RETURN BLOB DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION molconvertb(query BLOB, type VARCHAR2, tmpBlob BLOB := null)
                RETURN BLOB DETERMINISTIC PARALLEL_ENABLE;

    FUNCTION molconvertb(query VARCHAR2, inputFormat VARCHAR2, type VARCHAR2,
                         options VARCHAR2 := null, tmpBlob BLOB := null) 
				RETURN BLOB DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION molconvertb(query CLOB, inputFormat VARCHAR2, type VARCHAR2,
                         options VARCHAR2 := null, tmpBlob BLOB := null) 
				RETURN BLOB DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION molconvertb(query BLOB, inputFormat VARCHAR2, type VARCHAR2,
                         options VARCHAR2 := null, tmpBlob BLOB := null) 
				RETURN BLOB DETERMINISTIC PARALLEL_ENABLE;

    FUNCTION transform(reaction VARCHAR2, reactants VARCHAR2)
              RETURN VARCHAR2 DETERMINISTIC PARALLEL_ENABLE;

    FUNCTION react(reaction VARCHAR2, reactants VARCHAR2,
                             options VARCHAR2)
           RETURN VARCHAR2 DETERMINISTIC PARALLEL_ENABLE;

    FUNCTION react4(reaction VARCHAR2, reactant1 VARCHAR2,
                reactant2 VARCHAR2, reactant3 VARCHAR2, reactant4 VARCHAR2,
                options VARCHAR2)
        RETURN VARCHAR2 DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION react4(reaction CLOB, reactant1 CLOB, reactant2 CLOB,
                    reactant3 CLOB, reactant4 CLOB, options VARCHAR2,
                    tempClob CLOB := null)
        RETURN CLOB DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION react4(reaction BLOB, reactant1 BLOB, reactant2 BLOB,
                    reactant3 BLOB, reactant4 BLOB, options VARCHAR2,
                    tempBlob BLOB := null)
        RETURN BLOB DETERMINISTIC PARALLEL_ENABLE;

    FUNCTION t_react4(reaction VARCHAR2, reactant1 VARCHAR2, reactant2
                      VARCHAR2, reactant3 VARCHAR2, reactant4 VARCHAR2,
                      options VARCHAR2)
        RETURN CHAR_PRODUCT_ARRAY DETERMINISTIC PARALLEL_ENABLE PIPELINED;
    FUNCTION t_react4(reaction CLOB, reactant1 CLOB, reactant2 CLOB, reactant3
                      CLOB, reactant4 CLOB, options VARCHAR2)
        RETURN CLOB_PRODUCT_ARRAY DETERMINISTIC PARALLEL_ENABLE PIPELINED;
    FUNCTION t_react4(reaction BLOB, reactant1 BLOB, reactant2 BLOB, reactant3
                      BLOB, reactant4 BLOB, options VARCHAR2)
        RETURN BLOB_PRODUCT_ARRAY DETERMINISTIC PARALLEL_ENABLE PIPELINED;

    FUNCTION standardize(structure VARCHAR2, param VARCHAR2)
           RETURN VARCHAR2 DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION standardize(structure CLOB, param VARCHAR2, tempClob CLOB := null)
           RETURN CLOB DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION standardize(structure BLOB, param VARCHAR2, tempBlob BLOB := null)
           RETURN BLOB DETERMINISTIC PARALLEL_ENABLE;

    FUNCTION formula_eq(mol VARCHAR2, query VARCHAR2) RETURN NUMBER
      DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION formula_eq(mol CLOB, query VARCHAR2) RETURN NUMBER
      DETERMINISTIC PARALLEL_ENABLE;
    FUNCTION formula_eq(mol BLOB, query VARCHAR2) RETURN NUMBER
      DETERMINISTIC PARALLEL_ENABLE;

    FUNCTION insertable(structure VARCHAR2, schemaname VARCHAR2, tablename VARCHAR2, columnname VARCHAR2)
           RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;
	FUNCTION insertable(structure CLOB, schemaname VARCHAR2, tablename VARCHAR2, columnname VARCHAR2)
           RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;
	FUNCTION insertable(structure BLOB, schemaname VARCHAR2, tablename VARCHAR2, columnname VARCHAR2)
           RETURN NUMBER DETERMINISTIC PARALLEL_ENABLE;
		   
	FUNCTION fuse(first_structure IN VARCHAR2, second_structure IN VARCHAR2, inputformat IN VARCHAR2 := NULL)
      RETURN VARCHAR2 DETERMINISTIC PARALLEL_ENABLE;
	FUNCTION fuse(first_structure IN CLOB, second_structure IN CLOB, inputformat IN VARCHAR2 := NULL, tmpClob CLOB := null)
      RETURN CLOB DETERMINISTIC PARALLEL_ENABLE;
	FUNCTION fuse(first_structure IN BLOB, second_structure IN BLOB, inputformat IN VARCHAR2 := NULL, tmpBlob BLOB := null)
      RETURN BLOB DETERMINISTIC PARALLEL_ENABLE;

  function exec_v(options in varchar2, execPath in varchar2, param1 in
                    varchar2, param2 in varchar2) return varchar2 DETERMINISTIC PARALLEL_ENABLE
      as language java name
        'chemaxon.jchem.cartridge.oraresident.JFunctions.exec(java.lang.String,
        java.lang.String, java.lang.String, java.lang.String)
        return java.lang.String';

  function t_get_gmem_util return gmem_util_array DETERMINISTIC PARALLEL_ENABLE pipelined;

  function t_get_taskinfo return taskinfo_array DETERMINISTIC PARALLEL_ENABLE pipelined;

end jcf;
/
show errors

create or replace package body jcf as 
  FUNCTION contains(target VARCHAR2, query VARCHAR2)
      RETURN NUMBER AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF query is null or target is null THEN 
      return null; 
    END IF;
      error := jchem_core_pkg.exec_function('CONTAINS', target, query,
                           null, null, null, null, null, null, null, null,
                           null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION contains(target CLOB, query CLOB)
      RETURN NUMBER PARALLEL_ENABLE AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF query is null or target is null THEN 
      return null; 
    END IF;
      error := jchem_clob_pkg.exec_function('CONTAINS', target, query,
                           null, null, null, null, null, null, null, null,
                           null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION contains(target BLOB, query BLOB)
      RETURN NUMBER PARALLEL_ENABLE AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF query is null or target is null THEN 
      return null; 
    END IF;
      error := jchem_blob_pkg.exec_function('CONTAINS', target, query,
                       null, null, null, null, null, null, null, null,
                       null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION equals(target VARCHAR2, query VARCHAR2)
      RETURN NUMBER PARALLEL_ENABLE AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF query is null or target is null THEN 
      return null; 
    END IF;
      error := jchem_core_pkg.exec_function('EQUALS', target, query,
                       null, null, null, null, null, null, null, null, null,
                       null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION equals(target CLOB, query CLOB)
      RETURN NUMBER PARALLEL_ENABLE AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF query is null or target is null THEN 
      return null; 
    END IF;
      error := jchem_clob_pkg.exec_function('EQUALS', target, query,
                       null, null, null, null, null, null, null, null,
                       null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION equals(target BLOB, query BLOB)
      RETURN NUMBER PARALLEL_ENABLE AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF query is null or target is null THEN 
      return null; 
    END IF;
      error := jchem_blob_pkg.exec_function('EQUALS', target, query,
                       null, null, null, null, null, null, null, null,
                       null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION evaluate(target VARCHAR2, filter VARCHAR2)
                          RETURN NUMBER PARALLEL_ENABLE AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF target is null THEN 
      return null; 
    END IF;
    error := jchem_core_pkg.exec_function('EVALUATE', target, null,
                       filter, null, null, null, null, null, null, null,
                       null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION evaluate(target CLOB, filter VARCHAR2)
        RETURN NUMBER PARALLEL_ENABLE AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF target is null THEN 
      return null; 
    END IF;
      error := jchem_clob_pkg.exec_function('EVALUATE', target, null,
                       filter, null, null, null, null, null, null, null,
                       null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION evaluate(target BLOB, filter VARCHAR2)
        RETURN NUMBER PARALLEL_ENABLE AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF target is null THEN 
      return null; 
    END IF;
      error := jchem_blob_pkg.exec_function('EVALUATE', target, null,
                      filter, null, null, null, null, null, null, null,
                      null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION evaluate_x(target VARCHAR2, filter VARCHAR2)
        RETURN VARCHAR2 PARALLEL_ENABLE AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF target is null THEN 
      return null; 
    END IF;
    error := jchem_core_pkg.exec_function('EVALUATE_X', target, null,
                      filter, null, null, null, null, null, null, null,
                      null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION evaluate_x(target CLOB, filter VARCHAR2,
        tempClob CLOB) RETURN CLOB PARALLEL_ENABLE AS
    error varchar2(32767);
    result CLOB;
  BEGIN
    IF target is null THEN 
      return null; 
    END IF;
      error := jchem_clob_pkg.exec_function__c('EVALUATE_X', target, null,
                      filter, null, null, null, null, null, null, null,
                      null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION evaluate_x(target BLOB, filter VARCHAR2, tempBlob BLOB)
          RETURN BLOB PARALLEL_ENABLE AS
    error varchar2(32767);
    result BLOB;
  BEGIN
    IF target is null THEN 
      return null; 
    END IF;
      error := jchem_blob_pkg.exec_function__b('EVALUATE_X', target, null,
                      filter, null, null, null, null, null, null, null,
                      null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION t_evaluate(target VARCHAR2, filter VARCHAR2)
      RETURN COMP_CHAR_ARRAY PARALLEL_ENABLE PIPELINED AS
        chararr CHAR_ARRAY;
  BEGIN
    IF target is null or filter is null THEN 
      return; 
    END IF;

    chararr := jchem_core_pkg.evaluate_arr(target, filter, null,
                              null, null, null, null, null);
    for i in 1..chararr.count loop
      pipe row(COMPOSITE_CHAR(chararr(i)));
    end loop;

    return;
  END;

  FUNCTION t_evaluate(target CLOB, filter VARCHAR2)
          RETURN COMP_CLOB_ARRAY PARALLEL_ENABLE PIPELINED AS
      clobarr CLOB_ARRAY;
  BEGIN
      IF target is null or filter is null THEN 
        return;
      END IF;

      clobarr := jchem_clob_pkg.evaluate_arr(target, filter,
                      null, null, null, null, null, null);
      for i in 1..clobarr.count loop
        pipe row(COMPOSITE_CLOB(clobarr(i)));
      end loop;
    
      return;
  END;

  FUNCTION t_evaluate(target BLOB, filter VARCHAR2)
          RETURN COMP_BLOB_ARRAY PARALLEL_ENABLE PIPELINED AS
      blobarr BLOB_ARRAY;
  BEGIN
      IF target is null or filter is null THEN 
        return;
      END IF;

      blobarr := jchem_blob_pkg.evaluate_arr(target, filter,
                      null, null, null, null, null, null);
      for i in 1..blobarr.count loop
        pipe row(COMPOSITE_BLOB(blobarr(i)));
      end loop;
    
      return;
  END;

  FUNCTION compare(target VARCHAR2, query VARCHAR2, options VARCHAR2)
              RETURN NUMBER PARALLEL_ENABLE AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF query is null or target is null THEN 
      return null; 
    END IF;
    error := jchem_core_pkg.exec_function('COMPARE', target, query, options,
                                    null, null, null, null, null, null, null,
                                    null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION compare(target CLOB, query CLOB, options VARCHAR2)
      RETURN NUMBER PARALLEL_ENABLE AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF query is null or target is null THEN 
      return null; 
    END IF;
    error := jchem_clob_pkg.exec_function('COMPARE', target, query, options,
                                    null, null, null, null, null, null, null,
                                    null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION compare(target CLOB, query VARCHAR2, options VARCHAR2)
      RETURN NUMBER PARALLEL_ENABLE AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF query is null or target is null THEN 
      return null; 
    END IF;
    error := jchem_clob_pkg.exec_function_cv('COMPARE', target, query, options,
                                    null, null, null, null, null, null, null,
                                    null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION compare(target CLOB, query BLOB, options VARCHAR2)
      RETURN NUMBER PARALLEL_ENABLE AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF query is null or target is null THEN 
      return null; 
    END IF;
    error := jchem_clob_pkg.exec_function_cb('COMPARE', target, query, options,
                                    null, null, null, null, null, null, null,
                                    null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION compare(target BLOB, query BLOB, options VARCHAR2)
      RETURN NUMBER PARALLEL_ENABLE AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF query is null or target is null THEN 
      return null; 
    END IF;
    error := jchem_blob_pkg.exec_function('COMPARE', target, query, options,
                                    null, null, null, null, null, null, null,
                                    null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION compare(target VARCHAR2, query BLOB, options VARCHAR2)
              RETURN NUMBER PARALLEL_ENABLE AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF query is null or target is null THEN 
        return null; 
    END IF;
    error := jchem_blob_pkg.exec_function_vb('COMPARE', target, query, options,
                                    null, null, null, null, null, null, null,
                                    null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION hitColorAndAlign(tblSchema VARCHAR2, tblName VARCHAR2,
              colName VARCHAR2, query CLOB, rowids VARCHAR2,
              options VARCHAR2, hitColorAndAlignOptions VARCHAR2)
              RETURN COMP_CLOB_ARRAY PARALLEL_ENABLE PIPELINED AS
    clobarr CLOB_ARRAY;
  BEGIN
    IF query is null THEN 
      return; 
    END IF;

    clobarr := jchem_clob_pkg.hitColorAndAlign(tblSchema, tblName,
                  colName, query, rowids, options, hitColorAndAlignOptions,
                  jchem_defright_pkg.get_next_scanId());

    for i in 1..clobarr.count loop
      pipe row(COMPOSITE_CLOB(clobarr(i)));
    end loop;

    return;
  END;

  FUNCTION hitColorAndAlign(tblSchema VARCHAR2, tblName VARCHAR2,
              colName VARCHAR2, query BLOB, rowids VARCHAR2,
              options VARCHAR2, hitColorAndAlignOptions VARCHAR2)
              RETURN COMP_BLOB_ARRAY PARALLEL_ENABLE PIPELINED AS
    blobarr BLOB_ARRAY;
  BEGIN
    IF query is null THEN 
      return; 
    END IF;

    blobarr := jchem_blob_pkg.hitColorAndAlign(tblSchema, tblName,
                  colName, query, rowids, options, hitColorAndAlignOptions,
                  jchem_defright_pkg.get_next_scanId());

    for i in 1..blobarr.count loop
      pipe row(COMPOSITE_BLOB(blobarr(i)));
    end loop;

    return;
  END;


  FUNCTION matchcount(target VARCHAR2, query VARCHAR2)
            RETURN NUMBER PARALLEL_ENABLE AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF query is null or target is null THEN 
      return null; 
    END IF;
    error := jchem_core_pkg.exec_function('MATCHCOUNT', target, query, null,
                                    null, null, null, null, null, null, null,
                                    null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION matchcount(target CLOB, query CLOB) RETURN NUMBER PARALLEL_ENABLE AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF query is null or target is null THEN 
      return null; 
    END IF;
    error := jchem_clob_pkg.exec_function('MATCHCOUNT', target, query, null,
                                    null, null, null, null, null, null, null,
                                    null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION matchcount(target BLOB, query BLOB) RETURN NUMBER PARALLEL_ENABLE AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF query is null or target is null THEN 
      return null; 
    END IF;
    error := jchem_blob_pkg.exec_function('MATCHCOUNT', target, query, null,
                                    null, null, null, null, null, null, null,
                                    null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION tanimoto(target VARCHAR2, query VARCHAR2)
            RETURN NUMBER PARALLEL_ENABLE AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF query is null or target is null THEN 
      return null; 
    END IF;
    error := jchem_core_pkg.exec_function('TANIMOTO', target, query,
                                    'dissimilarityMetric:tanimoto',
                                    null, null, null, null, null, null, null,
                                    null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION tanimoto(target CLOB, query CLOB) RETURN NUMBER PARALLEL_ENABLE AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF query is null or target is null THEN 
      return null; 
    END IF;
    error := jchem_clob_pkg.exec_function('TANIMOTO', target, query,
                                    null, 'dissimilarityMetric:tanimoto',
                                    null, null, null, null, null, null,
                                    null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION tanimoto(target BLOB, query BLOB) RETURN NUMBER PARALLEL_ENABLE AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF query is null or target is null THEN 
      return null; 
    END IF;
    error := jchem_blob_pkg.exec_function('TANIMOTO', target, query,
                                    'dissimilarityMetric:tanimoto',
                                    null, null, null, null, null, null, null,
                                    null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION tversky(target VARCHAR2, query VARCHAR2, targetWeight number, queryWeight number)
            RETURN NUMBER PARALLEL_ENABLE AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF query is null or target is null THEN 
      return null; 
    END IF;
    error := jchem_core_pkg.exec_function('TVERSKY', target, query,
                              'tversky'
                              || jchem_core_pkg.simCalcSeparator 
                              || targetWeight || jchem_core_pkg.simCalcSeparator || queryWeight,
                              null, null, null, null, null, null, null,
                              null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION tversky(target CLOB, query CLOB, targetWeight number, queryWeight number)
            RETURN NUMBER PARALLEL_ENABLE AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF query is null or target is null THEN 
      return null; 
    END IF;
    error := jchem_clob_pkg.exec_function('TVERSKY', target, query,
                              'tversky'
                              || jchem_core_pkg.simCalcSeparator 
                              || targetWeight || jchem_core_pkg.simCalcSeparator || queryWeight,
                              null, null, null, null, null, null, null,
                              null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION tversky(target BLOB, query BLOB, targetWeight number, queryWeight number)
            RETURN NUMBER PARALLEL_ENABLE AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF query is null or target is null THEN 
      return null; 
    END IF;
    error := jchem_blob_pkg.exec_function('TVERSKY', target, query,
                              'tversky'
                              || jchem_core_pkg.simCalcSeparator 
                              || targetWeight || jchem_core_pkg.simCalcSeparator || queryWeight,
                              null, null, null, null, null, null, null,
                              null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION dissimilarity(target VARCHAR2, query VARCHAR2)
                        RETURN NUMBER PARALLEL_ENABLE AS 
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF query is null or target is null THEN 
      return null; 
    END IF;
    error := jchem_core_pkg.exec_function('DISSIMILARITY', target, query,
                              null, null, null, null, null, null, null, null,
                              null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION dissimilarity(target CLOB, query CLOB) RETURN NUMBER PARALLEL_ENABLE AS 
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF query is null or target is null THEN 
      return null; 
    END IF;
    error := jchem_clob_pkg.exec_function('DISSIMILARITY', target, query,
                              null, null, null, null, null, null, null, null,
                              null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION dissimilarity(target BLOB, query BLOB) RETURN NUMBER PARALLEL_ENABLE AS 
      res VARCHAR2(32767);
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF query is null or target is null THEN 
      return null; 
    END IF;
    error := jchem_blob_pkg.exec_function('DISSIMILARITY', target, query,
                              null, null, null, null, null, null, null, null,
                              null, null, null, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION molweight(query VARCHAR2) RETURN NUMBER PARALLEL_ENABLE AS
  BEGIN
    IF query is null THEN 
      return null; 
    END IF;
    return jchem_core_pkg.calc_molPropNum(query, 'molweight', null, null,
                                              null, null, null, null);
  END;

  FUNCTION molweight(query CLOB) RETURN NUMBER PARALLEL_ENABLE AS
  BEGIN
    IF query is null THEN 
      return null; 
    END IF;
    return jchem_clob_pkg.get_molweight(
        query, null, null, null, null, null, null);
  END;

  FUNCTION molweight(query BLOB) RETURN NUMBER PARALLEL_ENABLE AS
  BEGIN
    IF query is null THEN 
      return null; 
    END IF;
    return jchem_blob_pkg.get_molweight(
        query, null, null, null, null, null, null);
  END;

  FUNCTION formula(query VARCHAR2)
              RETURN VARCHAR2 PARALLEL_ENABLE AS
  BEGIN
    IF query is null THEN 
      return null; 
    END IF;
    return jchem_core_pkg.calc_molProp(query, 'molformula', null, 
                                       null, null, null, null, null);
  END;

  FUNCTION formula(query CLOB) RETURN VARCHAR2 PARALLEL_ENABLE AS
  BEGIN
      IF query is null THEN 
        return null; 
      END IF;
      return jchem_clob_pkg.get_molformula(query,
              null, null, null, null, null, null);
  END;

  FUNCTION formula(query BLOB) RETURN VARCHAR2 PARALLEL_ENABLE AS
  BEGIN
      IF query is null THEN 
        return null; 
      END IF;
      return jchem_blob_pkg.get_molformula(query,
              null, null, null, null, null, null);
  END;

  FUNCTION formula_search(target VARCHAR2, query VARCHAR2, searchType VARCHAR2)
              RETURN VARCHAR2 PARALLEL_ENABLE AS
    error varchar2(32767);
    result varchar2(32767);
  BEGIN
    IF target is null or query is null THEN 
      return null; 
    END IF;
    error := jchem_core_pkg.formula_search(target, query, searchType, result);
    if error is null then
      return result;
    else
      jchem_core_pkg.handle_java_error(error);
    end if;
  END;

  FUNCTION molconvert(query VARCHAR2, type VARCHAR2)
              RETURN VARCHAR2 PARALLEL_ENABLE AS
  BEGIN
      return molconvertv(query, type);
  END;

  FUNCTION molconvert(query CLOB, type VARCHAR2, tmpClob CLOB := null)
          RETURN CLOB PARALLEL_ENABLE AS
  BEGIN
      return molconvertc(query, type, tmpClob);
  END;

  FUNCTION molconvert(query BLOB, type VARCHAR2, tmpBlob BLOB := null)
          RETURN BLOB PARALLEL_ENABLE AS
  BEGIN
      return molconvertb(query, type, tmpBlob);
  END;

  FUNCTION molconvert(query VARCHAR2, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2 := null)
              RETURN VARCHAR2 PARALLEL_ENABLE AS
  BEGIN
      return molconvertv(query, inputFormat, type, options);
  END;

  FUNCTION molconvert(query CLOB, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2 := null, tmpClob CLOB := null)
          RETURN CLOB PARALLEL_ENABLE AS
  BEGIN
      return molconvertc(query, inputFormat, type, options, tmpClob);
  END;

  FUNCTION molconvert(query BLOB, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2 := null, tmpBlob BLOB := null)
          RETURN BLOB PARALLEL_ENABLE AS
  BEGIN
	  return molconvertb(query, inputFormat, type, options, tmpBlob);
  END;

  FUNCTION molconvertv(query VARCHAR2, type VARCHAR2)
              RETURN VARCHAR2 PARALLEL_ENABLE AS
  BEGIN
      return molconvertv(query, null, type, null);
  END;

  FUNCTION molconvertv(query CLOB, type VARCHAR2)
          RETURN VARCHAR2 PARALLEL_ENABLE AS
  BEGIN
      return molconvertv(query, null, type, null);
  END;

  FUNCTION molconvertv(query BLOB, type VARCHAR2)
          RETURN VARCHAR2 PARALLEL_ENABLE AS
  BEGIN
      return molconvertv(query, null, type, null);
  END;

  FUNCTION molconvertv(query VARCHAR2, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2 := null)
              RETURN VARCHAR2 PARALLEL_ENABLE AS
  BEGIN
      IF query is null THEN
          return null; 
      END IF;
      return jchem_core_pkg.molconvertv(query, inputFormat, type, options, null,
                          null, null, null, null, null);
  END;

  FUNCTION molconvertv(query CLOB, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2 := null)
          RETURN VARCHAR2 PARALLEL_ENABLE AS
  BEGIN
      IF query is null THEN
          return null; 
      END IF;
      return jchem_core_pkg.molconvertv(query, inputFormat, type, options, null,
                          null, null, null, null, null);
  END;

  FUNCTION molconvertv(query BLOB, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2 := null)
          RETURN VARCHAR2 PARALLEL_ENABLE AS
  BEGIN
      IF query is null THEN
          return null; 
      END IF;
      return jchem_core_pkg.molconvertv(query, inputFormat, type, options, null,
                          null, null, null, null, null);
  END;

  FUNCTION molconvertc(query VARCHAR2, type VARCHAR2, tmpClob CLOB := null)
              RETURN CLOB PARALLEL_ENABLE AS
  BEGIN
    return molconvertc(query, null, type, null, tmpClob);
  END;

  FUNCTION molconvertc(query CLOB, type VARCHAR2, tmpClob CLOB := null)
        RETURN CLOB PARALLEL_ENABLE AS
  BEGIN
    return molconvertc(query, null, type, null, tmpClob);
  END;

  FUNCTION molconvertc(query BLOB, type VARCHAR2, tmpClob CLOB := null)
        RETURN CLOB PARALLEL_ENABLE AS
  BEGIN
    return molconvertc(query, null, type, null, tmpClob);
  END;

  FUNCTION molconvertc(query VARCHAR2, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2 := null, tmpClob CLOB := null)
              RETURN CLOB PARALLEL_ENABLE AS
  BEGIN
    IF query is null THEN
        return NULL;
    END IF;
    return jchem_clob_pkg.molconvertc(
          query, inputFormat, type, options, null, null, null, null, null, null, tmpClob);
  END;

  FUNCTION molconvertc(query CLOB, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2 := null, tmpClob CLOB := null)
        RETURN CLOB PARALLEL_ENABLE AS
  BEGIN
      IF query is null THEN
          return NULL;
      END IF;
      return jchem_clob_pkg.molconvertc(
          query, inputFormat, type, options, null, null, null, null, null, null, tmpClob);
  END;

  FUNCTION molconvertc(query BLOB, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2 := null, tmpClob CLOB := null)
        RETURN CLOB PARALLEL_ENABLE AS
  BEGIN
      IF query is null THEN
          return NULL;
      END IF;
      return jchem_clob_pkg.molconvertc(
          query, inputFormat, type, options, null, null, null, null, null, null, tmpClob);
  END;

  FUNCTION molconvertb(query VARCHAR2, type VARCHAR2, tmpBlob BLOB := null)
              RETURN BLOB PARALLEL_ENABLE AS
  BEGIN
    return molconvertb(query, null, type, null, tmpBlob);
  END;

  FUNCTION molconvertb(query CLOB, type VARCHAR2, tmpBlob BLOB := null)
        RETURN BLOB PARALLEL_ENABLE AS
  BEGIN
    return molconvertb(query, null, type, null, tmpBlob);
  END;

  FUNCTION molconvertb(query BLOB, type VARCHAR2, tmpBlob BLOB := null)
        RETURN BLOB PARALLEL_ENABLE AS
  BEGIN
    return molconvertb(query, null, type, null, tmpBlob);
  END;

  FUNCTION molconvertb(query VARCHAR2, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2 := null, tmpBlob BLOB := null)
              RETURN BLOB PARALLEL_ENABLE AS
  BEGIN
    IF query is null THEN
        return NULL;
    END IF;
    return jchem_blob_pkg.molconvertb(
          query, inputFormat, type, options, null, null, null, null, null, null, tmpBlob);
  END;

  FUNCTION molconvertb(query CLOB, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2 := null, tmpBlob BLOB := null)
        RETURN BLOB PARALLEL_ENABLE AS
  BEGIN
      IF query is null THEN
          return NULL;
      END IF;
      return jchem_blob_pkg.molconvertb(
          query, inputFormat, type, options, null, null, null, null, null, null, tmpBlob);
  END;

  FUNCTION molconvertb(query BLOB, inputFormat VARCHAR2, type VARCHAR2, options VARCHAR2 := null, tmpBlob BLOB := null)
        RETURN BLOB PARALLEL_ENABLE AS
  BEGIN
      IF query is null THEN
          return NULL;
      END IF;
      return jchem_blob_pkg.molconvertb(
          query, inputFormat, type, options, null, null, null, null, null, null, tmpBlob);
  END;

  FUNCTION transform(reaction VARCHAR2, reactants VARCHAR2)
            RETURN VARCHAR2 PARALLEL_ENABLE AS
  BEGIN
      IF reaction IS NULL OR reactants IS NULL THEN
          return NULL;
      ELSE
          return jchem_core_pkg.react(reaction, reactants, null, null, null,
                              'method:n mappingstyle:d permuteReactants:y', 't');
      END IF;
  END;

  FUNCTION react(reaction VARCHAR2, reactants VARCHAR2,
                           options VARCHAR2)
         RETURN VARCHAR2 PARALLEL_ENABLE AS
  BEGIN
    IF reaction IS NULL OR reactants IS NULL THEN
      return NULL;
    ELSE
        return jchem_core_pkg.react(reaction, reactants, null, null, null,
                                    options, 't');
    END IF;
  END;

  FUNCTION react4(reaction VARCHAR2, reactant1 VARCHAR2,
              reactant2 VARCHAR2, reactant3 VARCHAR2, reactant4 VARCHAR2,
              options VARCHAR2)
      RETURN VARCHAR2 PARALLEL_ENABLE AS
  BEGIN
      IF reaction IS NULL OR reactant1 IS NULL THEN
            return NULL;
      ELSE
          return jchem_core_pkg.react(reaction, reactant1, reactant2, reactant3,
                reactant4, options, null);
      END IF;
  END;

  FUNCTION react4(reaction CLOB, reactant1 CLOB, reactant2 CLOB,
                  reactant3 CLOB, reactant4 CLOB, options VARCHAR2,
                  tempClob CLOB)
      RETURN CLOB PARALLEL_ENABLE AS
  BEGIN
      IF reaction IS NULL OR reactant1 IS NULL THEN
            return NULL;
      ELSE
          return jchem_clob_pkg.react(reaction, reactant1, reactant2,
                          reactant3, reactant4, options, tempClob, null);
      END IF;
  END;

  FUNCTION react4(reaction BLOB, reactant1 BLOB, reactant2 BLOB, reactant3
              BLOB, reactant4 BLOB, options VARCHAR2, tempBlob BLOB)
      RETURN BLOB PARALLEL_ENABLE AS
  BEGIN
      IF reaction IS NULL OR reactant1 IS NULL THEN
            return NULL;
      ELSE
          return jchem_blob_pkg.react(reaction, reactant1, reactant2,
                          reactant3, reactant4, options, tempBlob, null);
      END IF;
  END;

  FUNCTION t_react4(reaction VARCHAR2, reactant1 VARCHAR2, reactant2
          VARCHAR2, reactant3 VARCHAR2, reactant4 VARCHAR2, options VARCHAR2)
  RETURN CHAR_PRODUCT_ARRAY PARALLEL_ENABLE PIPELINED AS
    cp_arr CHAR_PRODUCT_ARRAY;
  BEGIN
    IF reaction is null or reactant1 is null THEN
      return;
    END IF;

    cp_arr := jchem_core_pkg.react_arr(reaction, reactant1, reactant2,
                                        reactant3, reactant4, options, null);
    for i in 1..cp_arr.count loop
      pipe row(cp_arr(i));
    end loop;
    return;
  END;

  FUNCTION t_react4(reaction CLOB, reactant1 CLOB, reactant2 CLOB, reactant3
                     CLOB, reactant4 CLOB, options VARCHAR2)
  RETURN CLOB_PRODUCT_ARRAY PARALLEL_ENABLE PIPELINED AS
    cp_arr CLOB_PRODUCT_ARRAY;
  BEGIN
    IF reaction is null or reactant1 is null THEN 
      return; 
    END IF;

    cp_arr := jchem_clob_pkg.react_arr(reaction, reactant1, reactant2,
                reactant3, reactant4, options, null);

    for i in 1..cp_arr.count loop
      pipe row(cp_arr(i));
    end loop;

    return;
  END;

  FUNCTION t_react4(reaction BLOB, reactant1 BLOB, reactant2 BLOB, reactant3
                     BLOB, reactant4 BLOB, options VARCHAR2)
  RETURN BLOB_PRODUCT_ARRAY PARALLEL_ENABLE PIPELINED AS
    bp_arr BLOB_PRODUCT_ARRAY;
  BEGIN
    IF reaction is null or reactant1 is null THEN 
      return; 
    END IF;

    bp_arr := jchem_blob_pkg.react_arr(reaction, reactant1, reactant2,
                reactant3, reactant4, options, null);

    for i in 1..bp_arr.count loop
      pipe row(bp_arr(i));
    end loop;

    return;
  END;

  FUNCTION standardize(structure VARCHAR2, param VARCHAR2)
         RETURN VARCHAR2 PARALLEL_ENABLE AS
  BEGIN
    IF structure IS NULL THEN
        return NULL;
    ELSE
        return jchem_core_pkg.standardize(structure, param);
    END IF;
  END;

  FUNCTION standardize(structure CLOB, param VARCHAR2, tempClob CLOB)
      RETURN CLOB PARALLEL_ENABLE AS
  BEGIN
      IF structure IS NULL THEN
          return NULL;
      ELSE
          return jchem_clob_pkg.standardize(structure, param, tempClob);
      END IF;
  END;

  FUNCTION standardize(structure BLOB, param VARCHAR2, tempBlob BLOB)
      RETURN BLOB PARALLEL_ENABLE AS
  BEGIN
      IF structure IS NULL THEN
          return NULL;
      ELSE
          return jchem_blob_pkg.standardize(structure, param, tempBlob);
      END IF;
  END;

  FUNCTION formula_eq(mol VARCHAR2, query VARCHAR2) RETURN NUMBER PARALLEL_ENABLE AS
  BEGIN
      IF mol is null THEN return null; END IF;
      IF query is null THEN return null; END IF;
      IF formula(mol) = query THEN 
          RETURN 1;
      ELSE
          RETURN 0;
      END IF;
  END;

  FUNCTION formula_eq(mol clob, query VARCHAR2) RETURN NUMBER PARALLEL_ENABLE AS
  BEGIN
      IF mol is null THEN return null; END IF;
      IF query is null THEN return null; END IF;
      IF formula(mol) = query THEN 
          RETURN 1;
      ELSE
          RETURN 0;
      END IF;
  END;

  FUNCTION formula_eq(mol blob, query VARCHAR2) RETURN NUMBER PARALLEL_ENABLE AS
  BEGIN
      IF mol is null THEN return null; END IF;
      IF query is null THEN return null; END IF;
      IF formula(mol) = query THEN 
          RETURN 1;
      ELSE
          RETURN 0;
      END IF;
  END;

  FUNCTION insertable(structure VARCHAR2, schemaname VARCHAR2, tablename VARCHAR2, columnname VARCHAR2)
           RETURN NUMBER PARALLEL_ENABLE
   as language java name
    'chemaxon.jchem.cartridge.oraresident.JCartDml.isInsertable(java.lang.String, java.lang.String, java.lang.String, java.lang.String) return int';  

  FUNCTION insertable(structure CLOB, schemaname VARCHAR2, tablename VARCHAR2, columnname VARCHAR2)
           RETURN NUMBER PARALLEL_ENABLE
   as language java name
    'chemaxon.jchem.cartridge.oraresident.JCFunctionsClob.isInsertable(oracle.sql.CLOB, java.lang.String, java.lang.String, java.lang.String) return int';  

  FUNCTION insertable(structure BLOB, schemaname VARCHAR2, tablename VARCHAR2, columnname VARCHAR2)
           RETURN NUMBER PARALLEL_ENABLE
   as language java name
    'chemaxon.jchem.cartridge.oraresident.JCFunctionsBlob.isInsertable(oracle.sql.BLOB, java.lang.String, java.lang.String, java.lang.String) return int';  
	
  FUNCTION fuse(first_structure IN VARCHAR2, second_structure IN VARCHAR2, inputformat IN VARCHAR2)
      RETURN VARCHAR2 PARALLEL_ENABLE AS
  BEGIN
    IF first_structure is null THEN
	  return second_structure;
	END IF;
	IF second_structure is null THEN
	  return first_structure;
	END IF;
    return jchem_core_pkg.fuse(first_structure, second_structure, inputformat);
  END;
  
  FUNCTION fuse(first_structure IN CLOB, second_structure IN CLOB, inputformat IN VARCHAR2, tmpClob CLOB)
      RETURN CLOB PARALLEL_ENABLE AS
  BEGIN
    IF first_structure is null THEN
	  return second_structure;
	END IF;
	IF second_structure is null THEN
	  return first_structure;
	END IF;
    return jchem_clob_pkg.fuse(first_structure, second_structure, inputformat, tmpClob);
  END;

  FUNCTION fuse(first_structure IN BLOB, second_structure IN BLOB, inputformat IN VARCHAR2, tmpBlob BLOB)
      RETURN BLOB DETERMINISTIC PARALLEL_ENABLE AS
  BEGIN
    IF first_structure is null THEN
	  return second_structure;
	END IF;
	IF second_structure is null THEN
	  return first_structure;
	END IF;
    return jchem_blob_pkg.fuse(first_structure, second_structure, inputformat, tmpBlob);
  END;

  function t_get_gmem_util return gmem_util_array pipelined as
    gmemutil_arr gmem_util_array;
  begin
    gmemutil_arr := jchem_core_pkg.get_gmemutil_arr();
    for i in 1..gmemutil_arr.count loop
      pipe row(gmemutil_arr(i));
    end loop;
  end;

  function t_get_taskinfo return taskinfo_array pipelined as
    taskinfo_arr taskinfo_array;
  begin
    taskinfo_arr := jchem_core_pkg.get_taskinfo_arr();
    for i in 1..taskinfo_arr.count loop
      pipe row(taskinfo_arr(i));
    end loop;
  end;
  
  end jcf;
/
show errors;

-- Should we rather declare this procedure to have definer privileges?
create or replace procedure jc_set_default_property(prop_name VARCHAR2,
                                                prop_value VARCHAR2)
   authid current_user as language java name
    'chemaxon.jchem.cartridge.oraresident.JFunctions.setDefaultProperty(java.lang.String, java.lang.String)';
/
show errors;

create or replace package rmi as
    procedure clear_directory_cache as language java name
	    'chemaxon.jchem.cartridge.sharedorajcsrv.rmi.client.RmiDirectory.clearRemoteRefCache()';
    procedure clear_directory_cache(serverName varchar2) as language java name
	    'chemaxon.jchem.cartridge.sharedorajcsrv.rmi.client.RmiDirectory.clearRemoteRef(java.lang.String)';
    procedure clear_locale as language java name
        'chemaxon.jchem.cartridge.sharedorajcsrv.rmi.client.RmiDirectory.clearLocale()';
  --call this procedure after changing some properties in jcart.properties file to load new values  
    procedure clear_property_cache as language java name
	    'chemaxon.jchem.cartridge.sharedorajcsrv.AdminProperties.clearPropertyCache()';
end rmi;
/
show errors;

create or replace package jcc_admin as
  function kill_task(task_id number) return number as language java name
      'chemaxon.jchem.cartridge.oraresident.AdminSp.killTask(long) return int';
  function remove_zombies(task_id number) return number as language java name
      'chemaxon.jchem.cartridge.oraresident.AdminSp.killTask(long) return int';
  procedure change_node(oldnodeconfig varchar2, newnodeconfig varchar2) as language java name
      'chemaxon.jchem.cartridge.oraresident.AdminSp.changeNode(java.lang.String,java.lang.String)';
  procedure restart_new_nodes as language java name
      'chemaxon.jchem.cartridge.oraresident.AdminSp.restartNewNodes()';
	  
  --only used for proof of concept
  procedure changecache(idxtablename VARCHAR2,cdid NUMBER,cdstructure VARCHAR2,event NUMBER) as language java name  
      'chemaxon.jchem.cartridge.oraresident.AdminSp.changeCache(java.lang.String,long,java.lang.String,long)';
end jcc_admin;
/
show errors;


create or replace package scanctx_cache authid current_user as
  function get_scan_context(scanId number) return SCAN_CONTEXT as language java name
        'chemaxon.jchem.cartridge.oraresident.ScanContext.getScanContextPlSql(long)
          return oracle.sql.STRUCT';
  function t_get_scanctxs return SCAN_CONTEXT_ARRAY parallel_enable pipelined;
  procedure clear as language java name
      'chemaxon.jchem.cartridge.oraresident.ScanContext.clear()';

  function get_scanctxs return SCAN_CONTEXT_ARRAY
      as language java name
      'chemaxon.jchem.cartridge.oraresident.ScanContext.getScanContexts()
        return oracle.sql.ARRAY';
  function get_scanctxs_checked return SCAN_CONTEXT_ARRAY;
end scanctx_cache;
/
show errors;

create or replace package body scanctx_cache as
  function get_scanctxs_checked return SCAN_CONTEXT_ARRAY as
    scanctxarr SCAN_CONTEXT_ARRAY;
  begin
    scanctxarr := get_scanctxs();
    if (scanctxarr is null) then
      jchem_core_pkg.handle_java_error(jchem_core_pkg.lastError);
    end if;
    return scanctxarr;
  end;

  function t_get_scanctxs return SCAN_CONTEXT_ARRAY
      parallel_enable pipelined as
    scanctxarr SCAN_CONTEXT_ARRAY;
  begin
    scanctxarr := get_scanctxs_checked();
    for i in 1..scanctxarr.count
    loop
      pipe row(scanctxarr(i));
    end loop;  

    return;
  end;
end scanctx_cache;
/
show errors

create or replace package error_cache authid current_user as
  function get_error_count(cnt out number) return varchar2
      as language java name
      'chemaxon.jchem.cartridge.oraresident.ErrorCacheSp.getErrorCount(int[])
        return java.lang.String';
  function get_errors(scanId number, cnt number) return ERROR_RECORD_ARRAY 
      as language java name
      'chemaxon.jchem.cartridge.oraresident.ErrorCacheSp.getErrors(java.lang.Long, int)
        return oracle.sql.ARRAY';
  function get_errors_checked(scanId number, chunk_size number) return
        ERROR_RECORD_ARRAY parallel_enable;

  function get_error_count return number;
  function t_get_errors(scanId number := null, chunk_size number) return
        ERROR_RECORD_ARRAY parallel_enable pipelined;
  function t_get_errors(scanId number := null) return ERROR_RECORD_ARRAY
      parallel_enable pipelined;

  procedure clear as language java name
      'chemaxon.jchem.cartridge.oraresident.ErrorCacheSp.clearCache()';
  procedure set_errrec_exptime(exptime number)
      as language java name
      'chemaxon.jchem.cartridge.oraresident.ErrorCacheSp.setErrorRecordExpirationTime(long)';
  procedure set_emptysess_exptime(exptime number)
      as language java name
      'chemaxon.jchem.cartridge.oraresident.ErrorCacheSp.setEmptySessionErrorsExpirationTime(long)';
  procedure set_remove_rec_on_access(do_remove number) as language java name
      'chemaxon.jchem.cartridge.oraresident.ErrorCacheSp.setRemoveRecordsOnAccess(int)';
      
end error_cache;
/
show errors;

create or replace package body error_cache as
  function get_error_count return number as
    errm varchar2(32767);
    cnt number;
  begin
    errm := get_error_count(cnt);
    if errm is not null then
      jchem_core_pkg.handle_java_error(errm);
    end if;
    return cnt;
  end;

  function get_errors_checked(scanId number, chunk_size number) return
      ERROR_RECORD_ARRAY parallel_enable as
    errarr ERROR_RECORD_ARRAY;
  begin
    errarr := get_errors(scanId, chunk_size);
    if (errarr is null) then
      jchem_core_pkg.handle_java_error(jchem_core_pkg.lastError);
    end if;
    return errarr;
  end;

  function t_get_errors(scanId number := null, chunk_size number) return
      ERROR_RECORD_ARRAY parallel_enable pipelined as
    errarr ERROR_RECORD_ARRAY;
  begin
    errarr := get_errors_checked(scanId, chunk_size);
    for i in 1..errarr.count
    loop
      pipe row(errarr(i));
      errarr := get_errors_checked(scanId, chunk_size);
    end loop;  

    return;
  end;

  function t_get_errors(scanId number := null)  return ERROR_RECORD_ARRAY
      parallel_enable pipelined as
    errarr ERROR_RECORD_ARRAY;
    chunk_size number := 1000;
  begin
    errarr := get_errors_checked(scanId, chunk_size);
    while (errarr.count > 0)
    loop
      for i in 1..errarr.count
      loop
        pipe row(errarr(i));
      end loop;  
      errarr := get_errors_checked(scanId, chunk_size);
    end loop;  

    return;
  end;

end error_cache;
/
show errors;

declare
  cursor get_prop_name_length is
    select * from user_tab_cols
    where table_name = 'JC_IDX_PROPERTY' and column_name = 'PROP_NAME';
begin
  for uu_rec in get_prop_name_length
  loop
    if uu_rec.data_length < 200 then
      execute immediate 'alter table jc_idx_property modify (prop_name varchar2(4000))';
    end if;
  end loop;
end;
/
show errors;

exit;
