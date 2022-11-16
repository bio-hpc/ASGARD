whenever sqlerror exit 1

create or replace package jchem_opti_pkg authid current_user as
  procedure raise_error(msg VARCHAR2);

  function collect_idxstats(idxSchema varchar2, idxName varchar2,
    idxPartition varchar2, tblSchema varchar2, tblName varchar2,
    colName varchar2, tblPartition varchar2)
    return number as language java name
    'chemaxon.jchem.cartridge.oraresident.costestim.JcOptCollectDeleteStats.collectIdxStats(
    java.lang.String, java.lang.String, java.lang.String,
    java.lang.String, java.lang.String, java.lang.String,
    java.lang.String) return int';

  function delete_idxstats(idxSchema varchar2, idxName varchar2,
    idxPartition varchar2)
    return number as language java name
    'chemaxon.jchem.cartridge.oraresident.costestim.JcOptCollectDeleteStats.deleteTableStats(
    java.lang.String, java.lang.String, java.lang.String) return int';

  function get_numrows(idxSchema varchar2, idxName varchar2,
            idxPartition varchar2) return number as language java name
    'chemaxon.jchem.cartridge.oraresident.costestim.JcOptCollectDeleteStats.getNumRows(
    java.lang.String, java.lang.String, java.lang.String) return int';
  
  function init_selectivity(
            pred    sys.ODCIPredInfo,
            args    sys.ODCIArgDescList,
            errid   in out number) return number;

  function get_selectivity(opSchema VARCHAR2, opName VARCHAR2,
        targetTableSchema VARCHAR2, targetTableName VARCHAR2,
        targetTableCol VARCHAR2, targetTablePartLower VARCHAR2,
        targetTablePartUpper VARCHAR2, rstart NUMBER, rstop NUMBER,
        opFlag number, query VARCHAR2, options VARCHAR2)
      return number as language java
    name
    'chemaxon.jchem.cartridge.oraresident.costestim.JcOptSelectivity.getSelectivity(
                            java.lang.String, java.lang.String,
                            java.lang.String, java.lang.String,
                            java.lang.String, java.lang.String,
                            java.lang.String, java.lang.Double, java.lang.Double, int,
                            java.lang.String, java.lang.String) return float';

  function get_selectivity(opSchema VARCHAR2, opName VARCHAR2,
          targetTableSchema VARCHAR2, targetTableName VARCHAR2,
          targetTableCol VARCHAR2, targetTablePartLower VARCHAR2,
          targetTablePartUpper VARCHAR2, rstart NUMBER, rstop NUMBER,
          opFlag number, query BLOB, options VARCHAR2)
      return number as language java name
    'chemaxon.jchem.cartridge.oraresident.costestim.JcOptSelectivity.getSelectivity(
                            java.lang.String, java.lang.String,
                            java.lang.String, java.lang.String,
                            java.lang.String, java.lang.String,
                            java.lang.String, java.lang.Double, java.lang.Double, int,
                            oracle.sql.BLOB, java.lang.String) return float';

  function func_cost(func sys.ODCIFuncInfo, 
                     cost OUT sys.ODCICost,
                     options varchar2) return number;

  function get_cost_factor(costType varchar2, opName varchar2, options varchar2,
                           resrc varchar2, idxSchema varchar2, idxName varchar2)
      return number as language java name
    'chemaxon.jchem.cartridge.oraresident.costestim.JcOptCosts.getCostFactor(
    java.lang.String, java.lang.String, java.lang.String, java.lang.String,
    java.lang.String, java.lang.String) return double';

  procedure set_volatile_cost_factors(opName varchar2,
                                    idxCpuCost number, idxIoCost number,
                                    idxNetworkCost number, funcCpuCost number,
                                    funcIoCost number, funcNetworkCost number)
          as language java name
    'chemaxon.jchem.cartridge.oraresident.costestim.JcOptCosts.setVolatileCostEstimation(
      java.lang.String, double, double, double, double, double, double)';


  procedure store_default_cost_factors as language java name
    'chemaxon.jchem.cartridge.oraresident.costestim.JcOptCosts.storeDefaultCostFactors()';

end jchem_opti_pkg;
/
show errors;

create or replace package body jchem_opti_pkg as
  function init_selectivity(
            pred    sys.ODCIPredInfo,
            args    sys.ODCIArgDescList,
            errid   in out number) return number as
  begin
    -- pred.MethodName: do not support package functions
    IF pred.MethodName IS NOT NULL THEN
      jcart_logger.error('ODCIStatsSelectivity(...)',
                         'No optimizer support for package functions...');
      jchem_opti_pkg.raise_error('No optimizer support for package functions');
      return ODCIConst.Error;
    END IF;

/*
    -- start value
    IF (args(1).ArgType != ODCIConst.ArgLit AND
            args(1).ArgType != ODCIConst.ArgNull) THEN
        jcart_logger.debug('ODCIStatsSelectivity(...)',
                'rstart=' || to_char(rstart)
                || ', argType=' || to_char(args(1).ArgType)
                || ', errid=' || to_char(errid)
                || ' (ArgOther=' || ODCIConst.ArgOther);
        RETURN ODCIConst.Error;
    ELSE
        errid := errid + 1;
    END IF;

    -- stop value
    IF (args(2).ArgType != ODCIConst.ArgLit AND
            args(2).ArgType != ODCIConst.ArgNull) THEN
        jcart_logger.debug('ODCIStatsSelectivity(...)', 'errid=' || errid);
        RETURN ODCIConst.Error;
    ELSE
        errid := errid + 1;
    END IF;

*/
    -- first argument of function
    IF (args(3).ArgType != ODCIConst.ArgCol) THEN
        jcart_logger.debug('ODCIStatsSelectivity(...)', 'errid=' || errid);
        RETURN ODCIConst.Error;
    ELSE
        jcart_logger.debug('ODCIStatsSelectivity(...)', 'First argument OK');
        errid := errid + 1;
    END IF;

/*
    -- second argument of function
    IF (args(4).ArgType != ODCIConst.ArgLit AND
            args(4).ArgType != ODCIConst.ArgNull) THEN
        jcart_logger.debug('ODCIStatsSelectivity(...)', 'errid=' || errid);
        RETURN ODCIConst.Error;
    ELSE
        errid := errid + 1;
    END IF;

    -- third argument of function
    IF (args.count >= 5) THEN
        IF (args(5).ArgType != ODCIConst.ArgLit AND
                args(5).ArgType != ODCIConst.ArgNull) THEN
            jcart_logger.debug('ODCIStatsSelectivity(...)', 'errid=' || errid);
            RETURN ODCIConst.Error;
        END IF;
    ELSE
        errid := errid + 1;
    END IF;
*/

    return ODCIConst.Success;
  end;

  function func_cost(func sys.ODCIFuncInfo, 
                     cost OUT sys.ODCICost,
                     options varchar2) return number as
    cpuCostFactor number(16, 2);
    ioCostFactor  number(16, 2);
    netCostFactor number(16, 2);
  BEGIN
    -- Uses estimation calculated on-the-fly and stored
    -- in "package variables". We may need to move on-the-fly
    -- estimation into statistics gathering later.

    jcart_logger.debug('ODCIStatsFunctionCost(...)', 'called');
    

    cpuCostFactor := jchem_opti_pkg.get_cost_factor('func', func.ObjectName,
                        options, 'cpu', null, null);
    if cpuCostFactor = -1 then
        cpuCostFactor := null;
    end if;

    ioCostFactor := jchem_opti_pkg.get_cost_factor('func', func.ObjectName,
                        options, 'io', null, null);
    if ioCostFactor = -1 then
        ioCostFactor := null;
    end if;
    netCostFactor := jchem_opti_pkg.get_cost_factor('func', func.ObjectName,
                        options, 'net', null, null);
    if netCostFactor = -1 then
        netCostFactor := null;
    end if;

    cost := sys.ODCICost(NULL, NULL, NULL, NULL);
    cost.CPUCost := ceil(cpuCostFactor);
    cost.IOCost := ceil(ioCostFactor);
    -- cost.NetworkCost := ceil(netCostFactor);

    jcart_logger.debug('ODCIStatsFunctionCost(...)',
        'cpu=' || to_char(cost.CPUCost)
        || ', io=' || to_char(cost.IOCost)
        || ', net=' || to_char(cost.NetworkCost)
        || '. Returning success...');
    RETURN ODCIConst.Success;
  END;

  procedure raise_error(msg VARCHAR2) IS
  begin
        raise_application_error(-20101, msg);
  end;
end jchem_opti_pkg;
/
show errors;

CREATE OR REPLACE TYPE jc_idxstat_typ AUTHID CURRENT_USER AS OBJECT
(
  dummy number,

  STATIC FUNCTION ODCIGetInterfaces(
                  ifclist OUT sys.ODCIObjectList) RETURN NUMBER,

  STATIC FUNCTION ODCIStatsCollect(
                  col        sys.ODCIColInfo,
                  options    sys.ODCIStatsOptions,
                  statistics OUT RAW,
                  env        sys.ODCIEnv) RETURN NUMBER,

  STATIC FUNCTION ODCIStatsCollect(
                  ia         sys.ODCIIndexInfo,
                  options    sys.ODCIStatsOptions,
                  statistics OUT RAW,
                  env        sys.ODCIEnv) return NUMBER,

  STATIC FUNCTION ODCIStatsDelete(
                  col        sys.ODCIColInfo,
                  statistics OUT RAW,
                  env        sys.ODCIEnv) return NUMBER,

  STATIC FUNCTION ODCIStatsDelete(
                  col        sys.ODCIIndexInfo,
                  statistics OUT RAW,
                  env        sys.ODCIEnv) return NUMBER,

  STATIC FUNCTION ODCIStatsSelectivity(
                  pred    sys.ODCIPredInfo,
                  sel     OUT NUMBER,
                  args    sys.ODCIArgDescList,
                  rstart  NUMBER,
                  rstop   NUMBER,
                  target  VARCHAR2,
                  query   VARCHAR2,
                  env     sys.ODCIEnv)
                  return NUMBER,

  STATIC FUNCTION ODCIStatsSelectivity(
                  pred    sys.ODCIPredInfo,
                  sel     OUT NUMBER,
                  args    sys.ODCIArgDescList,
                  rstart  NUMBER,
                  rstop   NUMBER,
                  target  BLOB,
                  query   BLOB,
                  env     sys.ODCIEnv)
                  return NUMBER,

  STATIC FUNCTION ODCIStatsSelectivity(
                  pred    sys.ODCIPredInfo,
                  sel     OUT NUMBER,
                  args    sys.ODCIArgDescList,
                  rstart  NUMBER,
                  rstop   NUMBER,
                  target  VARCHAR2,
                  query   VARCHAR2,
                  options VARCHAR2,
                  env     sys.ODCIEnv)
                  return NUMBER,

  STATIC FUNCTION ODCIStatsSelectivity(
                  pred    sys.ODCIPredInfo,
                  sel     OUT NUMBER,
                  args    sys.ODCIArgDescList,
                  rstart  NUMBER,
                  rstop   NUMBER,
                  target  BLOB,
                  query   BLOB,
                  options VARCHAR2,
                  env     sys.ODCIEnv)
                  return NUMBER,

  STATIC FUNCTION ODCIStatsSelectivity(
                  pred    sys.ODCIPredInfo,
                  sel     OUT NUMBER,
                  args    sys.ODCIArgDescList,
                  rstart  NUMBER,
                  rstop   NUMBER,
                  target  BLOB,
                  query   VARCHAR2,
                  options VARCHAR2,
                  env     sys.ODCIEnv)
                  return NUMBER,

  STATIC FUNCTION ODCIStatsIndexCost(
                  ia      sys.ODCIIndexInfo,
                  sel     NUMBER,
                  cost    OUT sys.ODCICost,
                  qi      sys.ODCIQueryInfo,
                  pred    sys.ODCIPredInfo,
                  args    sys.ODCIArgDescList,
                  rstart  NUMBER,
                  rstop   NUMBER,
                  query   VARCHAR2,
                  env sys.ODCIEnv) return NUMBER,

  STATIC FUNCTION ODCIStatsIndexCost(
                  ia      sys.ODCIIndexInfo,
                  sel     NUMBER,
                  cost    OUT sys.ODCICost,
                  qi      sys.ODCIQueryInfo,
                  pred    sys.ODCIPredInfo,
                  args    sys.ODCIArgDescList,
                  rstart  NUMBER,
                  rstop   NUMBER,
                  query   BLOB,
                  env sys.ODCIEnv) return NUMBER,

  STATIC FUNCTION ODCIStatsIndexCost(
                  ia      sys.ODCIIndexInfo,
                  sel     NUMBER,
                  cost    OUT sys.ODCICost,
                  qi      sys.ODCIQueryInfo,
                  pred    sys.ODCIPredInfo,
                  args    sys.ODCIArgDescList,
                  rstart  NUMBER,
                  rstop   NUMBER,
                  query   VARCHAR2,
                  options VARCHAR2,
                  env sys.ODCIEnv) return NUMBER,

  STATIC FUNCTION ODCIStatsIndexCost(
                  ia      sys.ODCIIndexInfo,
                  sel     NUMBER,
                  cost    OUT sys.ODCICost,
                  qi      sys.ODCIQueryInfo,
                  pred    sys.ODCIPredInfo,
                  args    sys.ODCIArgDescList,
                  rstart  NUMBER,
                  rstop   NUMBER,
                  query   BLOB,
                  options VARCHAR2,
                  env sys.ODCIEnv) return NUMBER,

  STATIC FUNCTION ODCIStatsFunctionCost(
                  func    sys.ODCIFuncInfo,
                  cost    OUT sys.ODCICost,
                  args    sys.ODCIArgDescList,
                  target  VARCHAR2,
                  query   VARCHAR2,
                  env     sys.ODCIEnv) return NUMBER,

  STATIC FUNCTION ODCIStatsFunctionCost(
                  func    sys.ODCIFuncInfo,
                  cost    OUT sys.ODCICost,
                  args    sys.ODCIArgDescList,
                  target  BLOB,
                  query   BLOB,
                  env     sys.ODCIEnv) return NUMBER,

  STATIC FUNCTION ODCIStatsFunctionCost(
                  func    sys.ODCIFuncInfo,
                  cost    OUT sys.ODCICost,
                  args    sys.ODCIArgDescList,
                  target  VARCHAR2,
                  query   VARCHAR2,
                  options VARCHAR2,
                  env     sys.ODCIEnv) return NUMBER,

  STATIC FUNCTION ODCIStatsFunctionCost(
                  func    sys.ODCIFuncInfo,
                  cost    OUT sys.ODCICost,
                  args    sys.ODCIArgDescList,
                  target  BLOB,
                  query   VARCHAR2,
                  options VARCHAR2,
                  env     sys.ODCIEnv) return NUMBER,

  STATIC FUNCTION ODCIStatsFunctionCost(
                  func    sys.ODCIFuncInfo,
                  cost    OUT sys.ODCICost,
                  args    sys.ODCIArgDescList,
                  target  BLOB,
                  query   BLOB,
                  options VARCHAR2,
                  env     sys.ODCIEnv) return NUMBER
);
/
show errors


create or replace
type body jc_idxstat_typ is

  STATIC FUNCTION ODCIGetInterfaces(ifclist OUT sys.ODCIObjectList)
							RETURN NUMBER IS
  BEGIN
    ifclist := sys.ODCIObjectList(sys.ODCIObject('SYS','ODCISTATS2'));
    RETURN ODCIConst.Success;
  END ODCIGetInterfaces;

  STATIC FUNCTION ODCIStatsCollect(
                  col        sys.ODCIColInfo,
                  options    sys.ODCIStatsOptions,
                  statistics OUT RAW,
                  env        sys.ODCIEnv) RETURN NUMBER IS
  BEGIN
    -- TODO:
    -- Check to see if it worthwhile doing something meaningful
    -- here.
    jcart_logger.debug('ODCIStatsCollect(col...)', 'called');
    statistics := NULL;
  END;

  STATIC FUNCTION ODCIStatsCollect(
                  ia         sys.ODCIIndexInfo,
                  options    sys.ODCIStatsOptions,
                  statistics OUT RAW,
                  env        sys.ODCIEnv) return NUMBER IS
    query_str varchar2(4000);
  BEGIN
    jcart_logger.debug('ODCIStatsCollect(ia...)', 'called');
    statistics := NULL;

    -- Currently:
    --   it only collects stats on the index table;
    --
    -- Would be nice:
    --   Calculate time based on a sample of structures
    --   using a pre-configured set of query structures.
    --   Use fingerprint pre-screening only without graph search
    --   Record average cost per hit
    --   TEST: Trace variation of costs per hit
    --   TEST: Use complete search
    --   TEST: Trace variation of costs of complete searches
    return jchem_opti_pkg.collect_idxstats(ia.IndexSchema, ia.IndexName,
        ia.IndexPartition,
        ia.IndexCols(1).TableSchema, ia.IndexCols(1).TableName,
        ia.IndexCols(1).ColName, ia.IndexCols(1).TablePartition);
  END;

  STATIC FUNCTION ODCIStatsDelete(
                  col        sys.ODCIColInfo,
                  statistics OUT RAW,
                  env        sys.ODCIEnv) return NUMBER IS
  BEGIN
    jcart_logger.debug('ODCIStatsDelete(col...)', 'called');
    statistics := null;
    RETURN 0;
  END;

  STATIC FUNCTION ODCIStatsDelete(
                  col        sys.ODCIIndexInfo,
                  statistics OUT RAW,
                  env        sys.ODCIEnv) return NUMBER IS
  BEGIN
    jcart_logger.debug('ODCIStatsDelete(ia...)', 'called');
--    jchem_opti_pkg.delete_idxstats(ia.IndexSchema, ia.IndexName,
--                    ia.IndexPartition);
    statistics := null;
    RETURN 0;
  END;
  
    STATIC FUNCTION ODCIStatsSelectivity(
                  pred    sys.ODCIPredInfo,
                  sel     OUT NUMBER,
                  args    sys.ODCIArgDescList,
                  rstart  NUMBER,
                  rstop   NUMBER,
                  target  VARCHAR2,
                  query   VARCHAR2,
                  env     sys.ODCIEnv)
                  return NUMBER IS
  BEGIN
    return ODCIStatsSelectivity(pred, sel, args, rstart, rstop, target, query,
                                null, env);
  END;

  STATIC FUNCTION ODCIStatsSelectivity(
                  pred    sys.ODCIPredInfo,
                  sel     OUT NUMBER,
                  args    sys.ODCIArgDescList,
                  rstart  NUMBER,
                  rstop   NUMBER,
                  target  BLOB,
                  query   BLOB,
                  env     sys.ODCIEnv)
                  return NUMBER IS
  BEGIN
    return ODCIStatsSelectivity(pred, sel, args, rstart, rstop, target, query,
                                null, env);
  END;

  STATIC FUNCTION ODCIStatsSelectivity(
                  pred    sys.ODCIPredInfo,
                  sel     OUT NUMBER,
                  args    sys.ODCIArgDescList,
                  rstart  NUMBER,
                  rstop   NUMBER,
                  target  VARCHAR2,
                  query   VARCHAR2,
                  options VARCHAR2,
                  env     sys.ODCIEnv)
                  return NUMBER IS
    errid integer := 0;
    retval integer;
  BEGIN
    jcart_logger.debug('ODCIStatsSelectivity(..., VARCHAR2, VARCHAR2, VARCHAR2)', 'started');

    retval := jchem_opti_pkg.init_selectivity(pred, args, errid);
    if retval <> ODCIConst.Success then
      return retval;
    end if;

    jcart_logger.debug('ODCIStatsSelectivity(...)',
                       'Calling jchem_opti_pkg.get_selectivity...');
    sel := jchem_opti_pkg.get_selectivity(pred.ObjectSchema, pred.ObjectName,
                    args(3).TableSchema, args(3).TableName, args(3).ColName,
                    args(3).TablePartitionLower, args(3).TablePartitionUpper,
                    rstart, rstop, pred.Flags, query, options);

        jcart_logger.debug('ODCIStatsSelectivity(...)',
            'returning=' || to_char(sel) || '...');
    RETURN ODCIConst.Success;
  END;

  STATIC FUNCTION ODCIStatsSelectivity(
                  pred    sys.ODCIPredInfo,
                  sel     OUT NUMBER,
                  args    sys.ODCIArgDescList,
                  rstart  NUMBER,
                  rstop   NUMBER,
                  target  BLOB,
                  query   BLOB,
                  options VARCHAR2,
                  env     sys.ODCIEnv)
                  return NUMBER IS
    errid integer := 0;
    retval integer;
  BEGIN
    jcart_logger.debug('ODCIStatsSelectivity(..., VARCHAR2, VARCHAR2, VARCHAR2)',
                        'called...');
    
    retval := jchem_opti_pkg.init_selectivity(pred, args, errid);
    if retval <> ODCIConst.Success then
      return retval;
    end if;

    jcart_logger.debug('ODCIStatsSelectivity(...)',
                       'Calling jchem_opti_pkg.get_selectivity...');
    sel := jchem_opti_pkg.get_selectivity(pred.ObjectSchema, pred.ObjectName,
                    args(3).TableSchema, args(3).TableName, args(3).ColName,
                    args(3).TablePartitionLower, args(3).TablePartitionUpper,
                    rstart, rstop, pred.Flags, query, options);

    jcart_logger.debug('ODCIStatsSelectivity(...)',
            'returning=' || to_char(sel) || '...');
    RETURN ODCIConst.Success;
  END;

  STATIC FUNCTION ODCIStatsSelectivity(
                  pred    sys.ODCIPredInfo,
                  sel     OUT NUMBER,
                  args    sys.ODCIArgDescList,
                  rstart  NUMBER,
                  rstop   NUMBER,
                  target  BLOB,
                  query   VARCHAR2,
                  options VARCHAR2,
                  env     sys.ODCIEnv)
                  return NUMBER IS
    errid integer := 0;
    retval integer;
  BEGIN
    jcart_logger.debug('ODCIStatsSelectivity(..., VARCHAR2, VARCHAR2, VARCHAR2)',
                        'called...');

    retval := jchem_opti_pkg.init_selectivity(pred, args, errid);
    if retval <> ODCIConst.Success then
      return retval;
    end if;

    jcart_logger.debug('ODCIStatsSelectivity(...)',
                       'Calling jchem_opti_pkg.get_selectivity...');
    sel := jchem_opti_pkg.get_selectivity(pred.ObjectSchema, pred.ObjectName,
                    args(3).TableSchema, args(3).TableName, args(3).ColName,
                    args(3).TablePartitionLower, args(3).TablePartitionUpper,
                    rstart, rstop, pred.Flags, query, options);

    jcart_logger.debug('ODCIStatsSelectivity(...)',
            'returning=' || to_char(sel) || '...');
    RETURN ODCIConst.Success;
  END;

  STATIC FUNCTION ODCIStatsIndexCost(
                  ia      sys.ODCIIndexInfo,
                  sel     NUMBER,
                  cost    OUT sys.ODCICost,
                  qi      sys.ODCIQueryInfo,
                  pred    sys.ODCIPredInfo,
                  args    sys.ODCIArgDescList,
                  rstart  NUMBER,
                  rstop   NUMBER,
                  query   VARCHAR2,
                  env sys.ODCIEnv) return NUMBER IS
  BEGIN
    return ODCIStatsIndexCost(ia, sel, cost, qi, pred, args, rstart, rstop,
                              query, null, env);
  END;

  STATIC FUNCTION ODCIStatsIndexCost(
                  ia      sys.ODCIIndexInfo,
                  sel     NUMBER,
                  cost    OUT sys.ODCICost,
                  qi      sys.ODCIQueryInfo,
                  pred    sys.ODCIPredInfo,
                  args    sys.ODCIArgDescList,
                  rstart  NUMBER,
                  rstop   NUMBER,
                  query   BLOB,
                  env sys.ODCIEnv) return NUMBER IS
  BEGIN
    return ODCIStatsIndexCost(ia, sel, cost, qi, pred, args, rstart, rstop,
                              query, null, env);
  END;

  STATIC FUNCTION ODCIStatsIndexCost(
                  ia      sys.ODCIIndexInfo,
                  sel     NUMBER,
                  cost    OUT sys.ODCICost,
                  qi      sys.ODCIQueryInfo,
                  pred    sys.ODCIPredInfo,
                  args    sys.ODCIArgDescList,
                  rstart  NUMBER,
                  rstop   NUMBER,
                  query   VARCHAR2,
                  options VARCHAR2,
                  env sys.ODCIEnv) return NUMBER IS
    ixowner VARCHAR2(70);
    ixtable VARCHAR2(70);
    dotpos integer;
    numblocks INTEGER;
    numrows INTEGER;
    CURSOR c1(ownr VARCHAR2, tab VARCHAR2) IS
        SELECT * FROM all_tables WHERE owner = ownr and table_name = tab;

    cpuCostFactor number(16, 2);
    ioCostFactor  number(16, 2);
    netCostFactor number(16, 2);
  BEGIN
    jcart_logger.debug('ODCIStatsIndexCost(...)', 'called (really)');

    IF sel IS NULL THEN
        jcart_logger.warning('ODCIStatsIndexCost(...)', 'sel is null!');
        RETURN ODCIConst.Error;
    END IF;

    jcart_logger.debug('ODCIStatsIndexCost(...)', 'Finding out index table name...');
    -- Get name of table implementing the domain index
    ixtable := jchem_core_pkg.get_idxtable_qname(ia.IndexSchema,
                                ia.IndexName, ia.IndexPartition, 0);
    jcart_logger.debug('ODCIStatsIndexCost(...)', 'ixtable=' || ixtable);

    if ixtable is null then -- JChem table
      ixowner := ia.IndexCols(1).TableSchema;
      ixtable := ia.IndexCols(1).TableName;
    else
      dotpos := instr(ixtable, '.');
      ixowner := substr(ixtable, 1, dotpos - 1);
      ixtable := substr(ixtable, dotpos+1, length(ixtable));
    end if;

    jcart_logger.debug('ODCIStatsIndexCost(...)', 'ixowner=' || ixowner);
    FOR get_table IN c1(ixowner, upper(ixtable)) LOOP
        numblocks := get_table.blocks;
        numrows := get_table.num_rows;
    END LOOP;

    IF numblocks IS NULL THEN
        jcart_logger.warning('ODCIStatsIndexCost(...)', 'no statistics for '
                || ixowner || '.' || ixtable);
        -- Exit if there are no user-defined statistics for the index
        RETURN ODCIConst.Error;
    END IF;

    jcart_logger.debug('ODCIStatsIndexCost(...)', 'numblocks='
        || to_char(numblocks) || ', numrows=' || to_char(numrows));

    cpuCostFactor := jchem_opti_pkg.get_cost_factor('index', pred.ObjectName,
                        options, 'cpu', ia.IndexSchema, ia.IndexName);
    jcart_logger.debug('ODCIStatsIndexCost(...)', 'cpuCostFactor=' || to_char(cpuCostFactor));
    if cpuCostFactor = -1 then
        cpuCostFactor := null;
    end if;

    ioCostFactor := jchem_opti_pkg.get_cost_factor('index', pred.ObjectName,
                        options, 'io', ia.IndexSchema, ia.IndexName);
    jcart_logger.debug('ODCIStatsIndexCost(...)', 'ioCostFactor=' || to_char(ioCostFactor));
    if ioCostFactor = -1 then
        ioCostFactor := null;
    end if;

    netCostFactor := jchem_opti_pkg.get_cost_factor('index', pred.ObjectName,
                        options, 'net', ia.IndexSchema, ia.IndexName);
    if netCostFactor = -1 then
        netCostFactor := null;
    end if;

    jcart_logger.debug('ODCIStatsIndexCost(...)',
                         'netCostFactor = ' || to_char(netCostFactor)
                         || ' ioCostFactor = ' || to_char(ioCostFactor)
                         || ' cpuCostFactor = ' || cpuCostFactor);

    if netCostFactor is null and ioCostFactor is null and cpuCostFactor is null then
        jcart_logger.warning('ODCIStatsIndexCost(...)',
                             'Everything is null, returning error...');
        RETURN ODCIConst.Error;
    end if;

    cost := sys.ODCICost(NULL, NULL, NULL, NULL);
    cost.CPUCost := ceil(cpuCostFactor * (sel/100)*numblocks);
    cost.IOCost := ceil(ioCostFactor * (sel/100)*numblocks);
    cost.NetworkCost := ceil(netCostFactor * (sel/100)*numblocks);

    jcart_logger.debug('ODCIStatsIndexCost(...)', 'Returning Success...');
    RETURN ODCIConst.Success;
  END;

  STATIC FUNCTION ODCIStatsIndexCost(
                  ia      sys.ODCIIndexInfo,
                  sel     NUMBER,
                  cost    OUT sys.ODCICost,
                  qi      sys.ODCIQueryInfo,
                  pred    sys.ODCIPredInfo,
                  args    sys.ODCIArgDescList,
                  rstart  NUMBER,
                  rstop   NUMBER,
                  query   BLOB,
                  options VARCHAR2,
                  env sys.ODCIEnv) return NUMBER IS
  BEGIN
    return ODCIStatsIndexCost(ia, sel, cost, qi, pred, args, rstart, rstop,
                              query, options, env);
  END;

  STATIC FUNCTION ODCIStatsFunctionCost(
                  func    sys.ODCIFuncInfo,
                  cost    OUT sys.ODCICost,
                  args    sys.ODCIArgDescList,
                  target  VARCHAR2,
                  query   VARCHAR2,
                  env sys.ODCIEnv) return NUMBER IS
  BEGIN
    return ODCIStatsFunctionCost(func, cost, args, target, query, null, env);
  END;

  STATIC FUNCTION ODCIStatsFunctionCost(
                  func    sys.ODCIFuncInfo,
                  cost    OUT sys.ODCICost,
                  args    sys.ODCIArgDescList,
                  target  BLOB,
                  query   BLOB,
                  env sys.ODCIEnv) return NUMBER IS
  BEGIN
    return ODCIStatsFunctionCost(func, cost, args, target, query, null, env);
  END;

  STATIC FUNCTION ODCIStatsFunctionCost(
                  func    sys.ODCIFuncInfo,
                  cost    OUT sys.ODCICost,
                  args    sys.ODCIArgDescList,
                  target  VARCHAR2,
                  query   VARCHAR2,
                  options VARCHAR2,
                  env sys.ODCIEnv) return NUMBER IS
  begin
    return jchem_opti_pkg.func_cost(func, cost, options);
  end;

  STATIC FUNCTION ODCIStatsFunctionCost(
                  func    sys.ODCIFuncInfo,
                  cost    OUT sys.ODCICost,
                  args    sys.ODCIArgDescList,
                  target  BLOB,
                  query   BLOB,
                  options VARCHAR2,
                  env sys.ODCIEnv) return NUMBER IS
  BEGIN
    return  jchem_opti_pkg.func_cost(func, cost, options);
  END;

  STATIC FUNCTION ODCIStatsFunctionCost(
                  func    sys.ODCIFuncInfo,
                  cost    OUT sys.ODCICost,
                  args    sys.ODCIArgDescList,
                  target  BLOB,
                  query   VARCHAR2,
                  options VARCHAR2,
                  env sys.ODCIEnv) return NUMBER IS
  BEGIN
    return jchem_opti_pkg.func_cost(func, cost, options);
  END;

END;
/
show errors;

quit
