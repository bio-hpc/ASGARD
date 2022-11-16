whenever sqlerror exit 1

ASSOCIATE STATISTICS WITH INDEXTYPES jc_idxtype USING jc_idxstat_typ
/

ASSOCIATE STATISTICS WITH FUNCTIONS
Compare_Func
USING jc_idxstat_typ
/

ASSOCIATE STATISTICS WITH FUNCTIONS
Compare_FuncC
USING jc_idxstat_typ
/

ASSOCIATE STATISTICS WITH FUNCTIONS
Compare_FuncCV
USING jc_idxstat_typ
/

ASSOCIATE STATISTICS WITH FUNCTIONS
Compare_FuncB
USING jc_idxstat_typ
/

ASSOCIATE STATISTICS WITH FUNCTIONS
Compare_FuncBV
USING jc_idxstat_typ
/

ASSOCIATE STATISTICS WITH FUNCTIONS
Contains_Func_deprecated
USING jc_idxstat_typ
/

ASSOCIATE STATISTICS WITH FUNCTIONS
Contains_FuncC_deprecated
USING jc_idxstat_typ
/

ASSOCIATE STATISTICS WITH FUNCTIONS
Contains_FuncB_deprecated
USING jc_idxstat_typ
/

ASSOCIATE STATISTICS WITH FUNCTIONS
Tanimoto_Func
USING jc_idxstat_typ
/

ASSOCIATE STATISTICS WITH FUNCTIONS
Tanimoto_FuncC
USING jc_idxstat_typ
/

ASSOCIATE STATISTICS WITH FUNCTIONS
Tanimoto_FuncB
USING jc_idxstat_typ
/

insert into jc_idx_property
values('cost.factors.COMPARE.default', '4000;0.45;1000.0;3796875.0;0.0;1000.0', null)
/

insert into jc_idx_property
values('cost.factors.CONTAINS.default', '4000;0.45;1000.0;3796875.0;0.0;1000.0', null)
/

insert into jc_idx_property
values('cost.factors.TANIMOTO.default', '4000;0.45;1000.0;3796875.0;0.0;1000.0', null)
/

quit
