DISASSOCIATE STATISTICS FROM INDEXTYPES jc_idxtype force
/

DISASSOCIATE STATISTICS FROM FUNCTIONS Compare_Func force
/
DISASSOCIATE STATISTICS FROM FUNCTIONS Compare_FuncC force
/
DISASSOCIATE STATISTICS FROM FUNCTIONS Compare_FuncCV force
/

DISASSOCIATE STATISTICS FROM FUNCTIONS Compare_FuncB force
/
DISASSOCIATE STATISTICS FROM FUNCTIONS Compare_FuncBV force
/

DISASSOCIATE STATISTICS FROM FUNCTIONS Contains_Func_deprecated force
/
DISASSOCIATE STATISTICS FROM FUNCTIONS Contains_FuncC_deprecated force
/
DISASSOCIATE STATISTICS FROM FUNCTIONS Contains_FuncB_deprecated force
/
DISASSOCIATE STATISTICS FROM FUNCTIONS Tanimoto_Func force
/
DISASSOCIATE STATISTICS FROM FUNCTIONS Tanimoto_FuncC force
/
DISASSOCIATE STATISTICS FROM FUNCTIONS Tanimoto_FuncB force
/

delete from jc_idx_property where prop_name = 'cost.factors.COMPARE.default'
/

delete from jc_idx_property where prop_name = 'cost.factors.CONTAINS.default'
/

delete from jc_idx_property where prop_name = 'cost.factors.TANIMOTO.default'
/
quit
