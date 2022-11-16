#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "../../cifparse/cifparse.h"

int main( int argc, char* argv[] ) {

	int iblock, stat, i, iCat;
    int col1, col2, col3, col4, col5, col6, col7, col8, col9;
	char blockId[50] = "";
    char potential_type[50];
    strncpy( blockId, argv[1], 50 );
    
    /*
     **   Use the cifparse framework to extract information from an cif file:
     */
    ndb_cif_init();

    stat = ndb_cif_read_file(stdin);

    /*  Check for target data block.  */
    iblock = -1;
    for (i = 0; i < cifFiles.numDatablock; i++) {
        if (!strcmp(blockId, cifFiles.datablocks[i].datablockName)) {
            iblock = i;
            break;
        }
    }
    if (iblock == -1) {
        fprintf(stderr, "Target data block %s is not in file\n", blockId);
        exit(1);
    }

    iCat = get_category_index(iblock, "distance_restraint_list");
    if (iCat == -1) {
        fprintf(stderr, "Target category %s is not in data block %d\n",
                "distance_restraint_list", iblock);
        exit(1);
    }

    col1 = get_column_index(iblock, iCat, "potential_type"); 
    assert( col1 >= 0 );
    strncpy( potential_type, 
      cifFiles.datablocks[iblock].categories[iCat].rows[0].columns[col1], 50);
    fprintf( stderr, "Potential type is %s\n", potential_type );
   
    iCat = get_category_index(iblock, "distance_restraint");
    if (iCat == -1) {
        fprintf(stderr, "Target category %s is not in data block %d\n",
                "distance_restraint", iblock);
        exit(1);
    }
    fprintf(stderr, "Read %d rows in category %d\n",
            cifFiles.datablocks[iblock].categories[iCat].numRow, iCat);

    col1 = get_column_index(iblock, iCat, "restraint_id"); /* id */
    assert( col1 >= 0 );
    col2 = get_column_index(iblock, iCat, "sequence_code_1");
    assert( col2 >= 0 );
    col3 = get_column_index(iblock, iCat, "residue_type_1");
    assert( col3 >= 0 );
    col4 = get_column_index(iblock, iCat, "atom_name_1");
    assert( col4 >= 0 );
    col6 = get_column_index(iblock, iCat, "sequence_code_2");
    assert( col6 >= 0 );
    col7 = get_column_index(iblock, iCat, "residue_type_2");
    assert( col7 >= 0 );
    col8 = get_column_index(iblock, iCat, "atom_name_2");
    assert( col8 >= 0 );
    col9 = get_column_index(iblock, iCat, "upper_limit");
    assert( col9 >= 0 );

    /*  --- need some better error processing here!  */

    for (i = 0; i < cifFiles.datablocks[iblock].categories[iCat].numRow;
         i++) {

		printf( "%d    %d %s %s     %d %s %s   %10.3f\n", 
    atoi(cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col1]),
    atoi(cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col2]),
         cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col3],
         cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col4],
    atoi(cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col6]),
         cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col7],
         cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col8],
    atof(cifFiles.datablocks[iblock].categories[iCat].rows[i].columns[col9])
		);

    }

	/* probably need some calls here to free up memory cifparse has used */
}
