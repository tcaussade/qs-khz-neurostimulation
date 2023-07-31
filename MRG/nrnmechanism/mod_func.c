#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _AXNODE_reg();
extern void _KV12_reg();
extern void _Nav6_reg();
extern void _internode_reg();
extern void _juxtaparanode_reg();
extern void _newNav6_reg();
extern void _node_reg();
extern void _node_ikil_reg();
extern void _xtra_reg();
extern void _xtra2_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," AXNODE.mod");
fprintf(stderr," KV12.mod");
fprintf(stderr," Nav6.mod");
fprintf(stderr," internode.mod");
fprintf(stderr," juxtaparanode.mod");
fprintf(stderr," newNav6.mod");
fprintf(stderr," node.mod");
fprintf(stderr," node_ikil.mod");
fprintf(stderr," xtra.mod");
fprintf(stderr," xtra2.mod");
fprintf(stderr, "\n");
    }
_AXNODE_reg();
_KV12_reg();
_Nav6_reg();
_internode_reg();
_juxtaparanode_reg();
_newNav6_reg();
_node_reg();
_node_ikil_reg();
_xtra_reg();
_xtra2_reg();
}
