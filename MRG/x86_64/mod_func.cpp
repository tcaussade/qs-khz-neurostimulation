#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _AXNODE_reg(void);
extern void _internode_reg(void);
extern void _juxtaparanode_reg(void);
extern void _KV12_reg(void);
extern void _Nav6_reg(void);
extern void _newNav6_reg(void);
extern void _node_ikil_reg(void);
extern void _node_reg(void);
extern void _xtra2_reg(void);
extern void _xtra_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"nrnmechanism/AXNODE.mod\"");
    fprintf(stderr, " \"nrnmechanism/internode.mod\"");
    fprintf(stderr, " \"nrnmechanism/juxtaparanode.mod\"");
    fprintf(stderr, " \"nrnmechanism/KV12.mod\"");
    fprintf(stderr, " \"nrnmechanism/Nav6.mod\"");
    fprintf(stderr, " \"nrnmechanism/newNav6.mod\"");
    fprintf(stderr, " \"nrnmechanism/node_ikil.mod\"");
    fprintf(stderr, " \"nrnmechanism/node.mod\"");
    fprintf(stderr, " \"nrnmechanism/xtra2.mod\"");
    fprintf(stderr, " \"nrnmechanism/xtra.mod\"");
    fprintf(stderr, "\n");
  }
  _AXNODE_reg();
  _internode_reg();
  _juxtaparanode_reg();
  _KV12_reg();
  _Nav6_reg();
  _newNav6_reg();
  _node_ikil_reg();
  _node_reg();
  _xtra2_reg();
  _xtra_reg();
}

#if defined(__cplusplus)
}
#endif
