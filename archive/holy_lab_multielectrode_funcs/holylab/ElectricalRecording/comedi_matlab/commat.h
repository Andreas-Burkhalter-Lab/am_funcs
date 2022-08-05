struct commat_chanspec_struct {
  unsigned int *chanlist;
  unsigned int *rangelist;
  unsigned int *areflist;
  unsigned int *packedchanlist;
  int nchans_chan;
  int nchans_range;
  int nchans_aref;
  int nchans_packed;
};

typedef struct commat_chanspec_struct commat_chanspec;

extern bool IsScalar(const mxArray *);
extern int matstr_len(const mxArray *);
extern char* matstr_to_cstr(const mxArray *,char *,int);
extern int intmatch(unsigned int,const unsigned int[],int);
extern int strmatch(const char *,const char *[],int);
extern int getvalfromstr(const mxArray *matstr,const char *strs[],int nvals);

extern comedi_t* get_comedi_device_safe(const char *devname);

extern void matstruct_to_cmd(const mxArray *,comedi_cmd *);
extern void cmd_to_matstruct(const comedi_cmd *,mxArray *);
extern void matstruct_to_insn(const mxArray *,comedi_insnlist *);
extern void prep_chanspec(commat_chanspec *);
extern int do_chan_packing(unsigned int **,commat_chanspec *);
extern void comedi_cmd_mxFree(comedi_cmd *);
extern double mygetscalar(const mxArray *);

/* 
extern void capture_stderr();
extern char* stderrcpy(char *);
*/

extern const char *arefnames[];
extern const unsigned int arefvals[];
extern const int n_aref;
extern const char *unitnames[];
extern const unsigned int unitvals[];
extern const int n_unit;
extern const char *cmdfields[];
extern const int n_cmdfields;
extern const char *subdevnames[];
extern const unsigned int subdevvals[];
extern const int n_subdevnames;
extern const char *dio_direction_names[];
extern const unsigned int dio_direction_vals[];
extern const int n_dio_directions;
extern const char *insn_fields[];
extern const int n_insn_fields;
extern const char *insn_names[];
extern const unsigned int insn_vals[];
extern const int n_insn;

extern void dump_cmd(FILE *out,comedi_cmd *cmd);

