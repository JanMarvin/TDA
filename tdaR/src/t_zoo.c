#include "tda_rhooks.h"
/****************************************************************************/
/*  t_zoo.c                                                                 */
/*                                                                          */
/*      The functions in this module are adapted froom ZOO,                 */
/*      an archive package written by Rahul Dhesi.                          */
/*      Copyright is due to Rahul Dhesi.                                    */
/*      See the copyright statement in the ZOO package.                     */
/*                                                                          */
/*      We have taken mainly two sorts of functions from the ZOO package.   */
/*      1) Functions to read the header and directories of a ZOO archive.   */
/*      2) Functions for decompressing files from a ZOO archive.            */
/*         Type 1 - lzd compression                                         */
/*         Type 2 - lzh compression                                         */
/*                                                                          */
/*      All constants and descriptions for these functions are contained    */
/*      in this source. NOTE: t_uint16 MUST be 16 bit unsigned.             */
/*                                                                          */

#include "tda.h"
#include "t_gen.h"
#include "t_rzoo.h"
#include "tda_context.h"

#define t_uint16 unsigned short     /* must be 16 bit unsigned  !!!!        */

#define  IN_BUF_SIZE      8192 
#define  OUT_BUF_SIZE     8192 

/* MEM_BLOCK_SIZE must be no less than (2 * DICSIZ + MAXMATCH)
(see ar.h and lzh.h for values).  The buffer of this size will
also hold an input buffer of IN_BUF_SIZE and an output buffer
of OUT_BUF_SIZE.  FUDGE is a fudge factor, to keep some spare and
avoid off-by-one errors. */

#define  FUDGE    8
#define  MEM_BLOCK_SIZE (8192 + 8192 + 256 + 8)

#define  INBUFSIZ    ((unsigned int)(IN_BUF_SIZE - 10))  /* avoid obo errors */
#define  OUTBUFSIZ   (OUT_BUF_SIZE - 10)
#define  MAXBITS     13
#define  CLEAR       256        /* clear code */
#define  Z_EOF       257        /* end of file marker */
#define  FIRST_FREE  258        /* first free code */
#define  MAXMAX      8192       /* max code + 1 */

#define CHAR_BITS 8
#define UCHAR_MX 255
#define BITBUFSIZ (CHAR_BITS * sizeof bitbuf)
#define MATCHBIT   8            /* bits for MAXMATCH - THRESHOLD */
#define MAXMATCH 256            /* formerly F (not more than UCHAR_MX + 1) */
#define THRESHOLD  3            /* choose optimal value */
#define PERC_FLAG ((unsigned) 0x8000)
#define NC (UCHAR_MX + MAXMATCH + 2 - THRESHOLD)
                                /* alphabet = {0, 1, 2, ..., NC - 1} */
#define CBIT 9                  /* $\lfloor \log_2 NC \rfloor + 1$ */
#define CODE_BIT  16            /* codeword length */
#define DICBIT    13            /* 12(-lh4-) or 13(-lh5-) */
#define DICSIZ ((unsigned) 1 << DICBIT)

/* should be 0xFDC4A7DCUL but many c compilers don't recognize UL at end */
#define ZOOTAG ((unsigned TLONG) 0xFDC4A7DCL)          /* A random choice */
#define TEXT "ZOO 2.10 Archive.\032"         /* Header text for archive. */

/*  Structures: zoo_head, dirent */

#define  SIZ_TEXT    20         /* Size of header text */
#define  PATHSIZE   256         /* Max length of pathname */
#define  FNAMESIZE   13         /* Size of DOS filename */
#define  LFNAMESIZE 256         /* Size of long filename */
#define  MINZOOHSIZ  34         /* minimum size of archive header */
#define  SIZ_ZOOH    42         /* length of current archive header */

/* offsets of items within the canonical zoo archive header */
#define  TEXT_I       0         /* text in header */
#define  ZTAG_I      20         /* zoo tag */
#define  ZST_I       24         /* start offset */
#define  ZSTM_I      28         /* negative of start offset */
#define  MAJV_I      32         /* major version */
#define  MINV_I      33         /* minor version */
#define  HTYPE_I     34         /* archive header type */
#define  ACMTPOS_I   35         /* position of archive comment */
#define  ACMTLEN_I   39         /* length of archive comment */
#define  HVDATA_I    41         /* version data */

/* offsets of items within the canonical directory entry structure */
#define  SIZ_DIR     51         /* length of type 1 directory entry */
#define  SIZ_DIRL    56         /* length of type 2 directory entry */
#define  DTAG_I       0         /* tag within directory entry */
#define  DTYP_I       4         /* type of directory entry */
#define  PKM_I        5         /* packing method */
#define  NXT_I        6         /* pos'n of next directory entry */
#define  OFS_I       10         /* position (offset) of this file */
#define  DAT_I       14         /* DOS format date */
#define  TIM_I       16         /* DOS format time */
#define  CRC_I       18         /* CRC of this file */
#define  ORGS_I      20         /* original size */
#define  SIZNOW_I    24         /* size now */
#define  DMAJ_I      28         /* major version number */
#define  DMIN_I      29         /* minor version number */
#define  DEL_I       30         /* deleted or not */
#define  STRUC_I     31         /* file structure */
#define  CMT_I       32         /* comment [offset] */
#define  CMTSIZ_I    36         /* comment size */
#define  FNAME_I     38         /* filename */
#define  VARDIRLEN_I 51         /* length of var. direntry */
#define  TZ_I        53         /* timezone */
#define  DCRC_I      54         /* CRC of directory entry */
#define  FNM_SIZ     13         /* size of stored filename */

/* Offsets within variable part of directory entry */
#define  NAMLEN_I   (SIZ_DIRL + 0)
#define  DIRLEN_I   (SIZ_DIRL + 1)
#define  LFNAME_I   (SIZ_DIRL + 2)
#define  DIRNAME_I  LFNAME_I    /* plus length of filename */

/*
Total size of fixed plus variable directory recognized currently:
One byte each for dirlen and namlen, 256 each for long filename and
directory name, 2 for system id, 3 for file attributes, 1 for 
version flag, 2 for version number, plus a fudge factor of 5.
*/

#define  MAXDIRSIZE  (SIZ_DIRL+1+1+256+256+2+3+1+2+5)

struct zoo_head {
    char text[SIZ_TEXT];        /* archive header text */
    unsigned TLONG zoo_tag;      /* identifies archives */
    TLONG zoo_start;             /* where the archive's data starts */
    TLONG zoo_minus;             /* for consistency checking of zoo_start */
    unsigned char major_ver;
    unsigned char minor_ver;    /* minimum version to extract all files   */
    unsigned char type;         /* type of archive header */
    TLONG acmt_pos;              /* position of archive comment */
    unsigned int acmt_len;      /* length of archive comment */
    unsigned int vdata;         /* byte in archive;  data about versions */
};

struct dirent {
    unsigned TLONG zoo_tag;      /* tag -- redundancy check */
    unsigned char type;         /* type of directory entry.  always 1 for now */
    unsigned char packing_method;   /* 0 = no packing, 1 = normal LZW */
    TLONG next;                  /* pos'n of next directory entry */
    TLONG offset;                /* position of this file */
    unsigned int date;          /* DOS format date */
    unsigned int time;          /* DOS format time */
    unsigned int file_crc;      /* CRC of this file */
    TLONG org_size;
    TLONG size_now;
    unsigned char major_ver;
    unsigned char minor_ver;    /* minimum version needed to extract */
    unsigned char deleted;      /* will be 1 if deleted, 0 if not */
    unsigned char struc;        /* file structure if any */
    TLONG comment;               /* points to comment;  zero if none */
    unsigned int cmt_size;      /* length of comment, 0 if none */
    char fname[FNAMESIZE];      /* filename */

    int var_dir_len;            /* length of variable part of dir entry */
    unsigned char tz;           /* timezone where file was archived */
    unsigned int dir_crc;       /* CRC of directory entry */

   /* fields for variable part of directory entry follow */
    unsigned char namlen;       /* length of long filename */
    unsigned char dirlen;       /* length of directory name */
    char lfname[LFNAMESIZE];    /* long filename */
    char dirname[PATHSIZE];     /* directory name */
    unsigned int system_id;     /* Filesystem ID */
    unsigned TLONG fattr;        /* File attributes -- 24 bits */
    unsigned int vflag;         /* version flag bits -- one byte in archive */
    unsigned int version_no;    /* file version number if any */
};

int get_zoo(TDAContext *ctx);
int frd_zooh(TDAContext *ctx, struct zoo_head *zoo_header);
void b_to_zooh(TDAContext *ctx, struct zoo_head *zoo_header, unsigned char bytes[]);
int frd_dir(TDAContext *ctx, struct dirent *direntry);
void b_to_dir(TDAContext *ctx, struct dirent *direntry, unsigned char bytes[]);
TLONG to_long(TDAContext *ctx, unsigned char data[]);
int dbf_init(TDAContext *ctx, int opt);
int lzs(TDAContext *ctx, int fn);
int lzd(TDAContext *ctx, int fn);
void init_dtab(TDAContext *ctx);
void wr_dchar(TDAContext *ctx, int ch);
unsigned int rd_dcode(TDAContext *ctx);
void ad_dcode(TDAContext *ctx);
int lzh_decode(TDAContext *ctx, int fn);
void fillbuf(TDAContext *ctx, int n);
unsigned int getbits(TDAContext *ctx, int n);
void init_getbits(TDAContext *ctx);
void make_table(TDAContext *ctx, int nchar, unsigned char *bitlen, int tablebits, unsigned short *xtable);
unsigned int decode_c(TDAContext *ctx);
unsigned int decode_p(TDAContext *ctx);
void huf_decode_start(TDAContext *ctx);
void decode_start(TDAContext *ctx);
int decode(TDAContext *ctx, unsigned int count, unsigned char buffer[]);

/* ------------------------------------------------------------------------ */
/*  get_zoo     Read directory of ZOO archive. Return 0 if OK, -1 if error. */

int get_zoo(TDAContext *ctx)
{
    struct zoo_head zoo_header;             /* header for archive */
    struct dirent direntry;                 /* directory entry    */
    register int i;
    long next_ptr;  
    char *name,full[PATHSIZE + LFNAMESIZE + 1];

    if (frd_zooh(ctx, &zoo_header))           /* read the archive header */
        return(-1);

    if ((zoo_header.zoo_start + zoo_header.zoo_minus) != 0L) {
        printf1(ctx, "Bad archive header.\n");
        return(-1);
    }
    fseek(ctx->ZOOFd,zoo_header.zoo_start,0); 

    /*  Read directory entries of the archive file. If a match is found     */
    /*  with a file defined in the archive description file, the relevant   */
    /*  information is saved.                                               */

    while (1) {
        if (frd_dir(ctx, &direntry))
            return(-1);

        if (direntry.zoo_tag != ZOOTAG) {
            printf1(ctx, "Invalid archive header.\n");
            return(-1);
        }
        if (direntry.next == 0L)    
            break;

        if (!direntry.deleted) {        /* ignore deleted files */

            /* The name a member is known by: its long name when the
               entry carries one (zoo 2.1's type-2 extension), else the
               DOS-style short field; prefixed by the stored directory.
               A description file may give either the full stored path or,
               as before, the bare name. */

            if (direntry.namlen > 1)
                name = direntry.lfname;
            else
                name = direntry.fname;
            if (direntry.dirlen > 1)
                snprintf(full,sizeof(full),"%s/%s",direntry.dirname,name);
            else
                snprintf(full,sizeof(full),"%s",name);

            for (i = 0; i < ctx->ZANF; ++i) {
                if (!strcmp(ctx->ZAFNam[i],full) ||
                    !strcmp(ctx->ZAFNam[i],name)) {
                    if (ctx->ZAZOfs[i]) {
                        printf1(ctx, "Found more than one version of %s.\n",ctx->ZAFNam[i]);
                        return(-1);
                    }
                    ctx->ZAZOfs[i] = (int)direntry.offset;
                    ctx->ZAZTyp[i] = (char)direntry.packing_method;
                    ctx->ZAZSiz[i] = (int)direntry.org_size;
                }
            }
        }
        next_ptr = direntry.next;                 /* ptr to next dir entry  */
        fseek(ctx->ZOOFd,next_ptr,0);                  /* seek to next dir entry */
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  frd_zooh.                                                               */
/*      Reads the header of a Zoo archive in a machine-independent manner.  */
/*      Return 0 if OK, -1 if error.                                        */

int frd_zooh(TDAContext *ctx, struct zoo_head *zoo_header)
{
    int status;
    unsigned char bytes[SIZ_ZOOH];       /* canonical header representation */

    status = (int)(fread((char *)bytes,1,SIZ_ZOOH,ctx->ZOOFd));

    b_to_zooh(ctx, zoo_header,bytes);             /* convert array to structure */

    if (status < MINZOOHSIZ)
        return (-1);
    else
        return (0);
}

/* ------------------------------------------------------------------------ */
/*  b_to_zooh.  Converts an array of unsigned char to a zoo_header struct.  */
/*              We only convert entries that are actually needed for        */
/*              this program.                                               */

void b_to_zooh(TDAContext *ctx, struct zoo_head *zoo_header, unsigned char bytes[])
{
    int i;

    for (i = 0; i < SIZ_TEXT; i++)                     /* copy text */
        zoo_header->text[i] = (char)(bytes[TEXT_I + i]);

    zoo_header->zoo_tag = (unsigned int)(to_long(ctx, &bytes[ZTAG_I]));     /* copy zoo_tag */
    zoo_header->zoo_start = to_long(ctx, &bytes[ZST_I]);    /* copy zoo_start */
    zoo_header->zoo_minus = to_long(ctx, &bytes[ZSTM_I]);
}

/* ------------------------------------------------------------------------ */
/*  frd_dir.   Reads a directory entry in a machine-independent manner.     */
/*             Return 0 if OK, -1 if error.                                 */

int frd_dir(TDAContext *ctx, struct dirent *direntry)
{
    int status;
    unsigned char bytes[MAXDIRSIZE]; /* big enough to hold variable part too */

    /* To simplify things, we read the maximum possible size of the
       directory entry including the variable size and discard what is not
       needed */

    memset(bytes,0,sizeof(bytes));
    status = (int)(fread((char *)bytes,1,MAXDIRSIZE,ctx->ZOOFd));
    if (status < SIZ_DIR)
        return (-1);
    b_to_dir(ctx, direntry,bytes);
    return (0);
}

/* ------------------------------------------------------------------------ */
/*  b_to_dir.  Converts bytes to directory entry structure.                 */

void b_to_dir(TDAContext *ctx, struct dirent *direntry, unsigned char bytes[])
{
    int i;

    direntry->zoo_tag = (unsigned int)(to_long(ctx, &bytes[DTAG_I]));
    direntry->type = bytes[DTYP_I];
    direntry->packing_method = bytes[PKM_I];
    direntry->next = to_long(ctx, &bytes[NXT_I]);
    direntry->offset = to_long(ctx, &bytes[OFS_I]);
    direntry->org_size = to_long(ctx, &bytes[ORGS_I]);
    direntry->deleted = bytes[DEL_I];

    for (i = 0; i < FNM_SIZ; i++)
        direntry->fname[i] = (char)(bytes[FNAME_I + i]);
    direntry->fname[FNM_SIZ - 1] = '\0';

    /* Type-2 entries (zoo 2.1) carry a variable part after the 56-byte
       fixed part: a byte each for the lengths of the long file name and
       the directory name (both counted with their terminating NUL), then
       the two names.  Type-1 entries have neither. */

    direntry->var_dir_len = 0;
    direntry->namlen = direntry->dirlen = 0;
    direntry->lfname[0] = direntry->dirname[0] = '\0';
    if (direntry->type == 2) {
        direntry->var_dir_len = (int)(bytes[VARDIRLEN_I] |
                                      (bytes[VARDIRLEN_I + 1] << 8));
        if (direntry->var_dir_len > MAXDIRSIZE - SIZ_DIRL)
            direntry->var_dir_len = MAXDIRSIZE - SIZ_DIRL;
        if (direntry->var_dir_len > 0)
            direntry->namlen = bytes[NAMLEN_I];
        if (direntry->var_dir_len > 1)
            direntry->dirlen = bytes[DIRLEN_I];
        if (direntry->namlen + direntry->dirlen + 2 > direntry->var_dir_len)
            direntry->namlen = direntry->dirlen = 0;
        for (i = 0; i < direntry->namlen; i++)
            direntry->lfname[i] = (char)(bytes[LFNAME_I + i]);
        direntry->lfname[direntry->namlen] = '\0';
        for (i = 0; i < direntry->dirlen; i++)
            direntry->dirname[i] = (char)(bytes[DIRNAME_I + direntry->namlen + i]);
        direntry->dirname[direntry->dirlen] = '\0';
    }
}

/* ------------------------------------------------------------------------ */
/*  to_long. Converts four consecutive bytes, in order of increasing        */
/*           significance, to a long integer.  It is used to make Zoo       */
/*           independent of the byte order of the system.                   */

TLONG to_long(TDAContext *ctx, unsigned char data[])
{
    (void)ctx;        /* unused: the signature is shared */
    return (TLONG) ((unsigned TLONG) data[0] | ((unsigned TLONG) data[1] << 8) |
         ((unsigned TLONG) data[2] << 16) | ((unsigned TLONG) data[3] << 24));
}

/* ------------------------------------------------------------------------ */
/*  Functions to get data from a ZOO archive. There are two possible        */
/*  methods, recorded for each file in ZAZTyp[fn].                          */
/*                                                                          */
/*  Typ 1   -   lzd compression, function: lzd                              */
/*  Typ 2   -   lzh compression, function: lzh_decode                       */
/*                                                                          */

/*  Static storage and definitions used by lzd */

#define STACKSIZE 4000      
struct tabentry {
   unsigned int next;
   char z_ch;
};

#define  push(x)  {  \
            stack[stack_pointer++] = (unsigned int)(x); \
            if (stack_pointer >= STACKSIZE)    \
                printfe(ctx, "STACK OVERFLOW IN LZD().\n"); \
         }
#define  pop() (stack[--stack_pointer])

struct tabentry *table;
static unsigned int stack_pointer;
static unsigned int *stack;
static unsigned int cur_code;
static unsigned int old_code;
static unsigned int in_code;
static unsigned int free_code;
static int nbits;
static unsigned int max_code;
static char fin_char;
static char kk;
static unsigned int bit_offset;
static unsigned int masks[] = {
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0x1ff, 0x3ff, 0x7ff, 0xfff, 0x1fff
};

/*  Definition of a macro assert() that causes an assertion error if the 
    assertion fails.
    Conditional compilation:
    If NDEBUG is defined then
      assert() is defined as null so all assertions vanish
    else
      if __FILE__ and __LINE__ are defined then
         assertions print message including filename and line number
      else
         assertions print a message but not the filename and line number
      endif
   endif
*/

#ifdef NDEBUG
# define assert(E)
#else
#undef LINE_FILE
#ifdef __LINE__
#ifdef __FILE__
#define LINE_FILE
#endif
#endif
#ifdef LINE_FILE
#undef LINE_FILE
#define assert(E) \
   { if (!(E)) \
       printf1(ctx, "\nassertion error \n"); \
   }
#else
#define assert(E) \
   { if (!(E)) \
       printf1(ctx, " \nassertion error \n"); \
   }
#endif
#endif /* NDEBUG */

/* static storage and definitions used by lzh_decode */

static char decoded;

#define NP (DICBIT + 1)
#define NT (CODE_BIT + 3)
#define PBIT 4          /* smallest integer such that (1U << PBIT) > NP */
#define TBIT 5          /* smallest integer such that (1U << TBIT) > NT */
#if NT > NP
# define NPT NT
#else
# define NPT NP
#endif

static unsigned short *left; 
static unsigned short *right;
static unsigned char *c_len; 
static unsigned char *pt_len;
static unsigned int blocksize;
static unsigned short *c_table;
static unsigned short *pt_table;

t_uint16 bitbuf;
static unsigned int subbitbuf;
static int bitcount;
static int jj;    

/* ------------------------------------------------------------------------ */
/*  dbf_init(opt)   if opt != 0 allocate, otherwise free memory.            */
/*                  if ZAAlloc is set, then initialization without alloc.   */
/*                  return 0 if OK, -1 if error.                            */

int dbf_init(TDAContext *ctx, int opt)
{
    int err = -1;

    if (opt == 0) {
        err = 0;
        if (ctx->ZAAlloc)
            goto DBFin10;
        else
            goto DBFin;
    }

    ctx->ZAEOF = 0;
    ctx->Out_Buf_Cnt = 0;
    ctx->ZAInit = 0;
    if (ctx->ZAAlloc) { 
        ctx->In_Buf_Adr = ctx->Out_Buf_Adr + OUT_BUF_SIZE + (FUDGE/2);
        return(0);
    }
    if (!(c_table = (unsigned short *) calloc(4096,sizeof(unsigned short))))
        goto DBFin;
    memrq(ctx, 4096,sizeof(unsigned short));
   
    if (!(pt_table = (unsigned short *) calloc(256,sizeof(unsigned short))))
        goto DBFin1;
    memrq(ctx, 256,sizeof(unsigned short));
  
    if (!(c_len = (unsigned char *) calloc(NC,sizeof(unsigned char))))  
        goto DBFin2;
    memrq(ctx, NC,sizeof(char));
  
    if (!(pt_len = (unsigned char *) calloc(NPT,sizeof(unsigned char))))  
        goto DBFin3;
    memrq(ctx, NPT,sizeof(unsigned char));
  
    if (!(left = (unsigned short *) calloc(2 * NC - 1,sizeof(unsigned short))))
        goto DBFin4;
    memrq(ctx, 2 * NC - 1,sizeof(unsigned short));
  
    if (!(right = (unsigned short *) calloc(2 * NC - 1,sizeof(unsigned short))))
        goto DBFin5;
    memrq(ctx, 2 * NC - 1,sizeof(unsigned short));
  
    if (!(ctx->Out_Buf_Adr = (char *) calloc(MEM_BLOCK_SIZE,sizeof(char))))  
        goto DBFin6;
    memrq(ctx, MEM_BLOCK_SIZE,sizeof(char));
   
    if (!(ctx->ZABuf = (char *) calloc((size_t)(ctx->ZABLen + 4),sizeof(char))))
        goto DBFin7;
    memrq(ctx, ctx->ZABLen + 4,sizeof(char));

    if (!(stack = (unsigned int *) calloc(STACKSIZE + 20,sizeof(unsigned int))))  
        goto DBFin8;
    memrq(ctx, STACKSIZE + 20,sizeof(unsigned int));

    if (!(table = (struct tabentry *) calloc(MAXMAX + 10,sizeof(struct tabentry))))  
        goto DBFin9;
    memrq(ctx, MAXMAX + 10,sizeof(struct tabentry));

    ctx->In_Buf_Adr = ctx->Out_Buf_Adr + OUT_BUF_SIZE + (FUDGE/2);
    ctx->ZAAlloc = 1;
    return(0);

DBFin10:
    free((char *)table);
    memrq(ctx, -MAXMAX - 10,sizeof(struct tabentry));
DBFin9:
    free((char *)stack);
    memrq(ctx, -STACKSIZE - 20,sizeof(unsigned int));
DBFin8:
    free((char *)ctx->ZABuf);
    memrq(ctx, -ctx->ZABLen - 4,sizeof(char));
DBFin7:
    free((char *)ctx->Out_Buf_Adr);
    memrq(ctx, -MEM_BLOCK_SIZE,sizeof(char));
DBFin6:
    free((char *)right);
    memrq(ctx, -2 * NC + 1,sizeof(unsigned short));
DBFin5:
    free((char *)left);
    memrq(ctx, -2 * NC + 1,sizeof(unsigned short));
DBFin4:
    free((char *)pt_len);
    memrq(ctx, -NPT,sizeof(unsigned char));
DBFin3: 
    free((char *)c_len);
    memrq(ctx, -NC,sizeof(char));
DBFin2:
    free((char *)pt_table);
    memrq(ctx, -256,sizeof(unsigned short));
DBFin1:
    free((char *)c_table);
    memrq(ctx, -4096,sizeof(unsigned short));
DBFin:
    ctx->ZAAlloc = 0;
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  lzd.    Decompression function from ZOO, packing method 1. The function */
/*          is modified to decompress several files simultaneously. The     */
/*          decompressed data are written into Out_Buf_Adr. The number      */
/*          of bytes writte into this buffer is recorded in Out_Buf_Cnt     */
/*          The function returns, if LZD_COUNT bytes are written. This      */
/*          value must be smaller than OUTBUFSIZ to avoid overflow.         */

#define LZD_COUNT INBUFSIZ/2

/*  Packing method 0: the member is stored, so "decoding" is a copy.  TDA
    never wrote archives and so never met one -- every archive it saw came
    from Dhesi's zoo, which compresses -- but a stored member is perfectly
    legal ZOO and every other reader takes it.  Mirrors lzd()'s own
    structure: seek on first call, then hand back one buffer at a time and
    set ZAEOF when the member is exhausted.  */
int lzs(TDAContext *ctx, int fn)
{
    long want;
    size_t got;

    ctx->Out_Buf_Cnt = 0;
    if (ctx->ZAEOF)
        return(0);
    if (!ctx->ZAInit) {
        ctx->FPtr = ctx->ZAZOfs[fn];
        if (fseek(ctx->ZOOFd,ctx->FPtr,0) == -1L) {
            printfe(ctx, "SEEK ERROR IN LZS. CAN'T CONTINUE.\n");
            gerr_exit(ctx, 18);
        }
        ctx->ZALeft = ctx->ZAZSiz[fn];
        ctx->ZAInit = 1;
    }
    want = ctx->ZALeft;
    if (want > OUT_BUF_SIZE)
        want = OUT_BUF_SIZE;
    if (want <= 0) {
        ctx->ZAEOF = 1;
        return(0);
    }
    got = fread(ctx->Out_Buf_Adr,1,(size_t)want,ctx->ZOOFd);
    if (got == 0) {
        ctx->ZAEOF = 1;
        return(0);
    }
    ctx->ZALeft -= (int)got;
    ctx->Out_Buf_Cnt = (int)got;
    if (ctx->ZALeft <= 0)
        ctx->ZAEOF = 1;
    return(ctx->Out_Buf_Cnt);
}

int lzd(TDAContext *ctx, int fn)
{
    ctx->Out_Buf_Cnt = 0;
    if (ctx->ZAEOF)
        return(0);

    if (!ctx->ZAInit) {
        ctx->FPtr = ctx->ZAZOfs[fn];
        if (fseek(ctx->ZOOFd,ctx->FPtr,0) == -1L) {
            printfe(ctx, "SEEK ERROR IN LZD. CAN't CONTINUE.\n");
            gerr_exit(ctx, 18);
        }
        nbits = 9;
        max_code = 512;
        free_code = FIRST_FREE;
        stack_pointer = 0;
        bit_offset = 0;
        /* fread() returns size_t and never -1: the == -1 this replaces
           came from read() and made the error branch dead code.  A short
           read at end of archive is legitimate, so the error is ferror(). */
        if (fread(ctx->In_Buf_Adr,1,INBUFSIZ,ctx->ZOOFd) < (size_t)INBUFSIZ &&
            ferror(ctx->ZOOFd)) {
            printfe(ctx, "READ ERROR IN LZD. CAN'T CONTINUE.\n");
            gerr_exit(ctx, 19);
        }
        ctx->FPtr = ftell(ctx->ZOOFd);
        init_dtab(ctx);           /* initialize table */
        ctx->ZAInit = 1;
    }
loop:
    if ((unsigned int)ctx->Out_Buf_Cnt >= LZD_COUNT)
        return(ctx->Out_Buf_Cnt);

    cur_code = rd_dcode(ctx);

goteof:         /* special case for CLEAR then Z_EOF, for 0-length files */

    if (cur_code == Z_EOF) {
        ctx->ZAEOF = 1;
        return(ctx->Out_Buf_Cnt);
    }
    assert(nbits >= 9 && nbits <= 13);
    if (cur_code == CLEAR) {
        init_dtab(ctx);
        cur_code = rd_dcode(ctx);
        old_code = cur_code;
        kk = (char)cur_code;
        fin_char = (char)kk;

        if (cur_code == Z_EOF)    /* special case for 0-length files */
            goto goteof;
        wr_dchar(ctx, kk);
        goto loop;
    }
    in_code = cur_code;
    if (cur_code >= free_code) {          /* if code not in table (k<w>k<w>k) */
        cur_code = old_code;              /* previous code becomes current */
        push(fin_char);
    }
    while (cur_code > 255) {                        /* if code, not character */
        push(table[cur_code].z_ch);                     /* push suffix char */
        cur_code = table[cur_code].next;               /* <w> := <w>.code */
    }
    assert(nbits >= 9 && nbits <= 13);
    fin_char = (char)cur_code;
            kk = fin_char;
    push(kk);
    while (stack_pointer != 0) {
        wr_dchar(ctx, (int)pop());
    }
    assert(nbits >= 9 && nbits <= 13);
    ad_dcode(ctx);
    old_code = in_code;
    assert(nbits >= 9 && nbits <= 13);
    goto loop;
}

/* ------------------------------------------------------------------------ */
/*  init_dtab.  Initialize table.                                           */

void init_dtab(TDAContext *ctx)
{
    (void)ctx;        /* unused: the signature is shared */
    nbits = 9;
    max_code = 512;
    free_code = FIRST_FREE;
}

/* ------------------------------------------------------------------------ */
/*  wr_dchar                                                                */

void wr_dchar(TDAContext *ctx, int ch)
{
    if (ctx->Out_Buf_Cnt >= OUTBUFSIZ) {      /* if buffer full */
        printfe(ctx, "LZD OUTPUT BUF OVERFLOW. CAN'T CONTINUE.\n");
        gerr_exit(ctx, 20);
    }
    assert(ctx->Out_Buf_Cnt < OUTBUFSIZ);
    ctx->Out_Buf_Adr[ctx->Out_Buf_Cnt++] = (char)(ch);            /* store character */
}

/* ------------------------------------------------------------------------ */
/*  rd_dcode. Reads a code from the input (compressed) file and returns     */
/*            its value.                                                    */

unsigned int rd_dcode(TDAContext *ctx)
{
    register char *ptra, *ptrb;                 /* miscellaneous pointers */
    unsigned int word;                          /* first 16 bits in buffer */
    unsigned int byte_offset;
    char nextch;                                /* next 8 bits in buffer */
    unsigned int ofs_inbyte;                    /* offset within byte */
    ofs_inbyte = bit_offset % 8;
    byte_offset = bit_offset / 8;
    bit_offset = bit_offset + (unsigned int)nbits;
    assert(nbits >= 9 && nbits <= 13);

    if (byte_offset >= INBUFSIZ - 5) {
        int space_left;

        assert(byte_offset >= INBUFSIZ - 5);
        bit_offset = ofs_inbyte + (unsigned int)nbits;
        space_left = (int)(INBUFSIZ - byte_offset);
        ptrb = byte_offset + ctx->In_Buf_Adr;          /* point to char */
        ptra = ctx->In_Buf_Adr;

        /* we now move the remaining characters down buffer beginning */
        while (space_left > 0) {
            *ptra++ = *ptrb++;
            space_left--;
        }
        assert(ptra - ctx->In_Buf_Adr == ptrb - (ctx->In_Buf_Adr + byte_offset));
        assert(space_left == 0);
        if (fseek(ctx->ZOOFd,ctx->FPtr,0) == -1L) {
            printfe(ctx, "SEEK ERROR IN LZD. CAN'T CONTINUE.\n");
            gerr_exit(ctx, 21);
        }
        if (fread(ptra,1,(size_t)byte_offset,ctx->ZOOFd) < (size_t)byte_offset &&
            ferror(ctx->ZOOFd)) {                    /* see above */
            printfe(ctx, "READ ERROR IN LZD. CAN'T CONTINUE.\n");
            gerr_exit(ctx, 22);
        }
        ctx->FPtr = ftell(ctx->ZOOFd);
        byte_offset = 0;
    }
    ptra = byte_offset + ctx->In_Buf_Adr;

    /* NOTE: "word = *((int *) ptra)" would not be independent of byte order. */

    word = (unsigned char) *ptra; ptra++;
    word = word | (unsigned int)(((unsigned char) *ptra) << 8); ptra++;
    nextch = *ptra;
    if (ofs_inbyte != 0) {
        /* shift nextch right by ofs_inbyte bits */
        /* and shift those bits right into word; */
        word = (word >> ofs_inbyte) | (((unsigned int)nextch) << (16-ofs_inbyte));
    }
    return (word & masks[nbits]); 
}

/* ------------------------------------------------------------------------ */
/*  ad_dcode.  Adds a code to table                                         */

void ad_dcode(TDAContext *ctx)
{
    (void)ctx;        /* unused: the signature is shared */
    assert(nbits >= 9 && nbits <= 13);
    assert(free_code <= MAXMAX+1);
    table[free_code].z_ch = kk;                        /* save suffix char */
    table[free_code].next = old_code;                  /* save prefix code */
    free_code++;
    assert(nbits >= 9 && nbits <= 13);
    if (free_code >= max_code) {
        if (nbits < MAXBITS) {
            nbits++;
            assert(nbits >= 9 && nbits <= 13);
            max_code = max_code << 1;            /* double max_code */
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  lzd_decode. Decompression function from ZOO, packing method 2. The      */
/*       function is modified to decompress several files simultaneously.   */
/*       The decompressed data are written into Out_Buf_Adr. The number     */
/*       of bytes writte into this buffer is recorded in Out_Buf_Cnt        */
/*       The function returns, if LZD_COUNT bytes are written. This value   */
/*       must be smaller than OUTBUFSIZ to avoid overflow.                  */

int lzh_decode(TDAContext *ctx, int fn)
{

    if (ctx->ZAEOF) {
        ctx->Out_Buf_Cnt = 0;
        return(0);
    }
    if (!ctx->ZAInit) {
        ctx->FPtr = ctx->ZAZOfs[fn];
        if (fseek(ctx->ZOOFd,ctx->FPtr,0) == -1L) {
            printfe(ctx, "SEEK ERROR IN LZD. CAN'T CONTINUE.\n");
            gerr_exit(ctx, 23);
        }
        decode_start(ctx);
        ctx->ZAInit = 1;
    }
    else {
        if (fseek(ctx->ZOOFd,ctx->FPtr,0) == -1L) {
            printfe(ctx, "SEEK ERROR IN LZD. CAN'T CONTINUE.\n");
            gerr_exit(ctx, 24);
        }
    }
    if (!decoded)  
        ctx->Out_Buf_Cnt = decode(ctx, (unsigned int)DICSIZ,(unsigned char *)ctx->Out_Buf_Adr);

    ctx->FPtr = ftell(ctx->ZOOFd);
    if (decoded)
        ctx->ZAEOF = 1;
    return(ctx->Out_Buf_Cnt);
}

/* ------------------------------------------------------------------------ */
/*  fillbuf    Shift bitbuf n bits left, read n bits                        */

void fillbuf(TDAContext *ctx, int n)
{
    bitbuf <<= n;
    while (n > bitcount) {
        bitbuf = (unsigned short)(bitbuf | (subbitbuf << (n -= bitcount)));
        if (feof(ctx->ZOOFd))
            subbitbuf = 0;
        else
            subbitbuf = (unsigned char)getc(ctx->ZOOFd);
        bitcount = CHAR_BITS;
    }
    bitbuf = (unsigned short)(bitbuf | (subbitbuf >> (bitcount -= n)));
}

/* ------------------------------------------------------------------------ */
/*  getbits                                                                 */

unsigned int getbits(TDAContext *ctx, int n)
{
    unsigned int x;

    x = bitbuf >> (unsigned int)((int)BITBUFSIZ - n); 
    fillbuf(ctx, n);
    return(x);
}

/* ------------------------------------------------------------------------ */
/*  init_getbits                                                            */

void init_getbits(TDAContext *ctx)
{
    bitbuf = 0; 
    subbitbuf = 0;
    bitcount = 0;
    fillbuf(ctx, BITBUFSIZ);
}

/* ------------------------------------------------------------------------ */
/*  maketable. Make table for decoding                                      */

void make_table(TDAContext *ctx, int nchar, unsigned char *bitlen, int tablebits, unsigned short *xtable)
{
    unsigned short count[17], weight[17], start[18], *p;
    unsigned int i, k, len, ch, jutbits, avail, nextcode, mask;

    for (i = 1; i <= 16; i++)
        count[i] = 0;
    for (i = 0; i < (unsigned int)nchar; i++) 
        count[bitlen[i]]++;

    start[1] = 0;
    for (i = 1; i <= 16; i++)
        start[i + 1] = (unsigned short)(start[i] + (count[i] << (16 - i)));
    if (start[17] != (unsigned short)((unsigned int) 1 << 16)) {
        printfe(ctx, "LZH BAD DECODE TABLE. CAN'T CONTINUE.\n");
        gerr_exit(ctx, 25);
    }
    jutbits = (unsigned int)(16 - tablebits);
    for (i = 1; i <= (unsigned int)tablebits; i++) {
        start[i] >>= jutbits;
        weight[i] = (unsigned short)(1U << ((unsigned int)tablebits - i));
    }
    while (i <= 16) {
        weight[i] = (unsigned short)((unsigned int) 1 << (16 - i));
        i++;
    }
    i = start[tablebits + 1] >> jutbits;
    if (i != (unsigned short)((unsigned int) 1 << 16)) {
        k = 1 << tablebits;
        while (i != k)
            xtable[i++] = 0;
    }
    avail = (unsigned int)(nchar);
    mask = (unsigned int) 1 << (15 - tablebits);
    for (ch = 0; ch < (unsigned int)nchar; ch++) {
        if ((len = bitlen[ch]) == 0)
            continue;
        nextcode = start[len] + weight[len];
        if (len <= (unsigned int)tablebits) {
            for (i = start[len]; i < nextcode; i++) xtable[i] = (unsigned short)ch;
        }
        else {
            k = start[len];
            p = &xtable[k >> jutbits];
            i = len - (unsigned int)tablebits;
            while (i != 0) {
                if (*p == 0) {
                    right[avail] = left[avail] = 0;
                    *p = (unsigned short)(avail++);
                }
                if (k & mask)
                    p = &right[*p];
                else    
                    p = &left[*p];
                k <<= 1;  i--;
            }
            *p = (unsigned short)(ch);
        }
        start[len] = (unsigned short)(nextcode);
    }
}

/* ------------------------------------------------------------------------ */
/*  read_pt_len                                                             */

static void read_pt_len(TDAContext *ctx, int nn, int nbit, int i_special)
{
    int i, c, n;
    unsigned int mask;

    n = (int)(getbits(ctx, nbit));
    if (n == 0) {
        c = (int)(getbits(ctx, nbit));
        for (i = 0; i < nn; i++) pt_len[i] = 0;
        for (i = 0; i < 256; i++) pt_table[i] = (unsigned short)c;
    }
    else {
        i = 0;
        while (i < n) {
            c = bitbuf >> (BITBUFSIZ - 3);
            if (c == 7) {
                mask = (unsigned int) 1 << (BITBUFSIZ - 1 - 3);
                while (mask & bitbuf) {  mask >>= 1;  c++;  }
            }
            fillbuf(ctx, (c < 7) ? 3 : c - 3);
            pt_len[i++] = (unsigned char)(c);
            if (i == i_special) {
                c = (int)(getbits(ctx, 2));
                while (--c >= 0) pt_len[i++] = 0;
            }
        }
        while (i < nn) pt_len[i++] = 0;
        make_table(ctx, nn,pt_len,8,pt_table);
    }
}

/* ------------------------------------------------------------------------ */
/*  read_c_len                                                              */

static void read_c_len(TDAContext *ctx)
{
    int i, c, n;
    unsigned int mask;
    n = (int)(getbits(ctx, CBIT));
    if (n == 0) {
        c = (int)(getbits(ctx, CBIT));
        for (i = 0; i < NC; i++) c_len[i] = 0;
        for (i = 0; i < 4096; i++) c_table[i] = (unsigned short)c;
    }
    else {
        i = 0;
        while (i < n) {
            c = pt_table[bitbuf >> (BITBUFSIZ - 8)];
            if (c >= NT) {
                mask = (unsigned int) 1 << (BITBUFSIZ - 1 - 8);
                do {
                    if (bitbuf & mask)
                        c = right[c];
                    else
                        c = left[c];
                    mask >>= 1;
                } while (c >= NT);
            }
            fillbuf(ctx, (int) pt_len[c]);
            if (c <= 2) {
                if      (c == 0) c = 1;
                else if (c == 1) c = (int)getbits(ctx, 4) + 3;
                else             c = (int)(getbits(ctx, CBIT) + 20);
                while (--c >= 0) c_len[i++] = 0;
            } else c_len[i++] = (unsigned char)(c - 2);
        }
        while (i < NC) c_len[i++] = 0;
        make_table(ctx, NC,c_len,12,c_table);
    }
}

/* ------------------------------------------------------------------------ */
/*  decode_c                                                                */

unsigned int decode_c(TDAContext *ctx)
{
    unsigned int j, mask;
    if (blocksize == 0) {
        blocksize = getbits(ctx, 16);
        if (blocksize == 0) {
            decoded = 1;
            return 0;
        }
        read_pt_len(ctx, NT,TBIT,3);
        read_c_len(ctx);
        read_pt_len(ctx, NP,PBIT,-1);
    }
    blocksize--;
    j = c_table[bitbuf >> (BITBUFSIZ - 12)];
    if (j >= NC) {
        mask = (unsigned int) 1 << (BITBUFSIZ - 1 - 12);
        do {
            if (bitbuf & mask)
                j = right[j];
            else
                j = left[j];
            mask >>= 1;
        } while (j >= NC);
    }
    fillbuf(ctx, (int) c_len[j]);
    return(j);
}

/* ------------------------------------------------------------------------ */
/*  decode_p                                                                */

unsigned int decode_p(TDAContext *ctx)
{
    unsigned int j, mask;
    j = pt_table[bitbuf >> (BITBUFSIZ - 8)];
    if (j >= NP) {
        mask = (unsigned int) 1 << (BITBUFSIZ - 1 - 8);
        do {
            if (bitbuf & mask)
                j = right[j];
            else
                j = left[j];
            mask >>= 1;
        } while (j >= NP);
    }
    fillbuf(ctx, (int)pt_len[j]);
    if (j != 0) j = ((unsigned int) 1 << (j - 1)) + getbits(ctx, (int) (j - 1));
    return(j);
}

/* ------------------------------------------------------------------------ */
/*  huf_decode_start                                                        */

void huf_decode_start(TDAContext *ctx)
{
    init_getbits(ctx); 
    blocksize = 0;
}

/* ------------------------------------------------------------------------ */
/*  decode_start.                                                           */

void decode_start(TDAContext *ctx)
{
    huf_decode_start(ctx);
    jj = 0;
    decoded = 0;
}

/* ------------------------------------------------------------------------ */
/*  decode(). returns no. of chars decoded                                  */
/*                                                                          */
/*  The calling function must keep the number of bytes to be processed.     */
/*  This function decodes either 'count' bytes or 'DICSIZ' bytes, whichever */
/*  is smaller, into the array 'buffer[]' of size 'DICSIZ' or more.         */
/*  Call decode_start() once for each new file before calling this function */

int decode(TDAContext *ctx, unsigned int count, unsigned char buffer[])
{
    unsigned int r, c;
    r = 0;  

    while (--jj >= 0) {
        buffer[r] = buffer[ctx->s_decode_i];
        ctx->s_decode_i = (ctx->s_decode_i + 1) & (DICSIZ - 1);
        if (++r == count)
            return ((int)(r));
    }
    for ( ; ; ) {
        c = decode_c(ctx);
        if (decoded)  
            return ((int)(r));
        if (c <= UCHAR_MX) {
            buffer[r] = (unsigned char)(c);
            if (++r == count)    
                return ((int)(r));
        }
        else {
            jj = (int)(c - (UCHAR_MX + 1 - THRESHOLD));
            ctx->s_decode_i = (r - decode_p(ctx) - 1) & (DICSIZ - 1);
            while (--jj >= 0) {
                buffer[r] = buffer[ctx->s_decode_i];
                ctx->s_decode_i = (ctx->s_decode_i + 1) & (DICSIZ - 1);
                if (++r == count)     
                    return ((int)(r));
            }
        }
    }
}




void tda_reset_t_zoo(void)
{
    stack_pointer = cur_code = old_code = in_code = free_code = 0;
    nbits = 0; max_code = 0; fin_char = kk = 0; bit_offset = 0;
    stack = NULL; decoded = 0;
    left = right = NULL; c_len = pt_len = NULL; blocksize = 0;
    c_table = pt_table = NULL; subbitbuf = 0; bitcount = 0; jj = 0;
}
