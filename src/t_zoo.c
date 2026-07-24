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

#define  INBUFSIZ    (IN_BUF_SIZE - 10)   /* avoid obo errors */
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

int get_zoo(void);
int frd_zooh(struct zoo_head *zoo_header);
void b_to_zooh(struct zoo_head *zoo_header, unsigned char bytes[]);
int frd_dir(struct dirent *direntry);
void b_to_dir(struct dirent *direntry, unsigned char bytes[]);
TLONG to_long(unsigned char data[]);
int dbf_init(int opt);
int lzd(int fn);
void init_dtab(void);
void wr_dchar(int ch);
unsigned int rd_dcode(void);
void ad_dcode(void);
int lzh_decode(int fn);
void fillbuf(int n);
unsigned int getbits(int n);
void init_getbits(void);
void make_table(int nchar, unsigned char *bitlen, int tablebits, unsigned short *xtable);
unsigned int decode_c(void);
unsigned int decode_p(void);
void huf_decode_start(void);
void decode_start(void);
int decode(unsigned int count, unsigned char buffer[]);

/* ------------------------------------------------------------------------ */
/*  get_zoo     Read directory of ZOO archive. Return 0 if OK, -1 if error. */

int get_zoo(void)
{
    struct zoo_head zoo_header;             /* header for archive */
    struct dirent direntry;                 /* directory entry    */
    register int i;
    long next_ptr;  

    if (frd_zooh(&zoo_header))           /* read the archive header */
        return(-1);

    if ((zoo_header.zoo_start + zoo_header.zoo_minus) != 0L) {
        printf1("Bad archive header.\n");
        return(-1);
    }
    fseek(ZOOFd,zoo_header.zoo_start,0); 

    /*  Read directory entries of the archive file. If a match is found     */
    /*  with a file defined in the archive description file, the relevant   */
    /*  information is saved.                                               */

    while (1) {
        if (frd_dir(&direntry))
            return(-1);

        if (direntry.zoo_tag != ZOOTAG) {
            printf1("Invalid archive header.\n");
            return(-1);
        }
        if (direntry.next == 0L)    
            break;

        if (!direntry.deleted) {        /* ignore deleted files */

            /* check for matching files */

            for (i = 0; i < ZANF; ++i) {
                if (!strcmp(ZAFNam[i],direntry.fname)) {
                    if (ZAZOfs[i]) {
                        printf1("Found more than one version of %s.\n",ZAFNam[i]);
                        return(-1);
                    }
                    ZAZOfs[i] = (int)direntry.offset;
                    ZAZTyp[i] = (char)direntry.packing_method;
                    ZAZSiz[i] = (int)direntry.org_size;
                }
            }
        }
        next_ptr = direntry.next;                 /* ptr to next dir entry  */
        fseek(ZOOFd,next_ptr,0);                  /* seek to next dir entry */
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  frd_zooh.                                                               */
/*      Reads the header of a Zoo archive in a machine-independent manner.  */
/*      Return 0 if OK, -1 if error.                                        */

int frd_zooh(struct zoo_head *zoo_header)
{
    int status;
    unsigned char bytes[SIZ_ZOOH];       /* canonical header representation */

    status = fread((char *)bytes,1,SIZ_ZOOH,ZOOFd);

    b_to_zooh(zoo_header,bytes);             /* convert array to structure */

    if (status < MINZOOHSIZ)
        return (-1);
    else
        return (0);
}

/* ------------------------------------------------------------------------ */
/*  b_to_zooh.  Converts an array of unsigned char to a zoo_header struct.  */
/*              We only convert entries that are actually needed for        */
/*              this program.                                               */

void b_to_zooh(struct zoo_head *zoo_header, unsigned char bytes[])
{
    int i;

    for (i = 0; i < SIZ_TEXT; i++)                     /* copy text */
        zoo_header->text[i] = bytes[TEXT_I + i];

    zoo_header->zoo_tag = to_long(&bytes[ZTAG_I]);     /* copy zoo_tag */
    zoo_header->zoo_start = to_long(&bytes[ZST_I]);    /* copy zoo_start */
    zoo_header->zoo_minus = to_long(&bytes[ZSTM_I]);
}

/* ------------------------------------------------------------------------ */
/*  frd_dir.   Reads a directory entry in a machine-independent manner.     */
/*             Return 0 if OK, -1 if error.                                 */

int frd_dir(struct dirent *direntry)
{
    int status;
    unsigned char bytes[MAXDIRSIZE]; /* big enough to hold variable part too */

    /* To simplify things, we read the maximum possible size of the
       directory entry including the variable size and discard what is not
       needed */

    status = fread((char *)bytes,1,MAXDIRSIZE,ZOOFd);
    if (status < SIZ_DIR)
        return (-1);
    b_to_dir(direntry,bytes);
    return (0);
}

/* ------------------------------------------------------------------------ */
/*  b_to_dir.  Converts bytes to directory entry structure.                 */

void b_to_dir(struct dirent *direntry, unsigned char bytes[])
{
    int i;

    direntry->zoo_tag = to_long(&bytes[DTAG_I]);
    direntry->type = bytes[DTYP_I];
    direntry->packing_method = bytes[PKM_I];
    direntry->next = to_long(&bytes[NXT_I]);
    direntry->offset = to_long(&bytes[OFS_I]);
    direntry->org_size = to_long(&bytes[ORGS_I]);
    direntry->deleted = bytes[DEL_I];

    for (i = 0; i < FNM_SIZ; i++)
        direntry->fname[i] = bytes[FNAME_I + i];
}

/* ------------------------------------------------------------------------ */
/*  to_long. Converts four consecutive bytes, in order of increasing        */
/*           significance, to a long integer.  It is used to make Zoo       */
/*           independent of the byte order of the system.                   */

TLONG to_long(unsigned char data[])
{
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
            stack[stack_pointer++] = (x); \
            if (stack_pointer >= STACKSIZE)    \
                printfe("STACK OVERFLOW IN LZD().\n"); \
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
       printf1("\nassertion error \n"); \
   }
#else
#define assert(E) \
   { if (!(E)) \
       printf1(" \nassertion error \n"); \
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

int dbf_init(int opt)
{
    int err = -1;

    if (opt == 0) {
        err = 0;
        if (ZAAlloc)
            goto DBFin10;
        else
            goto DBFin;
    }

    ZAInit = ZAEOF = Out_Buf_Cnt = 0;
    if (ZAAlloc) { 
        In_Buf_Adr = Out_Buf_Adr + OUT_BUF_SIZE + (FUDGE/2);
        return(0);
    }
    if (!(c_table = (unsigned short *) calloc(4096,sizeof(unsigned short))))
        goto DBFin;
    memrq(4096,sizeof(unsigned short));
   
    if (!(pt_table = (unsigned short *) calloc(256,sizeof(unsigned short))))
        goto DBFin1;
    memrq(256,sizeof(unsigned short));
  
    if (!(c_len = (unsigned char *) calloc(NC,sizeof(unsigned char))))  
        goto DBFin2;
    memrq(NC,sizeof(char));
  
    if (!(pt_len = (unsigned char *) calloc(NPT,sizeof(unsigned char))))  
        goto DBFin3;
    memrq(NPT,sizeof(unsigned char));
  
    if (!(left = (unsigned short *) calloc(2 * NC - 1,sizeof(unsigned short))))
        goto DBFin4;
    memrq(2 * NC - 1,sizeof(unsigned short));
  
    if (!(right = (unsigned short *) calloc(2 * NC - 1,sizeof(unsigned short))))
        goto DBFin5;
    memrq(2 * NC - 1,sizeof(unsigned short));
  
    if (!(Out_Buf_Adr = (char *) calloc(MEM_BLOCK_SIZE,sizeof(char))))  
        goto DBFin6;
    memrq(MEM_BLOCK_SIZE,sizeof(char));
   
    if (!(ZABuf = (char *) calloc(ZABLen + 4,sizeof(char))))
        goto DBFin7;
    memrq(ZABLen + 4,sizeof(char));

    if (!(stack = (unsigned int *) calloc(STACKSIZE + 20,sizeof(unsigned int))))  
        goto DBFin8;
    memrq(STACKSIZE + 20,sizeof(unsigned int));

    if (!(table = (struct tabentry *) calloc(MAXMAX + 10,sizeof(struct tabentry))))  
        goto DBFin9;
    memrq(MAXMAX + 10,sizeof(struct tabentry));

    In_Buf_Adr = Out_Buf_Adr + OUT_BUF_SIZE + (FUDGE/2);
    ZAAlloc = 1;
    return(0);

DBFin10:
    free((char *)table);
    memrq(-MAXMAX - 10,sizeof(struct tabentry));
DBFin9:
    free((char *)stack);
    memrq(-STACKSIZE - 20,sizeof(unsigned int));
DBFin8:
    free((char *)ZABuf);
    memrq(-ZABLen - 4,sizeof(char));
DBFin7:
    free((char *)Out_Buf_Adr);
    memrq(-MEM_BLOCK_SIZE,sizeof(char));
DBFin6:
    free((char *)right);
    memrq(-2 * NC + 1,sizeof(unsigned short));
DBFin5:
    free((char *)left);
    memrq(-2 * NC + 1,sizeof(unsigned short));
DBFin4:
    free((char *)pt_len);
    memrq(-NPT,sizeof(unsigned char));
DBFin3: 
    free((char *)c_len);
    memrq(-NC,sizeof(char));
DBFin2:
    free((char *)pt_table);
    memrq(-256,sizeof(unsigned short));
DBFin1:
    free((char *)c_table);
    memrq(-4096,sizeof(unsigned short));
DBFin:
    ZAAlloc = 0;
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

int lzd(int fn)
{
    Out_Buf_Cnt = 0;
    if (ZAEOF)
        return(0);

    if (!ZAInit) {
        FPtr = ZAZOfs[fn];
        if (fseek(ZOOFd,FPtr,0) == -1L) {
            printfe("SEEK ERROR IN LZD. CAN't CONTINUE.\n");
            gerr_exit(18);
        }
        nbits = 9;
        max_code = 512;
        free_code = FIRST_FREE;
        stack_pointer = 0;
        bit_offset = 0;
        if (fread(In_Buf_Adr,1,INBUFSIZ,ZOOFd) == -1) {
            printfe("READ ERROR IN LZD. CAN'T CONTINUE.\n");
            gerr_exit(19);
        }
        FPtr = ftell(ZOOFd);
        init_dtab();           /* initialize table */
        ZAInit = 1;
    }
loop:
    if (Out_Buf_Cnt >= LZD_COUNT)
        return(Out_Buf_Cnt);

    cur_code = rd_dcode();

goteof:         /* special case for CLEAR then Z_EOF, for 0-length files */

    if (cur_code == Z_EOF) {
        ZAEOF = 1;
        return(Out_Buf_Cnt);
    }
    assert(nbits >= 9 && nbits <= 13);
    if (cur_code == CLEAR) {
        init_dtab();
        fin_char = kk = old_code = cur_code = rd_dcode();

        if (cur_code == Z_EOF)    /* special case for 0-length files */
            goto goteof;
        wr_dchar(kk);
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
    kk = fin_char = cur_code;
    push(kk);
    while (stack_pointer != 0) {
        wr_dchar(pop());
    }
    assert(nbits >= 9 && nbits <= 13);
    ad_dcode();
    old_code = in_code;
    assert(nbits >= 9 && nbits <= 13);
    goto loop;
}

/* ------------------------------------------------------------------------ */
/*  init_dtab.  Initialize table.                                           */

void init_dtab(void)
{
    nbits = 9;
    max_code = 512;
    free_code = FIRST_FREE;
}

/* ------------------------------------------------------------------------ */
/*  wr_dchar                                                                */

void wr_dchar(int ch)
{
    if (Out_Buf_Cnt >= OUTBUFSIZ) {      /* if buffer full */
        printfe("LZD OUTPUT BUF OVERFLOW. CAN'T CONTINUE.\n");
        gerr_exit(20);
    }
    assert(Out_Buf_Cnt < OUTBUFSIZ);
    Out_Buf_Adr[Out_Buf_Cnt++] = ch;            /* store character */
}

/* ------------------------------------------------------------------------ */
/*  rd_dcode. Reads a code from the input (compressed) file and returns     */
/*            its value.                                                    */

unsigned int rd_dcode(void)
{
    register char *ptra, *ptrb;                 /* miscellaneous pointers */
    unsigned int word;                          /* first 16 bits in buffer */
    unsigned int byte_offset;
    char nextch;                                /* next 8 bits in buffer */
    unsigned int ofs_inbyte;                    /* offset within byte */
    ofs_inbyte = bit_offset % 8;
    byte_offset = bit_offset / 8;
    bit_offset = bit_offset + nbits;
    assert(nbits >= 9 && nbits <= 13);

    if (byte_offset >= INBUFSIZ - 5) {
        int space_left;

        assert(byte_offset >= INBUFSIZ - 5);
        bit_offset = ofs_inbyte + nbits;
        space_left = INBUFSIZ - byte_offset;
        ptrb = byte_offset + In_Buf_Adr;          /* point to char */
        ptra = In_Buf_Adr;

        /* we now move the remaining characters down buffer beginning */
        while (space_left > 0) {
            *ptra++ = *ptrb++;
            space_left--;
        }
        assert(ptra - In_Buf_Adr == ptrb - (In_Buf_Adr + byte_offset));
        assert(space_left == 0);
        if (fseek(ZOOFd,FPtr,0) == -1L) {
            printfe("SEEK ERROR IN LZD. CAN'T CONTINUE.\n");
            gerr_exit(21);
        }
        if (fread(ptra,1,byte_offset,ZOOFd) == -1) {
            printfe("READ ERROR IN LZD. CAN'T CONTINUE.\n");
            gerr_exit(22);
        }
        FPtr = ftell(ZOOFd);
        byte_offset = 0;
    }
    ptra = byte_offset + In_Buf_Adr;

    /* NOTE: "word = *((int *) ptra)" would not be independent of byte order. */

    word = (unsigned char) *ptra; ptra++;
    word = word | ( ((unsigned char) *ptra) << 8 ); ptra++;
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

void ad_dcode(void)
{
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

int lzh_decode(int fn)
{

    if (ZAEOF) {
        Out_Buf_Cnt = 0;
        return(0);
    }
    if (!ZAInit) {
        FPtr = ZAZOfs[fn];
        if (fseek(ZOOFd,FPtr,0) == -1L) {
            printfe("SEEK ERROR IN LZD. CAN'T CONTINUE.\n");
            gerr_exit(23);
        }
        decode_start();
        ZAInit = 1;
    }
    else {
        if (fseek(ZOOFd,FPtr,0) == -1L) {
            printfe("SEEK ERROR IN LZD. CAN'T CONTINUE.\n");
            gerr_exit(24);
        }
    }
    if (!decoded)  
        Out_Buf_Cnt = decode((unsigned int)DICSIZ,(unsigned char *)Out_Buf_Adr);

    FPtr = ftell(ZOOFd);
    if (decoded)
        ZAEOF = 1;
    return(Out_Buf_Cnt);
}

/* ------------------------------------------------------------------------ */
/*  fillbuf    Shift bitbuf n bits left, read n bits                        */

void fillbuf(int n)
{
    bitbuf <<= n;
    while (n > bitcount) {
        bitbuf |= subbitbuf << (n -= bitcount);
        if (feof(ZOOFd))
            subbitbuf = 0;
        else
            subbitbuf = (unsigned char)getc(ZOOFd);
        bitcount = CHAR_BITS;
    }
    bitbuf |= subbitbuf >> (bitcount -= n);
}

/* ------------------------------------------------------------------------ */
/*  getbits                                                                 */

unsigned int getbits(int n)
{
    unsigned int x;

    x = bitbuf >> (BITBUFSIZ - n); 
    fillbuf(n);
    return(x);
}

/* ------------------------------------------------------------------------ */
/*  init_getbits                                                            */

void init_getbits(void)
{
    bitbuf = 0; 
    subbitbuf = 0;
    bitcount = 0;
    fillbuf(BITBUFSIZ);
}

/* ------------------------------------------------------------------------ */
/*  maketable. Make table for decoding                                      */

void make_table(int nchar, unsigned char *bitlen, int tablebits,
     unsigned short *xtable)
{
    unsigned short count[17], weight[17], start[18], *p;
    unsigned int i, k, len, ch, jutbits, avail, nextcode, mask;

    for (i = 1; i <= 16; i++)
        count[i] = 0;
    for (i = 0; i < nchar; i++) 
        count[bitlen[i]]++;

    start[1] = 0;
    for (i = 1; i <= 16; i++)
        start[i + 1] = start[i] + (count[i] << (16 - i));
    if (start[17] != (unsigned short)((unsigned int) 1 << 16)) {
        printfe("LZH BAD DECODE TABLE. CAN'T CONTINUE.\n");
        gerr_exit(25);
    }
    jutbits = 16 - tablebits;
    for (i = 1; i <= tablebits; i++) {
        start[i] >>= jutbits;
        weight[i] = (unsigned int) 1 << (tablebits - i);
    }
    while (i <= 16) {
        weight[i] = (unsigned int) 1 << (16 - i);
        i++;
    }
    i = start[tablebits + 1] >> jutbits;
    if (i != (unsigned short)((unsigned int) 1 << 16)) {
        k = 1 << tablebits;
        while (i != k)
            xtable[i++] = 0;
    }
    avail = nchar;
    mask = (unsigned int) 1 << (15 - tablebits);
    for (ch = 0; ch < nchar; ch++) {
        if ((len = bitlen[ch]) == 0)
            continue;
        nextcode = start[len] + weight[len];
        if (len <= tablebits) {
            for (i = start[len]; i < nextcode; i++) xtable[i] = ch;
        }
        else {
            k = start[len];
            p = &xtable[k >> jutbits];
            i = len - tablebits;
            while (i != 0) {
                if (*p == 0) {
                    right[avail] = left[avail] = 0;
                    *p = avail++;
                }
                if (k & mask)
                    p = &right[*p];
                else    
                    p = &left[*p];
                k <<= 1;  i--;
            }
            *p = ch;
        }
        start[len] = nextcode;
    }
}

/* ------------------------------------------------------------------------ */
/*  read_pt_len                                                             */

static void read_pt_len(int nn, int nbit, int i_special)
{
    int i, c, n;
    unsigned int mask;

    n = getbits(nbit);
    if (n == 0) {
        c = getbits(nbit);
        for (i = 0; i < nn; i++) pt_len[i] = 0;
        for (i = 0; i < 256; i++) pt_table[i] = c;
    }
    else {
        i = 0;
        while (i < n) {
            c = bitbuf >> (BITBUFSIZ - 3);
            if (c == 7) {
                mask = (unsigned int) 1 << (BITBUFSIZ - 1 - 3);
                while (mask & bitbuf) {  mask >>= 1;  c++;  }
            }
            fillbuf((c < 7) ? 3 : c - 3);
            pt_len[i++] = c;
            if (i == i_special) {
                c = getbits(2);
                while (--c >= 0) pt_len[i++] = 0;
            }
        }
        while (i < nn) pt_len[i++] = 0;
        make_table(nn,pt_len,8,pt_table);
    }
}

/* ------------------------------------------------------------------------ */
/*  read_c_len                                                              */

static void read_c_len(void)
{
    int i, c, n;
    unsigned int mask;
    n = getbits(CBIT);
    if (n == 0) {
        c = getbits(CBIT);
        for (i = 0; i < NC; i++) c_len[i] = 0;
        for (i = 0; i < 4096; i++) c_table[i] = c;
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
            fillbuf((int) pt_len[c]);
            if (c <= 2) {
                if      (c == 0) c = 1;
                else if (c == 1) c = getbits(4) + 3;
                else             c = getbits(CBIT) + 20;
                while (--c >= 0) c_len[i++] = 0;
            } else c_len[i++] = c - 2;
        }
        while (i < NC) c_len[i++] = 0;
        make_table(NC,c_len,12,c_table);
    }
}

/* ------------------------------------------------------------------------ */
/*  decode_c                                                                */

unsigned int decode_c(void)
{
    unsigned int j, mask;
    if (blocksize == 0) {
        blocksize = getbits(16);
        if (blocksize == 0) {
            decoded = 1;
            return 0;
        }
        read_pt_len(NT,TBIT,3);
        read_c_len();
        read_pt_len(NP,PBIT,-1);
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
    fillbuf((int) c_len[j]);
    return(j);
}

/* ------------------------------------------------------------------------ */
/*  decode_p                                                                */

unsigned int decode_p(void)
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
    fillbuf((int)pt_len[j]);
    if (j != 0) j = ((unsigned int) 1 << (j - 1)) + getbits((int) (j - 1));
    return(j);
}

/* ------------------------------------------------------------------------ */
/*  huf_decode_start                                                        */

void huf_decode_start(void)
{
    init_getbits(); 
    blocksize = 0;
}

/* ------------------------------------------------------------------------ */
/*  decode_start.                                                           */

void decode_start(void)
{
    huf_decode_start();
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

int decode(unsigned int count, unsigned char buffer[])
{
    static unsigned int i;
    unsigned int r, c;
    r = 0;  

    while (--jj >= 0) {
        buffer[r] = buffer[i];
        i = (i + 1) & (DICSIZ - 1);
        if (++r == count)
            return r;
    }
    for ( ; ; ) {
        c = decode_c();
        if (decoded)  
            return(r);
        if (c <= UCHAR_MX) {
            buffer[r] = c;
            if (++r == count)    
                return(r);
        }
        else {
            jj = c - (UCHAR_MX + 1 - THRESHOLD);
            i = (r - decode_p() - 1) & (DICSIZ - 1);
            while (--jj >= 0) {
                buffer[r] = buffer[i];
                i = (i + 1) & (DICSIZ - 1);
                if (++r == count)     
                    return(r);
            }
        }
    }
}


