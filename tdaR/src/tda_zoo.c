/****************************************************************************/
/*  tda_zoo.c                                                              */
/*                                                                          */
/*  Read and write Zoo archives (.zoo) -- tdaR's own code, NOT part of     */
/*  TDA. TDA itself only ever read archives, and only through commands     */
/*  (arcd/arcv/...) that require a companion .zad description file and     */
/*  load data straight into TDA's own data matrix; it has no "extract to   */
/*  a file" or "create an archive" capability at all. This gives R that:   */
/*  zoo() and unzoo(), mimicking base R's zip()/unzip().                   */
/*                                                                          */
/*  The container format (header/directory layout) and both packing        */
/*  methods -- 1 (LZD) and 2 (LZH) -- below are a from-scratch,            */
/*  self-contained port, cross-checked against t_zoo.c's own decoders      */
/*  (the algorithm TDA itself has used since the 1990s) and against the    */
/*  reference zoo/unzoo tools (https://github.com/troglobit/zoo).          */
/*  Archives this writes use method 1 (LZD) unless the caller asks for     */
/*  method 0 (stored); TDA's own arcd reads only methods 1 and 2.          */
/*                                                                          */
/*  Directory entries are written as zoo 2.1's type 2: the 56-byte fixed   */
/*  part is followed by a variable part carrying a long file name (when    */
/*  the name is not already a DOS-style 8.3 name) and a directory name,    */
/*  exactly as zoo itself stores them, so any zoo reader shows the same    */
/*  paths.  On read, both entry types are understood.                      */
/*                                                                          */
/*  Field offsets are the canonical Zoo 2.x layout (Rahul Dhesi's zoo,     */
/*  see zoo.h in the repository above).                                     */
/****************************************************************************/

#include <R.h>
#include <Rinternals.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ---- little-endian byte <-> integer helpers (zoo files are always LE) -- */

static uint32_t rd_u32(const unsigned char *p) {
    return (uint32_t) p[0] | ((uint32_t) p[1] << 8) |
           ((uint32_t) p[2] << 16) | ((uint32_t) p[3] << 24);
}
static void wr_u32(unsigned char *p, uint32_t v) {
    p[0] = (unsigned char) v; p[1] = (unsigned char) (v >> 8);
    p[2] = (unsigned char) (v >> 16); p[3] = (unsigned char) (v >> 24);
}
static void wr_u16(unsigned char *p, uint16_t v) {
    p[0] = (unsigned char) v; p[1] = (unsigned char) (v >> 8);
}
static uint16_t rd_u16(const unsigned char *p) {
    return (uint16_t) (p[0] | (p[1] << 8));
}

/* ---- canonical Zoo 2.x field layout ------------------------------------ */

#define ZOOTAG     0xFDC4A7DCu
#define SIZ_ZOOH   42          /* archive header length */
#define SIZ_TEXT   20          /* "ZOO x.xx Archive." text field */
#define ZTAG_I     20
#define ZST_I      24
#define ZSTM_I     28
#define MAJV_I     32
#define MINV_I     33
#define HTYPE_I    34
#define HVDATA_I   41
#define FILE_LEADER "@)#("      /* zoo writes these 5 bytes before each member */
#define SIZ_FLDR   5

#define SIZ_DIR    51          /* minimal (type-1) directory entry length */
#define DTAG_I     0
#define DTYP_I     4
#define PKM_I      5
#define NXT_I      6
#define OFS_I      10
#define DAT_I      14
#define TIM_I      16
#define CRC_I      18
#define ORGS_I     20
#define SIZNOW_I   24
#define DMAJ_I     28
#define DMIN_I     29
#define DEL_I      30
#define STRUC_I    31
#define CMT_I      32
#define CMTSIZ_I   36
#define FNAME_I    38
#define FNM_SIZ    13          /* zoo's short-name field: 12 chars + NUL */
#define VARDIRLEN_I 51         /* type 2: length of the variable part */
#define TZ_I       53
#define DCRC_I     54          /* type 2: CRC-16 of the whole entry */
#define SIZ_DIRL   56          /* type 2 fixed part */
#define NAMLEN_I   (SIZ_DIRL + 0)
#define DIRLEN_I   (SIZ_DIRL + 1)
#define LFNAME_I   (SIZ_DIRL + 2)
#define LFNAMESIZE 256
#define PATHSIZE   256
#define NO_TZ      127
#define VFL_ON     0x80
#define VFL_LAST   0x40

/* ------------------------------------------------------------------------ */
/*  LZD decoder (packing method 1): a self-contained port of lzd() in      */
/*  t_zoo.c, operating on a full in-memory buffer rather than TDA's        */
/*  sliding file-backed window (which exists only to bound TDA's own       */
/*  memory use; not needed when the whole entry is already in RAM).        */
/*  Verified byte-for-byte against both t_zoo.c's own output and           */
/*  tools/unzoo.py's independent reference decoder.                        */
/* ------------------------------------------------------------------------ */

enum { LZD_CLEAR = 256, LZD_EOF = 257, LZD_FIRST_FREE = 258,
       LZD_MAXBITS = 13, LZD_MAXMAX = 8192 };

static unsigned int lzd_read_code(const unsigned char *raw, long rawlen,
                                  long bit_offset, unsigned int nbits) {
    static const unsigned int masks[] =
        {0,0,0,0,0,0,0,0,0,0x1ff,0x3ff,0x7ff,0xfff,0x1fff};
    long byte_off = bit_offset / 8;
    unsigned int ofs_inbyte = (unsigned int) (bit_offset % 8);
    unsigned int word = (unsigned int) raw[byte_off];
    word |= ((unsigned int) (byte_off + 1 < rawlen ? raw[byte_off + 1] : 0)) << 8;
    unsigned int nextch = (byte_off + 2 < rawlen ? raw[byte_off + 2] : 0);
    if (ofs_inbyte)
        word = (word >> ofs_inbyte) | (nextch << (16 - ofs_inbyte));
    return word & masks[nbits];
}

/* Returns a malloc'd buffer of *outlen bytes (<= orig_size), or NULL on
   allocation failure. A truncated/corrupt stream simply yields fewer bytes
   than orig_size rather than reading out of bounds. */
static unsigned char *lzd_decode(const unsigned char *raw, long rawlen,
                                 long orig_size, long *outlen) {
    unsigned int *nextc = (unsigned int *) malloc(LZD_MAXMAX * sizeof(unsigned int));
    unsigned char *zch  = (unsigned char *) malloc(LZD_MAXMAX);
    unsigned char *out  = (unsigned char *) malloc((size_t) orig_size + 32);
    unsigned char *stackbuf = (unsigned char *) malloc(LZD_MAXMAX);
    if (!nextc || !zch || !out || !stackbuf) {
        free(nextc); free(zch); free(out); free(stackbuf);
        *outlen = 0;
        return NULL;
    }

    long opos = 0, bit_offset = 0, total_bits = rawlen * 8;
    unsigned int free_code = LZD_FIRST_FREE, nbits = 9, max_code = 512;
    unsigned int old_code = 0;
    unsigned char fin_char = 0, kk = 0;
    int have_old = 0;

    for (;;) {
        if (bit_offset + (long) nbits > total_bits)
            break;
        unsigned int cur_code = lzd_read_code(raw, rawlen, bit_offset, nbits);
        bit_offset += nbits;

        if (cur_code == LZD_EOF)
            break;

        if (cur_code == LZD_CLEAR) {
            free_code = LZD_FIRST_FREE; nbits = 9; max_code = 512;
            if (bit_offset + (long) nbits > total_bits)
                break;
            cur_code = lzd_read_code(raw, rawlen, bit_offset, nbits);
            bit_offset += nbits;
            if (cur_code == LZD_EOF)
                break;
            fin_char = kk = (unsigned char) cur_code;
            old_code = cur_code;
            have_old = 1;
            if (opos < orig_size + 16) out[opos++] = kk;
            continue;
        }

        unsigned int in_code = cur_code;
        int sp = 0;
        if (cur_code >= free_code) {
            if (!have_old)
                break;                       /* corrupt stream */
            cur_code = old_code;
            stackbuf[sp++] = fin_char;
        }
        while (cur_code > 255) {
            stackbuf[sp++] = zch[cur_code];
            cur_code = nextc[cur_code];
        }
        kk = fin_char = (unsigned char) cur_code;
        stackbuf[sp++] = kk;
        while (sp > 0) {
            if (opos >= orig_size + 16) break;
            out[opos++] = stackbuf[--sp];
        }

        if (have_old && free_code < LZD_MAXMAX) {
            zch[free_code] = kk;
            nextc[free_code] = old_code;
            free_code++;
            if (free_code >= max_code && nbits < LZD_MAXBITS) {
                nbits++;
                max_code <<= 1;
            }
        }
        old_code = in_code;
        have_old = 1;
    }

    free(nextc); free(zch); free(stackbuf);
    *outlen = opos;
    return out;
}

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/*  LZD ENCODER (packing method 1).                                         */
/*                                                                          */
/*  Port of lzc() from zoo's own lzc.c -- "Lempel-Ziv compression, mostly   */
/*  based on Tom Pfau's assembly language code", Rahul Dhesi 1986/12/31,    */
/*  released by him to the public domain.  The algorithm is unchanged; the  */
/*  only adaptation is I/O.  lzc.c reaches the outside world through two    */
/*  macros over an opaque handle (BLOCKREAD/BLOCKWRITE over BLOCKFILE), so  */
/*  reading from and writing to memory here needs no change to the          */
/*  compression itself.                                                     */
/*                                                                          */
/*  Why it is needed: TDA's arcd rejects any member whose packing method is */
/*  outside 1..2, so an archive of STORED members -- which every other zoo  */
/*  reader accepts -- is exactly the one TDA will not open.  Method 1 is    */
/*  the original and the one t_zoo.c's own decoder was written for.         */
/*                                                                          */
/*  The ratio check of the original is deliberately left out: it emits a    */
/*  CLEAR when compression is going badly, which is legal but makes output  */
/*  depend on a heuristic.  Omitting it costs a little compression on       */
/*  incompressible input and makes the encoder deterministic.               */

typedef struct { int first, next; unsigned char z_ch; } lzc_ent;

typedef struct {
    unsigned char *out;
    long           cap, len;          /* bytes written */
    unsigned long  bit_offset;
    int            nbits;
    unsigned int   free_code, max_code;
    lzc_ent       *tab;
} lzc_state;

static int lzc_grow(lzc_state *st, long need) {
    if (st->len + need <= st->cap) return 0;
    long cap = st->cap ? st->cap * 2 : 4096;
    while (cap < st->len + need) cap *= 2;
    unsigned char *q = (unsigned char *) realloc(st->out, (size_t) cap);
    if (!q) return -1;
    /* the bit writer ORs into the current byte, so new space starts zeroed */
    memset(q + st->cap, 0, (size_t) (cap - st->cap));
    st->out = q; st->cap = cap;
    return 0;
}

static int lzc_wr_code(lzc_state *st, int code) {
    long byte_offset = (long) (st->bit_offset / 8);
    unsigned int ofs = (unsigned int) (st->bit_offset % 8);
    if (lzc_grow(st, byte_offset + 4 - st->len)) return -1;
    st->bit_offset += (unsigned long) st->nbits;
    code &= 0xffff;
    if (ofs == 0) st->out[byte_offset]  = (unsigned char) (code & 0xff);
    else          st->out[byte_offset] |= (unsigned char) ((code << ofs) & 0xff);
    {
        unsigned int hib = ((unsigned int) code) >> (8 - ofs);
        st->out[byte_offset + 1] = (unsigned char) (hib & 0xff);
        st->out[byte_offset + 2] = (unsigned char) ((hib >> 8) & 0xff);
    }
    {
        long used = (long) ((st->bit_offset + 7) / 8);
        if (used > st->len) st->len = used;
    }
    return 0;
}

static void lzc_init_tab(lzc_state *st) {
    int i;
    for (i = 0; i < LZD_MAXMAX + 1; ++i) {
        st->tab[i].z_ch = 0;
        st->tab[i].first = st->tab[i].next = -1;
    }
    st->free_code = LZD_FIRST_FREE;
    st->nbits = 9;
    st->max_code = 512;
}

/* returns FOUND(0) / FIRST_USE(2) / NEXT_USE(1), *where = last entry seen */
static int lzc_lookup(lzc_state *st, int index, int ch, int *where) {
    *where = index;
    index = st->tab[index].first;
    if (index == -1) return 2;
    for (;;) {
        if ((st->tab[index].z_ch & 0xff) == (ch & 0xff)) { *where = index; return 0; }
        *where = index;
        index = st->tab[index].next;
        if (index == -1) return 1;
    }
}

static void lzc_add(lzc_state *st, int status, int ch, int index) {
    if (status == 1) st->tab[index].next  = (st->free_code >= LZD_MAXMAX ? -1 : (int) st->free_code);
    else             st->tab[index].first = (st->free_code >= LZD_MAXMAX ? -1 : (int) st->free_code);
    if (st->free_code <= LZD_MAXMAX) {
        st->tab[st->free_code].first = st->tab[st->free_code].next = -1;
        st->tab[st->free_code].z_ch = (unsigned char) (ch & 0xff);
        st->free_code++;
    }
}

/*  Compress in[0..inlen) into a freshly allocated buffer.  Returns the
    buffer and sets *outlen, or NULL on allocation failure.  */
static unsigned char *lzd_encode(const unsigned char *in, long inlen,
                                 long *outlen) {
    lzc_state st;
    long pos = 0;
    int nextch, prefix_code, k, status, where;

    memset(&st, 0, sizeof st);
    st.tab = (lzc_ent *) malloc((LZD_MAXMAX + 10) * sizeof(lzc_ent));
    if (!st.tab) return NULL;
    lzc_init_tab(&st);

    if (lzc_wr_code(&st, LZD_CLEAR)) goto fail;

    if (pos >= inlen) {                      /* empty member */
        if (lzc_wr_code(&st, LZD_EOF)) goto fail;
        free(st.tab); *outlen = st.len; return st.out;
    }
    nextch = in[pos++];

    for (;;) {
        nextch &= 0xff;
        for (;;) {
            prefix_code = nextch;
            if (pos >= inlen) {
                if (lzc_wr_code(&st, prefix_code)) goto fail;
                if (lzc_wr_code(&st, LZD_EOF)) goto fail;
                free(st.tab); *outlen = st.len; return st.out;
            }
            nextch = in[pos++] & 0xff;
            k = nextch;
            status = lzc_lookup(&st, prefix_code, nextch, &where);
            if (status != 0) break;
            nextch = where;                  /* FOUND: extend the match */
        }
        if (lzc_wr_code(&st, prefix_code)) goto fail;
        lzc_add(&st, status, k, where);
        nextch = k;
        if (st.free_code > st.max_code) {
            if (st.nbits >= LZD_MAXBITS) {
                if (lzc_wr_code(&st, LZD_CLEAR)) goto fail;
                lzc_init_tab(&st);
                continue;
            }
            st.nbits++;
            st.max_code <<= 1;
        }
    }
fail:
    free(st.tab); free(st.out);
    return NULL;
}

/*  LZH decoder (packing method 2): a self-contained port of lzh_decode()  */
/*  and its supporting functions in t_zoo.c -- a -lh5--like Huffman+LZSS   */
/*  variant with a 13-bit dictionary, unique to Zoo (not read by generic   */
/*  LHA/LZH libraries). Verified against a real archive: every compressed  */
/*  entry in a 67-file, real-world TDA .zoo archive decoded to exactly its */
/*  recorded original size.                                                */
/* ------------------------------------------------------------------------ */

enum {
    LZH_UCHAR_MX = 255,
    LZH_MAXMATCH = 256,
    LZH_THRESHOLD = 3,
    LZH_CBIT = 9,
    LZH_CODE_BIT = 16,
    LZH_DICBIT = 13,
    LZH_DICSIZ = 1 << LZH_DICBIT,
    LZH_NC = LZH_UCHAR_MX + LZH_MAXMATCH + 2 - LZH_THRESHOLD,  /* 510 */
    LZH_NP = LZH_DICBIT + 1,                                   /* 14 */
    LZH_NT = LZH_CODE_BIT + 3,                                 /* 19 */
    LZH_PBIT = 4,
    LZH_TBIT = 5,
    LZH_NPT = 19,                                              /* max(NT,NP) */
    LZH_BITBUFSIZ = 8 * 2                                      /* CHAR_BITS * sizeof(unsigned short) */
};

typedef struct {
    const unsigned char *in;
    long inlen, inpos;

    unsigned short bitbuf;
    unsigned int subbitbuf;
    int bitcount;

    unsigned short left[2 * LZH_NC - 1];
    unsigned short right[2 * LZH_NC - 1];
    unsigned char c_len[LZH_NC];
    unsigned char pt_len[LZH_NPT];
    unsigned int blocksize;
    unsigned short c_table[4096];
    unsigned short pt_table[256];

    int jj;
    char decoded;
    unsigned int s_decode_i;

    unsigned char dict[LZH_DICSIZ];  /* the ring-buffer / sliding-window */
} lzh_state;

static int lzh_getbyte(lzh_state *s) {
    return s->inpos < s->inlen ? s->in[s->inpos++] : -1;
}

static void lzh_fillbuf(lzh_state *s, int n) {
    s->bitbuf = (unsigned short) (s->bitbuf << n);
    while (n > s->bitcount) {
        s->bitbuf = (unsigned short) (s->bitbuf | (s->subbitbuf << (n -= s->bitcount)));
        int c = lzh_getbyte(s);
        s->subbitbuf = (c < 0) ? 0 : (unsigned char) c;
        s->bitcount = 8;
    }
    s->bitbuf = (unsigned short) (s->bitbuf | (s->subbitbuf >> (s->bitcount -= n)));
}

static unsigned int lzh_getbits(lzh_state *s, int n) {
    unsigned int x = (unsigned int) (s->bitbuf >> (LZH_BITBUFSIZ - n));
    lzh_fillbuf(s, n);
    return x;
}

static void lzh_init_getbits(lzh_state *s) {
    s->bitbuf = 0; s->subbitbuf = 0; s->bitcount = 0;
    lzh_fillbuf(s, LZH_BITBUFSIZ);
}

static void lzh_make_table(lzh_state *s, int nchar, const unsigned char *bitlen,
                           int tablebits, unsigned short *xtable) {
    unsigned short count[17], weight[17], start[18], *p;
    unsigned int i, k, len, ch, jutbits, avail, nextcode, mask;

    for (i = 1; i <= 16; i++) count[i] = 0;
    for (i = 0; i < (unsigned) nchar; i++) count[bitlen[i]]++;

    start[1] = 0;
    for (i = 1; i <= 16; i++)
        start[i + 1] = (unsigned short) (start[i] + (count[i] << (16 - i)));

    jutbits = 16 - (unsigned) tablebits;
    for (i = 1; i <= (unsigned) tablebits; i++) {
        start[i] = (unsigned short) (start[i] >> jutbits);
        weight[i] = (unsigned short)((unsigned int) 1 << ((unsigned int)tablebits - i));
    }
    while (i <= 16) { weight[i] = (unsigned short)(1U << (16U - (unsigned int)i)); i++; }

    i = start[tablebits + 1] >> jutbits;
    if (i != (unsigned short) ((unsigned int) 1 << 16)) {
        k = (unsigned int) 1 << tablebits;
        while (i != k) xtable[i++] = 0;
    }
    avail = (unsigned) nchar;
    mask = (unsigned int) 1 << (15 - tablebits);
    for (ch = 0; ch < (unsigned) nchar; ch++) {
        if ((len = bitlen[ch]) == 0) continue;
        nextcode = start[len] + weight[len];
        if (len <= (unsigned) tablebits) {
            for (i = start[len]; i < nextcode; i++) xtable[i] = (unsigned short) ch;
        } else {
            k = start[len];
            p = &xtable[k >> jutbits];
            i = len - (unsigned) tablebits;
            while (i != 0) {
                if (*p == 0) {
                    s->right[avail] = s->left[avail] = 0;
                    *p = (unsigned short) avail++;
                }
                p = (k & mask) ? &s->right[*p] : &s->left[*p];
                k <<= 1; i--;
            }
            *p = (unsigned short) ch;
        }
        start[len] = (unsigned short) nextcode;
    }
}

static void lzh_read_pt_len(lzh_state *s, int nn, int nbit, int i_special) {
    int i, c, n;
    unsigned int mask;
    n = (int) lzh_getbits(s, nbit);
    if (n == 0) {
        c = (int) lzh_getbits(s, nbit);
        for (i = 0; i < nn; i++) s->pt_len[i] = 0;
        for (i = 0; i < 256; i++) s->pt_table[i] = (unsigned short) c;
    } else {
        i = 0;
        while (i < n) {
            c = s->bitbuf >> (LZH_BITBUFSIZ - 3);
            if (c == 7) {
                mask = (unsigned int) 1 << (LZH_BITBUFSIZ - 1 - 3);
                while (mask & s->bitbuf) { mask >>= 1; c++; }
            }
            lzh_fillbuf(s, (c < 7) ? 3 : c - 3);
            s->pt_len[i++] = (unsigned char) c;
            if (i == i_special) {
                c = (int) lzh_getbits(s, 2);
                while (--c >= 0) s->pt_len[i++] = 0;
            }
        }
        while (i < nn) s->pt_len[i++] = 0;
        lzh_make_table(s, nn, s->pt_len, 8, s->pt_table);
    }
}

static void lzh_read_c_len(lzh_state *s) {
    int i, c, n;
    unsigned int mask;
    n = (int) lzh_getbits(s, LZH_CBIT);
    if (n == 0) {
        c = (int) lzh_getbits(s, LZH_CBIT);
        for (i = 0; i < LZH_NC; i++) s->c_len[i] = 0;
        for (i = 0; i < 4096; i++) s->c_table[i] = (unsigned short) c;
    } else {
        i = 0;
        while (i < n) {
            c = s->pt_table[s->bitbuf >> (LZH_BITBUFSIZ - 8)];
            if (c >= LZH_NT) {
                mask = (unsigned int) 1 << (LZH_BITBUFSIZ - 1 - 8);
                do {
                    c = (s->bitbuf & mask) ? s->right[c] : s->left[c];
                    mask >>= 1;
                } while (c >= LZH_NT);
            }
            lzh_fillbuf(s, s->pt_len[c]);
            if (c <= 2) {
                if (c == 0) c = 1;
                else if (c == 1) c = (int) lzh_getbits(s, 4) + 3;
                else c = (int) lzh_getbits(s, LZH_CBIT) + 20;
                while (--c >= 0) s->c_len[i++] = 0;
            } else {
                s->c_len[i++] = (unsigned char) (c - 2);
            }
        }
        while (i < LZH_NC) s->c_len[i++] = 0;
        lzh_make_table(s, LZH_NC, s->c_len, 12, s->c_table);
    }
}

static unsigned int lzh_decode_c(lzh_state *s) {
    unsigned int j, mask;
    if (s->blocksize == 0) {
        s->blocksize = lzh_getbits(s, 16);
        if (s->blocksize == 0) { s->decoded = 1; return 0; }
        lzh_read_pt_len(s, LZH_NT, LZH_TBIT, 3);
        lzh_read_c_len(s);
        lzh_read_pt_len(s, LZH_NP, LZH_PBIT, -1);
    }
    s->blocksize--;
    j = s->c_table[s->bitbuf >> (LZH_BITBUFSIZ - 12)];
    if (j >= LZH_NC) {
        mask = (unsigned int) 1 << (LZH_BITBUFSIZ - 1 - 12);
        do {
            j = (s->bitbuf & mask) ? s->right[j] : s->left[j];
            mask >>= 1;
        } while (j >= LZH_NC);
    }
    lzh_fillbuf(s, s->c_len[j]);
    return j;
}

static unsigned int lzh_decode_p(lzh_state *s) {
    unsigned int j, mask;
    j = s->pt_table[s->bitbuf >> (LZH_BITBUFSIZ - 8)];
    if (j >= LZH_NP) {
        mask = (unsigned int) 1 << (LZH_BITBUFSIZ - 1 - 8);
        do {
            j = (s->bitbuf & mask) ? s->right[j] : s->left[j];
            mask >>= 1;
        } while (j >= LZH_NP);
    }
    lzh_fillbuf(s, s->pt_len[j]);
    if (j != 0) j = ((unsigned int) 1 << (j - 1)) + lzh_getbits(s, (int) (j - 1));
    return j;
}

static unsigned int lzh_decode_chunk(lzh_state *s, unsigned int count,
                                     unsigned char *buffer) {
    unsigned int r = 0, c;
    while (--s->jj >= 0) {
        buffer[r] = buffer[s->s_decode_i];
        s->s_decode_i = (s->s_decode_i + 1) & (LZH_DICSIZ - 1);
        if (++r == count) return r;
    }
    for (;;) {
        c = lzh_decode_c(s);
        if (s->decoded) return r;
        if (c <= LZH_UCHAR_MX) {
            buffer[r] = (unsigned char) c;
            if (++r == count) return r;
        } else {
            s->jj = (int) (c - (LZH_UCHAR_MX + 1 - LZH_THRESHOLD));
            s->s_decode_i = (r - lzh_decode_p(s) - 1) & (LZH_DICSIZ - 1);
            while (--s->jj >= 0) {
                buffer[r] = buffer[s->s_decode_i];
                s->s_decode_i = (s->s_decode_i + 1) & (LZH_DICSIZ - 1);
                if (++r == count) return r;
            }
        }
    }
}

/* Driver: repeatedly decodes DICSIZ-size chunks, exactly as lzh_decode()'s
   caller does in t_zoo.c, copying each chunk to the full-size output since
   the dictionary buffer itself is reused/overwritten between chunks.
   Returns a malloc'd buffer of *outlen bytes (<= orig_size). */
static unsigned char *lzh_decode_all(const unsigned char *raw, long rawlen,
                                     long orig_size, long *outlen) {
    lzh_state *s = (lzh_state *) calloc(1, sizeof(lzh_state));
    unsigned char *out = (unsigned char *) malloc((size_t) orig_size + 16);
    if (!s || !out) { free(s); free(out); *outlen = 0; return NULL; }

    s->in = raw; s->inlen = rawlen; s->inpos = 0;
    lzh_init_getbits(s);
    s->blocksize = 0; s->jj = 0; s->decoded = 0;

    long opos = 0;
    while (!s->decoded && opos < orig_size) {
        unsigned int got = lzh_decode_chunk(s, LZH_DICSIZ, s->dict);
        if (got == 0) break;
        long take = got;
        if (opos + take > orig_size) take = orig_size - opos;
        memcpy(out + opos, s->dict, (size_t) take);
        opos += take;
        if (got < LZH_DICSIZ) break;
    }
    *outlen = opos;
    free(s);
    return out;
}

/* ------------------------------------------------------------------------ */
/*  CRC-16/ARC, as addbfcrc() in zoo's own crcdefs.c/addbfcrc.c computes    */
/*  it: poly 0xA001 (reflected), init 0, no final XOR, reset to 0 for each  */
/*  file. Writing a directory entry with this left as 0 (as an earlier     */
/*  version of this file did) produces a file every real zoo reader can    */
/*  still extract correctly, but flags as "CRC failed" -- content-correct  */
/*  but not what a conforming archive looks like, so it is computed here.  */
/*  Bit-by-bit rather than a hand-transcribed lookup table, to keep a      */
/*  256-entry copy-paste error off the table entirely; verified against    */
/*  the standard CRC-16/ARC test vector ("123456789" -> 0xBB3D).           */
/* ------------------------------------------------------------------------ */

static unsigned short zoo_crc16(const unsigned char *buf, long count) {
    unsigned int crc = 0;
    while (count-- > 0) {
        crc ^= *buf++;
        for (int b = 0; b < 8; b++)
            crc = (crc & 1) ? (crc >> 1) ^ 0xA001 : (crc >> 1);
    }
    return (unsigned short) crc;
}

/* ------------------------------------------------------------------------ */
/*  Reading an archive's directory                                         */
/* ------------------------------------------------------------------------ */

typedef struct {
    char name[PATHSIZE + LFNAMESIZE + 1];
    unsigned char method;
    long offset, orig_size, packed_size;
    int deleted;
} zoo_entry;

/* Parses the directory into *out (caller frees), returns the entry count,
   or -1 if this is not a Zoo archive at all. */
static int zoo_read_dir(const unsigned char *data, long len, zoo_entry **out) {
    if (len < SIZ_ZOOH || rd_u32(data + ZTAG_I) != ZOOTAG)
        return -1;
    int32_t start = (int32_t) rd_u32(data + ZST_I);
    uint32_t minus = rd_u32(data + ZSTM_I);
    if ((uint32_t) start + minus != 0)
        return -1;

    int cap = 16, n = 0;
    zoo_entry *entries = (zoo_entry *) malloc((size_t) cap * sizeof(zoo_entry));
    long off = start;
    while (off > 0 && off + SIZ_DIR <= len) {
        if (rd_u32(data + off + DTAG_I) != ZOOTAG) {
            free(entries);
            return -1;
        }
        long nxt = (int32_t) rd_u32(data + off + NXT_I);
        /* The record with next==0 is a terminal end-marker, not a real
           file (TDA's own get_zoo() and the independent tools/unzoo.py
           reference both stop before adding it) -- adding it here as a
           68th, empty, zero-size "file" was a real bug. */
        if (nxt == 0)
            break;
        if (n >= cap) {
            /* realloc returns NULL on failure without freeing the old
               block: assigning it straight back both leaked the entries
               read so far and left the &entries[n] below dereferencing
               NULL (cppcheck) */
            zoo_entry *grown;
            cap *= 2;
            grown = (zoo_entry *) realloc(entries,
                                          (size_t) cap * sizeof(zoo_entry));
            if (grown == NULL) {
                free(entries);
                return -1;
            }
            entries = grown;
        }
        zoo_entry *e = &entries[n];
        e->method = data[off + PKM_I];
        e->offset = (int32_t) rd_u32(data + off + OFS_I);
        e->orig_size = (int32_t) rd_u32(data + off + ORGS_I);
        e->packed_size = (int32_t) rd_u32(data + off + SIZNOW_I);
        e->deleted = data[off + DEL_I];
        /* The name: for a type-2 entry the long name from the variable
           part when there is one, prefixed by the stored directory;
           otherwise the 13-byte short field. Both lengths count the NUL. */
        char fname[FNM_SIZ];
        memcpy(fname, data + off + FNAME_I, FNM_SIZ);
        fname[FNM_SIZ - 1] = '\0';
        const char *lname = fname;
        const char *dname = "";
        char lbuf[LFNAMESIZE], dbuf[PATHSIZE];
        if (data[off + DTYP_I] == 2 && off + SIZ_DIRL <= len) {
            int vdl = rd_u16(data + off + VARDIRLEN_I);
            int namlen = vdl > 0 && off + NAMLEN_I < len ? data[off + NAMLEN_I] : 0;
            int dirlen = vdl > 1 && off + DIRLEN_I < len ? data[off + DIRLEN_I] : 0;
            if (namlen + dirlen + 2 <= vdl &&
                off + LFNAME_I + namlen + dirlen <= len) {
                if (namlen > 1) {
                    memcpy(lbuf, data + off + LFNAME_I, (size_t) namlen);
                    lbuf[namlen - 1] = '\0';
                    lname = lbuf;
                }
                if (dirlen > 1) {
                    memcpy(dbuf, data + off + LFNAME_I + namlen, (size_t) dirlen);
                    dbuf[dirlen - 1] = '\0';
                    dname = dbuf;
                }
            }
        }
        if (*dname)
            snprintf(e->name, sizeof(e->name), "%s/%s", dname, lname);
        else
            snprintf(e->name, sizeof(e->name), "%s", lname);
        if (!e->deleted)
            n++;
        off = nxt;
    }
    *out = entries;
    return n;
}

/* ------------------------------------------------------------------------ */
/*  .Call entry points                                                     */
/* ------------------------------------------------------------------------ */

SEXP C_tda_unzoo(SEXP path_sexp) {
    const char *path = CHAR(STRING_ELT(path_sexp, 0));
    FILE *f = fopen(path, "rb");
    if (!f)
        Rf_error("cannot open '%s'", path);
    fseek(f, 0, SEEK_END);
    long flen = ftell(f);
    fseek(f, 0, SEEK_SET);
    unsigned char *data = (unsigned char *) R_alloc((size_t) flen, 1);
    if (fread(data, 1, (size_t) flen, f) != (size_t) flen) {
        fclose(f);
        Rf_error("short read on '%s'", path);
    }
    fclose(f);

    zoo_entry *entries;
    int n = zoo_read_dir(data, flen, &entries);
    if (n < 0)
        Rf_error("'%s' is not a Zoo archive (bad header)", path);

    SEXP names_r = PROTECT(Rf_allocVector(STRSXP, n));
    SEXP methods_r = PROTECT(Rf_allocVector(INTSXP, n));
    SEXP sizes_r = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP data_r = PROTECT(Rf_allocVector(VECSXP, n));
    SEXP ok_r = PROTECT(Rf_allocVector(LGLSXP, n));

    for (int i = 0; i < n; ++i) {
        zoo_entry *e = &entries[i];
        SET_STRING_ELT(names_r, i, Rf_mkChar(e->name));
        INTEGER(methods_r)[i] = e->method;
        REAL(sizes_r)[i] = (double) e->orig_size;

        if (e->method == 0) {
            SEXP raw = PROTECT(Rf_allocVector(RAWSXP, e->orig_size));
            if (e->orig_size > 0 && e->offset + e->orig_size <= flen)
                memcpy(RAW(raw), data + e->offset, (size_t) e->orig_size);
            SET_VECTOR_ELT(data_r, i, raw);
            LOGICAL(ok_r)[i] = 1;
            UNPROTECT(1);
        } else if (e->method == 1) {
            long outlen = 0;
            unsigned char *dec = lzd_decode(data + e->offset, flen - e->offset,
                                            e->orig_size, &outlen);
            SEXP raw = PROTECT(Rf_allocVector(RAWSXP, outlen));
            if (dec && outlen > 0)
                memcpy(RAW(raw), dec, (size_t) outlen);
            free(dec);
            SET_VECTOR_ELT(data_r, i, raw);
            LOGICAL(ok_r)[i] = (outlen == e->orig_size);
            UNPROTECT(1);
        } else if (e->method == 2) {
            long outlen = 0;
            unsigned char *dec = lzh_decode_all(data + e->offset, flen - e->offset,
                                                e->orig_size, &outlen);
            SEXP raw = PROTECT(Rf_allocVector(RAWSXP, outlen));
            if (dec && outlen > 0)
                memcpy(RAW(raw), dec, (size_t) outlen);
            free(dec);
            SET_VECTOR_ELT(data_r, i, raw);
            LOGICAL(ok_r)[i] = (outlen == e->orig_size);
            UNPROTECT(1);
        } else {
            /* An unrecognized packing method (only 0, 1, 2 are defined). */
            SET_VECTOR_ELT(data_r, i, Rf_allocVector(RAWSXP, 0));
            LOGICAL(ok_r)[i] = 0;
        }
    }
    free(entries);

    SEXP result = PROTECT(Rf_allocVector(VECSXP, 5));
    SET_VECTOR_ELT(result, 0, names_r);
    SET_VECTOR_ELT(result, 1, methods_r);
    SET_VECTOR_ELT(result, 2, sizes_r);
    SET_VECTOR_ELT(result, 3, data_r);
    SET_VECTOR_ELT(result, 4, ok_r);
    SEXP nm = PROTECT(Rf_allocVector(STRSXP, 5));
    SET_STRING_ELT(nm, 0, Rf_mkChar("name"));
    SET_STRING_ELT(nm, 1, Rf_mkChar("method"));
    SET_STRING_ELT(nm, 2, Rf_mkChar("size"));
    SET_STRING_ELT(nm, 3, Rf_mkChar("data"));
    SET_STRING_ELT(nm, 4, Rf_mkChar("ok"));
    Rf_setAttrib(result, R_NamesSymbol, nm);

    UNPROTECT(7);
    return result;
}

/* The short (DOS 8.3) form of a name, as zoo's dosname() derives it: up
   to eight characters of the root and three of the extension, every
   character outside zoo's legal set mapped as zoo maps it, an empty root
   becoming "X" and a leading dot an underscore.  This is the fname field
   every zoo reader shows when an entry has no long name, and what TDA's
   own get_zoo() matched before it learned the long name. */
static void zoo_dosname(const char *base, char *out) {
    static const char legal[] =
        "tabcdefghijklmnopqrs_uvwxyz0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ@^`{}~!#$%&'()-";
    char root[9], ext[4];
    const char *dot = strrchr(base, '.');
    size_t rl = dot ? (size_t) (dot - base) : strlen(base);
    if (rl > 8) rl = 8;
    memcpy(root, base, rl);
    root[rl] = '\0';
    ext[0] = '\0';
    if (dot && dot[1]) {
        strncpy(ext, dot + 1, 3);
        ext[3] = '\0';
    } else if (strlen(base) > 8) {   /* no extension: zoo takes characters
                                        9-11 of the whole name instead */
        strncpy(ext, base + 8, 3);
        ext[3] = '\0';
    }
    char *parts[2] = { root, ext };
    for (int k = 0; k < 2; ++k) {
        char *p = parts[k];
        if (k == 0 && *p == '\0')
            strcpy(p, "X");
        if (*p == '.')
            *p = '_';
        for (; *p; ++p)
            if (strchr(legal, *p) == NULL)
                *p = legal[(*p & 0xff) % 26];
    }
    if (ext[0])
        snprintf(out, FNM_SIZ, "%s.%s", root, ext);
    else
        snprintf(out, FNM_SIZ, "%s", root);
}

/* Writes a new archive at `path` containing `names` (paths relative to
   the archive root, '/' separated) with the raw bytes in `contents` (a
   list of raw vectors), each as a type-2 directory entry: the long name
   is stored when it differs from its 8.3 form, the directory part
   always. */
SEXP C_tda_zoo(SEXP path_sexp, SEXP names_sexp, SEXP contents_sexp,
               SEXP method_sexp) {
    int method = Rf_asInteger(method_sexp);
    const char *path = CHAR(STRING_ELT(path_sexp, 0));
    int n = (int) Rf_length(names_sexp);
    if (Rf_length(contents_sexp) != n)
        Rf_error("names and contents must have the same length");

    FILE *f = fopen(path, "wb");
    if (!f)
        Rf_error("cannot create '%s'", path);

    unsigned char hdr[SIZ_ZOOH];
    memset(hdr, 0, sizeof(hdr));
    /* the header zoo 2.1 writes for a new archive: its text ends in
       Ctrl-Z, version 2.0 is enough to manipulate it, header type 1,
       and the version data byte carries the default generation count */
    memcpy(hdr, "ZOO 2.10 Archive.\032", 18);
    wr_u32(hdr + ZTAG_I, ZOOTAG);
    wr_u32(hdr + ZST_I, SIZ_ZOOH);
    wr_u32(hdr + ZSTM_I, (uint32_t) (-(int32_t) SIZ_ZOOH));
    hdr[MAJV_I] = 2; hdr[MINV_I] = 0;
    hdr[HTYPE_I] = 1;
    hdr[HVDATA_I] = 3;
    fwrite(hdr, 1, SIZ_ZOOH, f);

    long *offsets = (long *)(void *) R_alloc((size_t) n, sizeof(long));
    long *sizes = (long *)(void *) R_alloc((size_t) n, sizeof(long));
    long pos = SIZ_ZOOH;

    for (int i = 0; i < n; ++i) {
        SEXP raw = VECTOR_ELT(contents_sexp, i);
        long len = Rf_length(raw);
        /*  Method 1 (LZD) by default.  TDA's arcd refuses anything outside
            1..2, so a stored archive -- which every other zoo reader takes
            -- is exactly the one TDA will not open.  If the compressed form
            is no smaller the member is stored instead, as zoo itself does. */
        const unsigned char *payload = RAW(raw);
        long paylen = len;
        unsigned char *packed = NULL;
        int pkm = 0;
        if (method == 1 && len > 0) {
            long plen = 0;
            packed = lzd_encode(RAW(raw), len, &plen);
            if (packed && plen < len) { payload = packed; paylen = plen; pkm = 1; }
            else { free(packed); packed = NULL; }
        }
        const char *nm = CHAR(STRING_ELT(names_sexp, i));
        const char *slash = strrchr(nm, '/');
        const char *base = slash ? slash + 1 : nm;
        size_t dlen = slash ? (size_t) (slash - nm) : 0;
        size_t blen = strlen(base);
        if (blen == 0 || blen >= LFNAMESIZE || dlen >= PATHSIZE)
            Rf_error("member name '%s' is empty or too long", nm);
        char fname[FNM_SIZ];
        zoo_dosname(base, fname);
        int namlen = strcmp(base, fname) ? (int) blen + 1 : 0;
        int dirlen = dlen ? (int) dlen + 1 : 0;
        /* 1 namlen, 1 dirlen, 2 system id, 3 attributes, 1 version flag,
           2 version number: zoo's newdir() */
        int vdl = namlen + dirlen + 10;
        int entlen = SIZ_DIRL + vdl;

        offsets[i] = pos + entlen + SIZ_FLDR;   /* the entry, the leader, the data */
        sizes[i] = paylen;

        unsigned char dir[SIZ_DIRL + 2 + LFNAMESIZE + PATHSIZE + 8];
        memset(dir, 0, sizeof(dir));
        wr_u32(dir + DTAG_I, ZOOTAG);
        dir[DTYP_I] = 2;
        dir[PKM_I] = (unsigned char) pkm;
        long next = pos + entlen + SIZ_FLDR + paylen;  /* always points onward: a real
            archive's last file still points to a genuine terminal dummy
            record (written below), never next==0 on a real file itself --
            that is what marks a record as *not* real data on read. */
        wr_u32(dir + NXT_I, (uint32_t) next);
        wr_u32(dir + OFS_I, (uint32_t) offsets[i]);
        wr_u32(dir + ORGS_I, (uint32_t) len);
        wr_u32(dir + SIZNOW_I, (uint32_t) paylen);  /* size in archive */
        /* the CRC is of the ORIGINAL bytes, not the packed ones */
        wr_u16(dir + CRC_I, zoo_crc16(RAW(raw), len));
        dir[DMAJ_I] = 1; dir[DMIN_I] = 0;    /* zoo 1.0 extracts LZD and stored */
        dir[DEL_I] = 0;
        dir[STRUC_I] = 0;
        memcpy(dir + FNAME_I, fname, strlen(fname));
        wr_u16(dir + VARDIRLEN_I, (uint16_t) vdl);
        dir[TZ_I] = NO_TZ;
        dir[NAMLEN_I] = (unsigned char) namlen;
        dir[DIRLEN_I] = (unsigned char) dirlen;
        int at = LFNAME_I;
        if (namlen) { memcpy(dir + at, base, (size_t) namlen); at += namlen; }
        if (dirlen) { memcpy(dir + at, nm, dlen); at += dirlen; }
        wr_u16(dir + at, 0);                 /* system id: unix */
        at += 2;
        at += 3;                             /* file attributes: none */
        dir[at] = VFL_ON | VFL_LAST;
        wr_u16(dir + at + 1, 1);             /* version number */
        /* the entry's own CRC, computed with the CRC field zero */
        wr_u16(dir + DCRC_I, zoo_crc16(dir, entlen));

        fwrite(dir, 1, (size_t) entlen, f);
        fwrite(FILE_LEADER, 1, SIZ_FLDR, f);
        if (paylen > 0)
            fwrite(payload, 1, (size_t) paylen, f);
        free(packed);
        pos = next;
    }
    /* Terminal dummy record: next==0 marks the end of the directory chain
       and is not itself a file (see zoo_read_dir() in this same file). */
    unsigned char term[SIZ_DIRL];
    memset(term, 0, sizeof(term));
    wr_u32(term + DTAG_I, ZOOTAG);
    term[DTYP_I] = 2;
    wr_u32(term + NXT_I, 0);
    wr_u16(term + DCRC_I, zoo_crc16(term, SIZ_DIRL));
    fwrite(term, 1, SIZ_DIRL, f);
    fclose(f);
    return R_NilValue;
}
