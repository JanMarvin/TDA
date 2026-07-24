/*  Standalone sanitiser harness for the LZD encoder/decoder pair in
    tda_zoo.c.  R's own build is not instrumented, so a round trip run
    from R proves correctness but says nothing about memory safety --
    this compiles the two functions on their own and runs them under
    ASAN, UBSan and leak detection.

    Regenerate after changing either function (the bodies are copied
    out of tda_zoo.c), then:

      gcc -std=gnu99 -g -O1 -fsanitize=address,undefined \
          -o /tmp/lzdh tools/lzd_harness.c && ASAN_OPTIONS=detect_leaks=1 /tmp/lzdh
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
enum { LZD_CLEAR = 256, LZD_EOF = 257, LZD_FIRST_FREE = 258,
       LZD_MAXBITS = 13, LZD_MAXMAX = 8192 };
/*  LZD decoder (packing method 1): a self-contained port of lzd() in      */
/*  t_zoo.c, operating on a full in-memory buffer rather than TDA's        */
/*  sliding file-backed window (which exists only to bound TDA's own       */
/*  memory use; not needed when the whole entry is already in RAM).        */
/*  Verified byte-for-byte against both t_zoo.c's own output and           */
/*  tools/unzoo.py's independent reference decoder.                        */
/* ------------------------------------------------------------------------ */


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


int main(void) {
    srand(7);
    for (int t = 0; t < 400; ++t) {
        long n = rand() % 20000;
        unsigned char *in = malloc(n ? n : 1);
        int mode = t % 3;
        for (long i = 0; i < n; ++i)
            in[i] = mode == 0 ? 0 : mode == 1 ? (unsigned char)(rand() & 0xff)
                                              : (unsigned char)(i % 7);
        long clen = 0;
        unsigned char *c = lzd_encode(in, n, &clen);
        if (!c) { printf("encode fail n=%ld\n", n); return 1; }
        long dlen = 0;
        unsigned char *d = lzd_decode(c, clen, n, &dlen);
        if (!d || dlen != n || (n && memcmp(d, in, n))) {
            printf("MISMATCH n=%ld dlen=%ld\n", n, dlen); return 1;
        }
        free(in); free(c); free(d);
    }
    printf("400 encode/decode round-trips clean\n");
    return 0;
}
