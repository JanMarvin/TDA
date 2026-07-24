#!/usr/bin/env python3
"""Build a TDA data archive: a Zoo container plus its description file.

TDA reads data archives through t_rzoo.c, which accepts Zoo packing methods 1
(LZD) and 2 (LZH) and explicitly rejects anything else, so a stored entry is
not an option and the payload has to be real LZW.  The parameters are read out
of lzd() in t_zoo.c: 9-bit initial codes packed least-significant-bit first,
CLEAR 256, EOF 257, first free code 258, widening to a maximum of 13 bits when
free_code reaches max_code.

The description file is parsed at t_rzoo.c:168.  Its first non-comment line is
the archive's own file name; every line after that is

    <number> <name> <type> <record-length> <records> <variables>

with type 1 or 2, records at least 1.

Usage: python3 tools/make_archive.py OUTDIR   (normally tests/fixtures)
"""
import struct, sys, os
from pathlib import Path

CLEAR, Z_EOF, FIRST_FREE, MAXBITS = 256, 257, 258, 13
ZOOTAG = 0xFDC4A7DC
ZOO_TEXT = b"ZOO 2.10 Archive.\032"


class BitWriter:
    """LSB-first bit packing, matching rd_dcode() in t_zoo.c."""

    def __init__(self):
        self.out = bytearray()
        self.acc = 0
        self.nbits = 0

    def write(self, code, width):
        self.acc |= (code & ((1 << width) - 1)) << self.nbits
        self.nbits += width
        while self.nbits >= 8:
            self.out.append(self.acc & 0xFF)
            self.acc >>= 8
            self.nbits -= 8

    def flush(self):
        if self.nbits:
            self.out.append(self.acc & 0xFF)
            self.acc = 0
            self.nbits = 0
        return bytes(self.out)


def lzd_compress(data):
    """LZW as lzd() decodes it, including its table-widening rule."""
    bw = BitWriter()
    table = {bytes([i]): i for i in range(256)}
    free_code, nbits, max_code = FIRST_FREE, 9, 512

    bw.write(CLEAR, nbits)
    if not data:
        bw.write(Z_EOF, nbits)
        return bw.flush()

    w = bytes([data[0]])
    for ch in data[1:]:
        wk = w + bytes([ch])
        if wk in table:
            w = wk
            continue
        bw.write(table[w], nbits)
        if free_code < (1 << MAXBITS):
            table[wk] = free_code
            free_code += 1
            # lzd() widens after adding, so the encoder must too
            if free_code >= max_code and nbits < MAXBITS:
                nbits += 1
                max_code <<= 1
        w = bytes([ch])
    bw.write(table[w], nbits)
    bw.write(Z_EOF, nbits)
    return bw.flush()


def dir_entry(next_ptr, offset, orig_size, now_size, name):
    """51-byte type 1 directory entry; field offsets from t_zoo.c."""
    e = bytearray(51)
    struct.pack_into("<I", e, 0, ZOOTAG)      # DTAG_I
    e[4] = 1                                   # DTYP_I  type 1
    e[5] = 1                                   # PKM_I   packing method lzd
    struct.pack_into("<i", e, 6, next_ptr)     # NXT_I
    struct.pack_into("<i", e, 10, offset)      # OFS_I
    struct.pack_into("<H", e, 14, 0x2821)      # DAT_I   arbitrary DOS date
    struct.pack_into("<H", e, 16, 0)           # TIM_I
    struct.pack_into("<H", e, 18, 0)           # CRC_I   not verified by TDA
    struct.pack_into("<i", e, 20, orig_size)   # ORGS_I
    struct.pack_into("<i", e, 24, now_size)    # SIZNOW_I
    e[28], e[29] = 2, 10                       # DMAJ_I / DMIN_I
    e[30] = 0                                  # DEL_I   not deleted
    e[31] = 0                                  # STRUC_I
    struct.pack_into("<i", e, 32, 0)           # CMT_I
    struct.pack_into("<H", e, 36, 0)           # CMTSIZ_I
    e[38:38 + len(name)] = name.encode()[:12]  # FNAME_I 13 bytes
    return bytes(e)


def build(outdir):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    # One fixed-length data file: 20 records of 3 fields, 8 columns each.
    reclen = 24
    nrec = 20
    payload = b"".join(
        ("%8d%8d%8d" % (i + 1, (i * 7) % 23 + 1, (i * 13) % 17 + 1)).encode()
        for i in range(nrec))
    assert len(payload) == reclen * nrec

    # A type 2 file describes the variables: name, file number, offset in the
    # record, and width.dec format.  check_avar() at t_rzoo.c:790 parses it,
    # and arcd refuses an archive that has none.
    # Records are fixed length here too: arcc checks record length against the
    # count declared in the description file and reports a mismatch otherwise.
    vrl = 40
    vlines = ["V1 1 0 8.0 first variable",
              "V2 1 8 8.0 second variable",
              "V3 1 16 8.0 third variable"]
    vardesc = b"".join(l.ljust(vrl - 1).encode() + b"\n" for l in vlines)

    packed = lzd_compress(payload)
    vpacked = lzd_compress(vardesc)

    hdr = bytearray(42)
    hdr[0:len(ZOO_TEXT)] = ZOO_TEXT
    struct.pack_into("<I", hdr, 20, ZOOTAG)
    data_at = 42
    dir_at = data_at + len(packed)
    struct.pack_into("<i", hdr, 24, dir_at)     # ZST_I
    struct.pack_into("<i", hdr, 28, -dir_at)    # ZSTM_I, must sum to zero
    hdr[32], hdr[33] = 2, 10                    # MAJV_I / MINV_I
    hdr[34] = 0                                 # HTYPE_I
    struct.pack_into("<i", hdr, 35, 0)          # ACMTPOS_I
    struct.pack_into("<H", hdr, 39, 0)          # ACMTLEN_I

    vdata_at = data_at + len(packed)
    dir_at = vdata_at + len(vpacked)
    struct.pack_into("<i", hdr, 24, dir_at)
    struct.pack_into("<i", hdr, 28, -dir_at)

    entry2_at = dir_at + 51
    term_at = entry2_at + 51
    entry = dir_entry(entry2_at, data_at, len(payload), len(packed), "adata.dat")
    entry2 = dir_entry(term_at, vdata_at, len(vardesc), len(vpacked), "avar.dat")
    # get_zoo() stops at the first entry whose next pointer is zero, and does
    # not treat that entry as a file, so a terminator has to follow the real
    # entries.
    term = dir_entry(0, 0, 0, 0, "")

    (out / "tda.zoo").write_bytes(
        bytes(hdr) + packed + vpacked + entry + entry2 + term)
    (out / "tda.zad").write_text(
        "# TDA archive description\n"
        "tda.zoo\n"
        "1 adata.dat 1 %d %d 3\n"
        "2 avar.dat 2 %d %d 0\n"
        % (reclen, nrec, vrl, len(vlines)))

    print("tda.zoo  %d bytes (%d raw -> %d packed)"
          % (42 + len(packed) + 102, len(payload), len(packed)))
    print("tda.zad  description for 1 file, %d records" % nrec)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else ".")
