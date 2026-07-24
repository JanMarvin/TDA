#!/usr/bin/env python3
"""Manual post-fixes for cases the automated pipeline can't handle."""
from pathlib import Path

src = Path(__file__).parent.parent / "src"

# t_eval.c: orphaned block comment tail after removed global definition
text = (src/"t_eval.c").read_text()
text = text.replace("                           this flag is zero. */\n",
                    "/* this flag is zero. */\n", 1)
(src/"t_eval.c").write_text(text)

# t_tri.c: voronoi has function-pointer parameter, breaks FUNC_DEF_RE
text = (src/"t_tri.c").read_text()
text = text.replace("int voronoi(struct Site *(*nextsite)());",
                    "int voronoi(TDAContext *ctx, struct Site *(*nextsite)());")
text = text.replace("int voronoi(struct Site *(*nextsite)())\n",
                    "int voronoi(TDAContext *ctx, struct Site *(*nextsite)())\n")
text = text.replace("    if (voronoi(nextone))",
                    "    if (voronoi(ctx, nextone))")
(src/"t_tri.c").write_text(text)

# t_zoo.c: TLONG (typedef) return type not caught by is_decl -> double ctx on to_long
text = (src/"t_zoo.c").read_text()
text = text.replace("TLONG to_long(TDAContext *ctx, ctx, unsigned char data[]);",
                    "TLONG to_long(TDAContext *ctx, unsigned char data[]);")
text = text.replace("TLONG to_long(TDAContext *ctx, ctx, unsigned char data[])\n",
                    "TLONG to_long(TDAContext *ctx, unsigned char data[])\n")
(src/"t_zoo.c").write_text(text)

print("Post-fixes applied.")
