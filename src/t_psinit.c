/****************************************************************************/
/*  t_psinit                                                                */
/*                                                                          */
/*  TDA. Program for Transition Data Analysis, written by Goetz Rohwer.     */
/*  Copyright (C) 1989,1991-97 Goetz Rohwer. All rights reserved.           */
/*                                                                          */
/*  This file is part of TDA.                                               */
/*                                                                          */
/*  TDA is free software; you can redistribute it and/or modify             */
/*  it under the terms of the GNU General Public License as published by    */
/*  the Free Software Foundation; either version 2 of the License, or       */
/*  (at your option) any later version.                                     */
/*                                                                          */
/*  TDA is distributed in the hope that it will be useful,                  */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*  GNU General Public License for more details.                            */
/*                                                                          */
/*  You should have received a copy of the GNU General Public License       */
/*  along with this program (it should be in a file named COPYING),         */
/*  if not, write to the Free Software Foundation, Inc.,                    */
/*  675 Mass Ave, Cambridge, MA 02139, USA.                                 */
/*                                                                          */

#include "tda.h"
#include "t_gen.h"
#include "t_plot.h"

/*  functions in t_psinit.c */

void ps_init(void);
void ps_inis(int opt);

/* ------------------------------------------------------------------------ */
/*  Some text for the PostScript header                                     */

char *PS_Head =
"/reencsmalldict 12 dict def % obtain dict for reencoding\n\
/ReEncodeSmall\n\
{\n\
    reencsmalldict begin\n\
    /newcodesandnames exch def\n\
    /newfontname exch def\n\
    /basefontname exch def\n\
    /basefontdict basefontname findfont def\n\
    /newfont basefontdict maxlength dict def\n\
    basefontdict\n\
    { exch dup /FID ne\n\
    { dup /Encoding eq\n\
    { exch dup length array copy newfont 3 1 roll put }\n\
    { exch newfont 3 1 roll put }\n\
    ifelse }\n\
    {pop pop }\n\
    ifelse }\n\
    forall\n\
    newfont /FontName newfontname put\n\
    newcodesandnames aload pop\n\
    newcodesandnames length 2 idiv\n\
    { newfont /Encoding get 3 1 roll put}\n\
    repeat newfontname newfont definefont pop\n\
    end \n\
} def\n\
/spanvec [\n\
8#0 /NUL\n\
8#1 /Eth\n\
8#2 /eth\n\
8#3 /Lslash\n\
8#4 /lslash\n\
8#5 /Scaron\n\
8#6 /scaron\n\
8#7 /Yacute\n\
8#10 /yacute\n\
8#11 /HT\n\
8#12 /LF\n\
8#13 /Thorn\n\
8#14 /thorn\n\
8#15 /CR\n\
8#16 /Zcaron\n\
8#17 /zcaron\n\
8#20 /DLE\n\
8#21 /DC1\n\
8#22 /DC2\n\
8#23 /DC3\n\
8#24 /DC4\n\
8#25 /onehalf\n\
8#26 /onequarter\n\
8#27 /onesuperior\n\
8#30 /threequarters\n\
8#31 /threesuperior\n\
8#32 /twosuperior\n\
8#33 /brokenbar\n\
8#34 /minus\n\
8#35 /multiply\n\
8#36 /RS\n\
8#37 /US\n\
8#200 /Adieresis\n\
8#201 /Aring\n\
8#202 /Ccedilla\n\
8#203 /Eacute\n\
8#204 /Ntilde\n\
8#205 /Odieresis\n\
8#206 /Udieresis\n\
8#207 /aacute\n\
8#210 /agrave\n\
8#211 /acircumflex\n\
8#212 /adieresis\n\
8#213 /atilde\n\
8#214 /aring\n\
8#215 /ccedilla\n\
8#216 /eacute\n\
8#217 /egrave\n\
8#220 /ecircumflex\n\
8#221 /edieresis\n\
8#222 /iacute\n\
8#223 /igrave\n\
8#224 /icircumflex\n\
8#225 /idieresis\n\
8#226 /ntilde\n\
8#227 /oacute\n\
8#230 /ograve\n\
8#231 /ocircumflex\n\
8#232 /odieresis\n\
8#233 /otilde\n\
8#234 /uacute\n\
8#235 /ugrave\n\
8#236 /ucircumflex\n\
8#237 /udieresis\n\
8#240 /dagger\n\
8#241 /degree\n\
8#242 /cent\n\
8#243 /sterling\n\
8#244 /section\n\
8#245 /bullet\n\
8#246 /paragraph\n\
8#247 /germandbls\n\
8#250 /registered\n\
8#251 /copyright\n\
8#252 /trademark\n\
8#253 /acute\n\
8#254 /dieresis\n\
8#255 /notequal\n\
8#256 /AE\n\
8#257 /Oslash\n\
8#260 /infinity\n\
8#261 /plusminus\n\
8#262 /lessequal\n\
8#263 /greaterequal\n\
8#264 /yen\n\
8#265 /mu\n\
8#266 /Gamma\n\
8#267 /summation\n\
8#270 /product\n\
8#271 /pi\n\
8#272 /integral\n\
8#273 /ordfeminine\n\
8#274 /ordmasculine\n\
8#275 /Omega\n\
8#276 /ae\n\
8#277 /oslash\n\
8#300 /questiondown\n\
8#301 /exclamdown\n\
8#302 /logicalnot\n\
8#303 /radical\n\
8#304 /florin\n\
8#305 /approxequal\n\
8#306 /Delta\n\
8#307 /guillemotleft\n\
8#310 /guillemotright\n\
8#311 /ellipsis\n\
8#312 /nbspace\n\
8#313 /Agrave\n\
8#314 /Atilde\n\
8#315 /Otilde\n\
8#316 /OE\n\
8#317 /oe\n\
8#320 /endash\n\
8#321 /emdash\n\
8#322 /quotedblleft\n\
8#323 /quotedblright\n\
8#324 /quoteleft\n\
8#325 /quoteright\n\
8#326 /divide\n\
8#327 /lozenge\n\
8#330 /ydieresis\n\
8#331 /Ydieresis\n\
8#332 /fraction\n\
8#333 /currency\n\
8#334 /guilsinglleft\n\
8#335 /guilsinglright\n\
8#336 /fi\n\
8#337 /fl\n\
8#340 /daggerdbl\n\
8#341 /periodcentered\n\
8#342 /quotesinglbase\n\
8#343 /quotedblbase\n\
8#344 /perthousand\n\
8#345 /Acircumflex\n\
8#346 /Ecircumflex\n\
8#347 /Aacute\n\
8#350 /Edieresis\n\
8#351 /Egrave\n\
8#352 /Iacute\n\
8#353 /Icircumflex\n\
8#354 /Idieresis\n\
8#355 /Igrave\n\
8#356 /Oacute\n\
8#357 /Ocircumflex\n\
8#360 /apple\n\
8#361 /Ograve\n\
8#362 /Uacute\n\
8#363 /Ucircumflex\n\
8#364 /Ugrave\n\
8#365 /dotlessi\n\
8#366 /circumflex\n\
8#367 /tilde\n\
8#370 /macron\n\
8#371 /breve\n\
8#372 /dotaccent\n\
8#373 /ring\n\
8#374 /cedilla\n\
8#375 /hungarumlaut\n\
8#376 /ogonek\n\
8#377 /caron\n\
] def\n";

/*****************
/xsym {\n\
    gsave\n\
    /s ssiz 0.71 mul def\n\
    s neg s neg rm s 2 mul s 2 mul rl\n\
    s 2 mul neg 0 rm 0 s 2 mul s 2 mul neg rl\n\
    stroke\n    grestore\n\
} def\n\
*********************/


char *PS1_Head = 
"/invers {\n\
    1 setgray fill 0 setgray\n\
} def\n\
/circle {\n\
    newpath ssiz 0 360 arc\n\
} def\n\
/cross {\n\
    gsave\n\
    ssiz neg 0 rm ssiz 2 mul 0 rl\n\
    ssiz neg ssiz rm 0 ssiz 2 mul neg rl\n\
    stroke\n    grestore\n\
} def\n\
/xsym {\n\
    gsave\n\
    45 rotate\
    ssiz neg 0 rm ssiz 2 mul 0 rl\n\
    ssiz neg ssiz rm 0 ssiz 2 mul neg rl\n\
    stroke\n    grestore\n\
} def\n\
/square {\n\
    gsave\n\
    ssiz neg ssiz neg rm ssiz 2 mul 0 rl\n\
    0 ssiz 2 mul rl ssiz 2 mul neg 0 rl closepath\n\
    gsave 1 setgray fill 0 setgray grestore\n\
} def\n\
/trian1 {\n\
    gsave\n\
    ssiz neg ssiz neg rm ssiz 2 mul 0 rl\n\
    ssiz neg ssiz 2 mul rl closepath\n\
    gsave 1 setgray fill 0 setgray grestore\n\
} def\n\
/trian2 {\n\
    gsave\n\
    ssiz neg ssiz rm ssiz ssiz 2 mul neg rl\n\
    ssiz ssiz 2 mul rl closepath\n\
    gsave 1 setgray fill 0 setgray grestore\n\
} def\n\
/rhomb {\n\
    gsave\n\
    ssiz neg 0 rm ssiz ssiz neg rl\n\
    ssiz ssiz rl ssiz neg ssiz rl closepath\n\
    gsave 1 setgray fill 0 setgray grestore\n\
} def\n";

/* ------------------------------------------------------------------------ */
/*  ps_init()   Initialization sequence for PostScript output.              */

void ps_init(void)
{
    fprintf(PSFd,"%%!PS-Adobe-2.0 EPSF-1.2\n");
    fprintf(PSFd,"%%%%Title: %s\n",PSFName);
    fprintf(PSFd,"%%%%Creator: TDA %3.1f\n",TDA_Version);
    fprintf(PSFd,"%%%%CreationDate: ");
#if TIME_ON
    prn_time(PSFd);
#else
    fprintf(PSFd,"\n");
#endif
    fprintf(PSFd,"%%%%BoundingBox: 00000 00000 00000 00000\n");
    fprintf(PSFd,"%%%%BeginProlog\n");
    fprintf(PSFd,"/TDAdict 200 dict def %% define local dictionary\n");
    fprintf(PSFd,"TDAdict begin         %% push dictionary onto the dictionary stack\n");

    fprintf(PSFd,"%s",PS_Head);

    fprintf(PSFd,"/m  {moveto} bind def\n");
    fprintf(PSFd,"/rm {rmoveto} bind def\n");
    fprintf(PSFd,"/l  {lineto} bind def\n");
    fprintf(PSFd,"/rl {rlineto} bind def\n\n");
    /*************
    fprintf(PSFd,"/xorg 0 def %% X and Y coordinates of origin\n");
    fprintf(PSFd,"/yorg 0 def\n");
    fprintf(PSFd,"/sfx  0 def %% x scaling\n");
    fprintf(PSFd,"/sfy  0 def %% y scaling\n");
    fprintf(PSFd,"/rf   0 def %% rotation\n\n");
    *********/
    fprintf(PSFd,"/fsiz 0 def %% font size\n");
    fprintf(PSFd,"/FT {/TFont findfont fsiz scalefont setfont} bind def\n");
    fprintf(PSFd,"/FS {/Symbol findfont fsiz scalefont setfont} bind def\n\n");

    fprintf(PSFd,"/ssiz 1 def %% half size of symbol\n");
    fprintf(PSFd,"/lpt  0 def %% size of numerical labels\n");
    fprintf(PSFd,"/str 20 string def\n");
    fprintf(PSFd,"/xmin 0 def\n\n");

    fprintf(PSFd,"/sclear { %% clear area for string\n");
    fprintf(PSFd,"gsave currentpoint pop -2 add fsiz 4 div m\n");
    fprintf(PSFd,"dup stringwidth exch 4 add exch rl\n");
    fprintf(PSFd,"0 setlinecap fsiz 2 add setlinewidth\n");
    fprintf(PSFd,"1 setgray stroke grestore\n");
    fprintf(PSFd,"} def\n");

    fprintf(PSFd,"/center { %% center a string\n");
    fprintf(PSFd,"    dup stringwidth pop\n");
    fprintf(PSFd,"    2 div neg 0 rm\n");
    fprintf(PSFd,"} def\n");

    fprintf(PSFd,"/aright { %% right justified\n");
    fprintf(PSFd,"    dup stringwidth pop\n");
    fprintf(PSFd,"    neg 0 rm\n");
    fprintf(PSFd,"} def\n");

    fprintf(PSFd,"/adjust { %% adjust for y axis labels\n");
    fprintf(PSFd,"    dup stringwidth pop\n");
    fprintf(PSFd,"    neg lpt 3 div neg rm\n"); 
    fprintf(PSFd,"    currentpoint pop xmin lt {currentpoint pop /xmin exch def} if\n");
    fprintf(PSFd,"} def\n");

    fprintf(PSFd,"%s",PS1_Head); /* additional text for the header */

    fprintf(PSFd,"end         %% pop TDAdict off the dictionary stack\n");
    fprintf(PSFd,"%%%%EndProlog\n");
    fprintf(PSFd,"TDAdict begin\n");
    fprintf(PSFd,"/Times-Roman /TFont spanvec ReEncodeSmall\n");
}

/* ------------------------------------------------------------------------ */
/*  ps_inis(opt)   Set physical size etc of current coordinate system.      */

void ps_inis(int opt)
{
    fprintf(PSFd,"\n%% define current coordinate system\n");
    fprintf(PSFd,"%% xorg=%d yorg=%d pxlen=%5.2f pylen=%5.2f\n",
                                                XOrg,YOrg,PXLen,PYLen);

    fprintf(PSFd,"/xorg %d def\n",XOrg - PSXOrg);
    fprintf(PSFd,"/yorg %d def\n",YOrg - PSYOrg);
    fprintf(PSFd,"xorg yorg translate\n");

    /* update bounding box */

    if (opt) {
        if (BBLX > XOrg)
            BBLX = XOrg;
        if (BBLY > YOrg)
            BBLY = YOrg;
        if (BBUX < XOrg + (int)(PSXLen + 0.5)) 
            BBUX = XOrg + (int)(PSXLen + 0.5);
        if (BBUY < YOrg + (int)(PSYLen + 0.5)) 
            BBUY = YOrg + (int)(PSYLen + 0.5);
    }
    PSXOrg = XOrg;
    PSYOrg = YOrg;

    fprintf(PSFd,"/rf %5.2f def rf rotate\n",ROTFac);
    fprintf(PSFd,"/sfx %5.2f def\n",SCALXFac);
    fprintf(PSFd,"/sfy %5.2f def\n",SCALYFac);
    fprintf(PSFd,"sfx sfy scale\n\n");
    fprintf(PSFd,"%%#Parameter: %d %d %g %g %g %g %g %g %d %d %g %g %g\n",
        XOrg,YOrg,PXLen,PYLen,PA1[0],PA1[1],PA2[0],PA2[1],PSLog[0],PSLog[1],
        SCALXFac,SCALYFac,ROTFac);

    fprintf(PSFd,"1 setlinecap\n");
    fprintf(PSFd,"1 setlinejoin\n\n");
}


