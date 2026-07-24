/* t_gcmd.h */

#ifndef _TGCMD_H
#define _TGCMD_H

/*  functions in t_gcmd.c */

int mem_use(void);
int prntime(void);
int prn_txt(void);
int machp(int opt);
int t_parse(void);          
int data(int opt);

extern int MacroN;          /* number of macros                             */
extern char *MacroName[];   /* names of macros                              */
extern char *MacroDef[];    /* definition of macros                         */
extern char *MacroArgs[];   /* macro arguments                              */
extern short MacroNP[];     /* number of parameters                         */
extern int MacroNLen;       /* max length of macro names                    */
extern int MacroExLevel;    /* execution level of macro, used for local     */
extern int CurMacroNum;     /* current macor number, used for exists()      */

#endif /* _TGCMD_H */

/*  end of t_gcmd.h */











