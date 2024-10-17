# Algorithm 814: Fortran 90 software for floating-point multiple precision arithmetic, gamma and related functions

> David M. Smith. 2001.  
> Algorithm 814: Fortran 90 software for floating-point multiple precision arithmetic, gamma and related functions.  
> ACM Trans. Math. Softw. 27, 4 (December 2001), 377–387.  
> https://doi.org/10.1145/504210.504211

## FM Package

```for
!     FM 1.2                    David M. Smith                  8-17-01

!  The routines in this package perform multiple precision arithmetic and
!  functions on three kinds of numbers.
!  FM routines handle floating-point real multiple precision numbers,
!  IM routines handle integer multiple precision numbers, and
!  ZM routines handle floating-point complex multiple precision numbers.
```

### LIST OF ROUTINES

```for

!  These are the FM routines that are designed to be called by the user.
!  All are subroutines except logical function FMCOMP.
!  MA, MB, MC refer to FM format numbers.

!  In Fortran-90 and later versions of the Fortran standard, it is potentially
!  unsafe to use the same array more than once in the calling sequence.  The
!  operation MA = MA + MB should not be written as
!        CALL FMADD(MA,MB,MA)
!  since the compiler is allowed to pass the three arguments with a copy in /
!  copy out mechanism.  This means the third argument, containing the result,
!  might not be copied out last, and then a later copy out of the original
!  input MA could destroy the computed result.

!  One solution is to use a third array and then put the result back in MA:
!        CALL FMADD(MA,MB,MC)
!        CALL FMEQ(MC,MA)

!  When the first call is doing one of the "fast" operations like addition,
!  the extra call to move the result back to MA can cause a noticeable loss in
!  efficiency.  To avoid this, separate routines are provided for the basic
!  arithmetic operations when the result is to be returned in the same array
!  as one of the inputs.

!  A routine name with a suffix of  "_R1" returns the result in the first
!  input array, and a suffix of "_R2" returns the result in the second input
!  array.  The example above would then be:
!        CALL FMADD_R1(MA,MB)

!  These routines each have one less argument than the original version, since
!  the output is re-directed to one of the inputs.  The result array should
!  not be the same as any input array when the original version of the routine
!  is used.

!  The routines that can be used this way are listed below.  For others, like
!        CALL FMEXP(MA,MA)
!  the relative cost of doing an extra copy is small.  This one should become
!        CALL FMEXP(MA,MB)
!        CALL FMEQ(MB,MA)

!  If the derived-type interface is used, as in
!        TYPE (FM) A,B
!        ...
!        A = A + B
!  there is no problem putting the result back into A, since the interface
!  routine creates a temporary scratch array for the result of A + B, allowing
!  copy in / copy out to work.

!  For each of these routines there is also a version available for which the
!  argument list is the same but all FM numbers are in packed format.  The
!  routines using packed numbers have the same names except 'FM' is replaced by
!  'FP' at the start of each name.


!  FMABS(MA,MB)         MB = ABS(MA)

!  FMACOS(MA,MB)        MB = ACOS(MA)

!  FMADD(MA,MB,MC)      MC = MA + MB

!  FMADD_R1(MA,MB)      MA = MA + MB

!  FMADD_R2(MA,MB)      MB = MA + MB

!  FMADDI(MA,IVAL)      MA = MA + IVAL   Increment an FM number by a one word
!                                        integer.  Note this call does not have
!                                        an "MB" result like FMDIVI and FMMPYI.

!  FMASIN(MA,MB)        MB = ASIN(MA)

!  FMATAN(MA,MB)        MB = ATAN(MA)

!  FMATN2(MA,MB,MC)     MC = ATAN2(MA,MB)

!  FMBIG(MA)            MA = Biggest FM number less than overflow.

!  FMCHSH(MA,MB,MC)     MB = COSH(MA),  MC = SINH(MA).
!                            Faster than making two separate calls.

!  FMCOMP(MA,LREL,MB)        Logical comparison of MA and MB.
!                            LREL is a CHARACTER*2 value identifying
!                            which of the six comparisons is to be made.
!                            Example:  IF (FMCOMP(MA,'GE',MB)) ...
!                            Also can be:  IF (FMCOMP(MA,'>=',MB)) ...
!                            CHARACTER*1 is ok:  IF (FMCOMP(MA,'>',MB)) ...

!  FMCONS                    Set several saved constants that depend on MBASE,
!                            the base being used.  FMCONS should be called
!                            immediately after changing MBASE.

!  FMCOS(MA,MB)         MB = COS(MA)

!  FMCOSH(MA,MB)        MB = COSH(MA)

!  FMCSSN(MA,MB,MC)     MB = COS(MA),  MC = SIN(MA).
!                            Faster than making two separate calls.

!  FMDIG(NSTACK,KST)         Find a set of precisions to use during Newton
!                            iteration for finding a simple root starting with
!                            about double precision accuracy.

!  FMDIM(MA,MB,MC)      MC = DIM(MA,MB)

!  FMDIV(MA,MB,MC)      MC = MA / MB

!  FMDIV_R1(MA,MB)      MA = MA / MB

!  FMDIV_R2(MA,MB)      MB = MA / MB

!  FMDIVI(MA,IVAL,MB)   MB = MA/IVAL   IVAL is a one word integer.

!  FMDIVI_R1(MA,IVAL)   MA = MA/IVAL

!  FMDP2M(X,MA)         MA = X    Convert from double precision to FM.

!  FMDPM(X,MA)          MA = X    Convert from double precision to FM.
!                                 Faster than FMDP2M, but MA agrees with X only
!                                 to D.P. accuracy.  See the comments in the
!                                 two routines.

!  FMEQ(MA,MB)          MB = MA   Both have precision NDIG.
!                                 This is the version to use for standard
!                                 B = A  statements.

!  FMEQU(MA,MB,NA,NB)   MB = MA   Version for changing precision.
!                                 MA has NA digits (i.e., MA was computed
!                                 using NDIG = NA), and MB will be defined
!                                 having NB digits.
!                                 MB is rounded if NB < NA
!                                 MB is zero-padded if NB > NA

!  FMEXP(MA,MB)         MB = EXP(MA)

!  FMFLAG(K)            K = KFLAG  get the value of the FM condition
!                                  flag -- stored in the internal FM
!                                  variable KFLAG in module FMVALS.

!  FMFORM(FORM,MA,STRING)    MA is converted to a character string using format
!                               FORM and returned in STRING.  FORM can represent
!                               I, F, E, or 1PE formats.  Example:
!                               CALL FMFORM('F60.40',MA,STRING)

!  FMFPRT(FORM,MA)           Print MA on unit KW using FORM format.

!  FMI2M(IVAL,MA)       MA = IVAL   Convert from one word integer to FM.

!  FMINP(LINE,MA,LA,LB) MA = LINE   Input conversion.
!                                   Convert LINE(LA) through LINE(LB)
!                                   from characters to FM.

!  FMINT(MA,MB)         MB = INT(MA)    Integer part of MA.

!  FMIPWR(MA,IVAL,MB)   MB = MA**IVAL   Raise an FM number to a one word
!                                       integer power.

!  FMLG10(MA,MB)        MB = LOG10(MA)

!  FMLN(MA,MB)          MB = LOG(MA)

!  FMLNI(IVAL,MA)       MA = LOG(IVAL)   Natural log of a one word integer.

!  FMM2DP(MA,X)         X  = MA     Convert from FM to double precision.

!  FMM2I(MA,IVAL)       IVAL = MA   Convert from FM to integer.

!  FMM2SP(MA,X)         X  = MA     Convert from FM to single precision.

!  FMMAX(MA,MB,MC)      MC = MAX(MA,MB)

!  FMMIN(MA,MB,MC)      MC = MIN(MA,MB)

!  FMMOD(MA,MB,MC)      MC = MA mod MB

!  FMMPY(MA,MB,MC)      MC = MA * MB

!  FMMPY_R1(MA,MB)      MA = MA * MB

!  FMMPY_R2(MA,MB)      MB = MA * MB

!  FMMPYI(MA,IVAL,MB)   MB = MA*IVAL    Multiply by a one word integer.

!  FMMPYI_R1(MA,IVAL)   MA = MA*IVAL

!  FMNINT(MA,MB)        MB = NINT(MA)   Nearest FM integer.

!  FMOUT(MA,LINE,LB)    LINE = MA   Convert from FM to character.
!                                   LINE is a character array of length LB.

!  FMPI(MA)             MA = pi

!  FMPRNT(MA)                Print MA on unit KW using current format.

!  FMPWR(MA,MB,MC)      MC = MA**MB

!  FM_RANDOM_NUMBER(X)  X    is returned as a double precision random number,
!                            uniform on (0,1).  High-quality, long-period
!                            generator.
!                            Note that X is double precision, unlike the similar
!                            Fortran intrinsic random number routine, which
!                            returns a single-precision result.
!                            See the comments in section 10 below and also those
!                            in the routine for more details.

!  FMREAD(KREAD,MA)     MA   is returned after reading one (possibly multi-line)
!                            FM number on unit KREAD.  This routine reads
!                            numbers written by FMWRIT.

!  FMRPWR(MA,K,J,MB)    MB = MA**(K/J)  Rational power.
!                            Faster than FMPWR for functions like the cube root.

!  FMSET(NPREC)              Set the internal FM variables so that the precision
!                            is at least NPREC base 10 digits plus three base 10
!                            guard digits.

!  FMSETVAR(STRING)          Define a new value for one of the internal FM
!                            variables in module FMVALS that controls one of the
!                            FM options.  STRING has the form  variable = value.
!                            Example:  To change the screen width for FM output:
!                                  CALL FMSETVAR(' KSWIDE = 120 ')
!                            The variables that can be changed and the options
!                            they control are listed in sections 2 through 6
!                            above.  Only one variable can be set per call.
!                            The variable name in STRING must have no embedded
!                            blanks.  The value part of STRING can be in any
!                            numerical format, except in the case of variable
!                            CMCHAR, which is character type.  To set CMCHAR to
!                            'E', don't use any quotes in STRING:
!                                  CALL FMSETVAR(' CMCHAR = E ')

!  FMSIGN(MA,MB,MC)     MC = SIGN(MA,MB)   Sign transfer.

!  FMSIN(MA,MB)         MB = SIN(MA)

!  FMSINH(MA,MB)        MB = SINH(MA)

!  FMSP2M(X,MA)         MA = X   Convert from single precision to FM.

!  FMSQR(MA,MB)         MB = MA * MA   Faster than FMMPY.

!  FMSQR_R1(MA)         MA = MA * MA

!  FMSQRT(MA,MB)        MB = SQRT(MA)

!  FMSQRT_R1(MA)        MA = SQRT(MA)

!  FMST2M(STRING,MA)    MA = STRING
!                            Convert from character string to FM.
!                            STRING may be in any numerical format.
!                            Often more convenient than FMINP, which converts
!                            an array of CHARACTER*1 values.  Example:
!                                  CALL FMST2M('123.4',MA)

!  FMSUB(MA,MB,MC)      MC = MA - MB

!  FMSUB_R1(MA,MB)      MA = MA - MB

!  FMSUB_R2(MA,MB)      MB = MA - MB

!  FMTAN(MA,MB)         MB = TAN(MA)

!  FMTANH(MA,MB)        MB = TANH(MA)

!  FMULP(MA,MB)         MB = One Unit in the Last Place of MA.

!  FMVARS                    Write the current values of the internal FM
!                            variables on unit KW.

!  FMWRIT(KWRITE,MA)         Write MA on unit KWRITE.
!                            Multi-line numbers will have '&' as the last
!                            nonblank character on all but the last line.  These
!                            numbers can then be read easily using FMREAD.



!  These are the Gamma and Related Functions.

!  FMBERN(N,MA,MB)      MB = MA*B(N)  Multiply by Nth Bernoulli number

!  FMBETA(MA,MB,MC)     MC = Beta(MA,MB)

!  FMCOMB(MA,MB,MC)     MC = Combination MA choose MB  (Binomial coeff.)

!  FMEULR(MA)           MA = Euler's constant ( 0.5772156649... )

!  FMFACT(MA,MB)        MB = MA Factorial  (Gamma(MA+1))

!  FMGAM(MA,MB)         MB = Gamma(MA)

!  FMIBTA(MX,MA,MB,MC)  MC = Incomplete Beta(MX,MA,MB)

!  FMIGM1(MA,MB,MC)     MC = Incomplete Gamma(MA,MB).  Lower case Gamma(a,x)

!  FMIGM2(MA,MB,MC)     MC = Incomplete Gamma(MA,MB).  Upper case Gamma(a,x)

!  FMLNGM(MA,MB)        MB = Ln(Gamma(MA))

!  FMPGAM(N,MA,MB)      MB = Polygamma(N,MA)  (Nth derivative of Psi)

!  FMPOCH(MA,N,MB)      MB = MA*(MA+1)*(MA+2)*...*(MA+N-1)  (Pochhammer)

!  FMPSI(MA,MB)         MB = Psi(MA)      (Derivative of Ln(Gamma(MA))



!  These are the integer routines that are designed to be called by the user.
!  All are subroutines except logical function IMCOMP.  MA, MB, MC refer to IM
!  format numbers.  In each case the version of the routine to handle packed IM
!  numbers has the same name, with 'IM' replaced by 'IP'.

!  IMABS(MA,MB)         MB = ABS(MA)

!  IMADD(MA,MB,MC)      MC = MA + MB

!  IMBIG(MA)            MA = Biggest IM number less than overflow.

!  IMCOMP(MA,LREL,MB)        Logical comparison of MA and MB.
!                            LREL is a CHARACTER*2 value identifying which of
!                            the six comparisons is to be made.
!                            Example:  IF (IMCOMP(MA,'GE',MB)) ...
!                            Also can be:  IF (IMCOMP(MA,'>=',MB))
!                            CHARACTER*1 is ok:  IF (IMCOMP(MA,'>',MB)) ...

!  IMDIM(MA,MB,MC)      MC = DIM(MA,MB)

!  IMDIV(MA,MB,MC)      MC = int(MA/MB)
!                            Use IMDIVR if the remainder is also needed.

!  IMDIVI(MA,IVAL,MB)   MB = int(MA/IVAL)
!                            IVAL is a one word integer.
!                            Use IMDVIR to get the remainder also.

!  IMDIVR(MA,MB,MC,MD)  MC = int(MA/MB),   MD = MA mod MB
!                            When both the quotient and remainder are needed,
!                            this routine is twice as fast as calling both
!                            IMDIV and IMMOD.

!  IMDVIR(MA,IVAL,MB,IREM)   MB = int(MA/IVAL),   IREM = MA mod IVAL
!                            IVAL and IREM are one word integers.

!  IMEQ(MA,MB)          MB = MA

!  IMFM2I(MAFM,MB)      MB = MAFM  Convert from real (FM) format to
!                                  integer (IM) format.

!  IMFORM(FORM,MA,STRING)    MA is converted to a character string using format
!                               FORM and returned in STRING.  FORM can represent
!                               I, F, E, or 1PE formats.  Example:
!                               CALL IMFORM('I70',MA,STRING)

!  IMFPRT(FORM,MA)           Print MA on unit KW using FORM format.

!  IMGCD(MA,MB,MC)      MC = greatest common divisor of MA and MB.

!  IMI2FM(MA,MBFM)      MBFM = MA  Convert from integer (IM) format to
!                                  real (FM) format.

!  IMI2M(IVAL,MA)       MA = IVAL   Convert from one word integer to IM.

!  IMINP(LINE,MA,LA,LB) MA = LINE   Input conversion.
!                                   Convert LINE(LA) through LINE(LB)
!                                   from characters to IM.

!  IMM2DP(MA,X)         X  = MA     Convert from IM to double precision.

!  IMM2I(MA,IVAL)       IVAL = MA   Convert from IM to one word integer.

!  IMMAX(MA,MB,MC)      MC = MAX(MA,MB)

!  IMMIN(MA,MB,MC)      MC = MIN(MA,MB)

!  IMMOD(MA,MB,MC)      MC = MA mod MB

!  IMMPY(MA,MB,MC)      MC = MA*MB

!  IMMPYI(MA,IVAL,MB)   MB = MA*IVAL    Multiply by a one word integer.

!  IMMPYM(MA,MB,MC,MD)  MD = MA*MB mod MC
!                            Slightly faster than calling IMMPY and IMMOD
!                            separately, and it works for cases where IMMPY
!                            would return OVERFLOW.

!  IMOUT(MA,LINE,LB)    LINE = MA   Convert from IM to character.
!                                   LINE is a character array of length LB.

!  IMPMOD(MA,MB,MC,MD)  MD = MA**MB mod MC

!  IMPRNT(MA)                Print MA on unit KW.

!  IMPWR(MA,MB,MC)      MC = MA**MB

!  IMREAD(KREAD,MA)     MA   is returned after reading one (possibly multi-line)
!                            IM number on unit KREAD.
!                            This routine reads numbers written by IMWRIT.

!  IMSIGN(MA,MB,MC)     MC = SIGN(MA,MB)   Sign transfer.

!  IMSQR(MA,MB)         MB = MA*MA   Faster than IMMPY.

!  IMST2M(STRING,MA)    MA = STRING
!                            Convert from character string to IM.
!                            Often more convenient than IMINP, which converts an
!                            array of CHARACTER*1 values.  Example:
!                                 CALL IMST2M('12345678901',MA)

!  IMSUB(MA,MB,MC)      MC = MA - MB

!  IMWRIT(KWRITE,MA)         Write MA on unit KWRITE.
!                            Multi-line numbers will have '&' as the last
!                            nonblank character on all but the last line.
!                            These numbers can then be read easily using IMREAD.



!  These are the complex routines that are designed to be called by the user.
!  All are subroutines, and in each case the version of the routine to handle
!  packed ZM numbers has the same name, with 'ZM' replaced by 'ZP'.

!  MA, MB, MC refer to ZM format complex numbers.
!  MAFM, MBFM, MCFM refer to FM format real numbers.
!  INTEG is a Fortran INTEGER variable.
!  ZVAL is a Fortran COMPLEX variable.

!  ZMABS(MA,MBFM)       MBFM = ABS(MA)    Result is real.

!  ZMACOS(MA,MB)        MB = ACOS(MA)

!  ZMADD(MA,MB,MC)      MC = MA + MB

!  ZMADDI(MA,INTEG)     MA = MA + INTEG  Increment an ZM number by a one word
!                                        integer.  Note this call does not have
!                                        an "MB" result like ZMDIVI and ZMMPYI.

!  ZMARG(MA,MBFM)       MBFM = Argument(MA)    Result is real.

!  ZMASIN(MA,MB)        MB = ASIN(MA)

!  ZMATAN(MA,MB)        MB = ATAN(MA)

!  ZMCHSH(MA,MB,MC)     MB = COSH(MA),  MC = SINH(MA).
!                            Faster than 2 calls.

!  ZMCMPX(MAFM,MBFM,MC) MC = CMPLX(MAFM,MBFM)

!  ZMCONJ(MA,MB)        MB = CONJG(MA)

!  ZMCOS(MA,MB)         MB = COS(MA)

!  ZMCOSH(MA,MB)        MB = COSH(MA)

!  ZMCSSN(MA,MB,MC)     MB = COS(MA),  MC = SIN(MA).
!                            Faster than 2 calls.

!  ZMDIV(MA,MB,MC)      MC = MA / MB

!  ZMDIVI(MA,INTEG,MB)  MB = MA / INTEG

!  ZMEQ(MA,MB)          MB = MA

!  ZMEQU(MA,MB,NDA,NDB) MB = MA    Version for changing precision.
!                                  (NDA and NDB are as in FMEQU)

!  ZMEXP(MA,MB)         MB = EXP(MA)

!  ZMFORM(FORM1,FORM2,MA,STRING)   STRING = MA
!                       MA is converted to a character string using format
!                       FORM1 for the real part and FORM2 for the imaginary
!                       part.  The  result is returned in STRING.  FORM1 and
!                       FORM2 can represent I, F, E, or 1PE formats.  Example:
!                             CALL ZMFORM('F20.10','F15.10',MA,STRING)
!                       A 1PE in the first format does not carry over to the
!                       other format descriptor, as it would in an ordinary
!                       FORMAT statement.

!  ZMFPRT(FORM1,FORM2,MA)    Print MA on unit KW using formats FORM1 and FORM2.

!  ZMI2M(INTEG,MA)           MA = CMPLX(INTEG,0)

!  ZM2I2M(INTEG1,INTEG2,MA)  MA = CMPLX(INTEG1,INTEG2)

!  ZMIMAG(MA,MBFM)           MBFM = IMAG(MA)    Imaginary part.

!  ZMINP(LINE,MA,LA,LB)      MA = LINE   Input conversion.
!                                 Convert LINE(LA) through LINE(LB) from
!                                 characters to ZM.  LINE is a character array
!                                 of length at least LB.

!  ZMINT(MA,MB)         MB = INT(MA)        Integer part of both Real
!                                           and Imaginary parts of MA.

!  ZMIPWR(MA,INTEG,MB)  MB = MA ** INTEG    Integer power function.

!  ZMLG10(MA,MB)        MB = LOG10(MA)

!  ZMLN(MA,MB)          MB = LOG(MA)

!  ZMM2I(MA,INTEG)      INTEG = INT(REAL(MA))

!  ZMM2Z(MA,ZVAL)       ZVAL = MA

!  ZMMPY(MA,MB,MC)      MC = MA * MB

!  ZMMPYI(MA,INTEG,MB)  MB = MA * INTEG

!  ZMNINT(MA,MB)        MB = NINT(MA)   Nearest integer of both Real
!                                       and Imaginary.

!  ZMOUT(MA,LINE,LB,LAST1,LAST2)        LINE = MA
!                       Convert from FM to character.
!                       LINE  is the returned character*1 array.
!                       LB    is the dimensioned size of LINE.
!                       LAST1 is returned as the position in LINE of
!                             the last character of REAL(MA).
!                       LAST2 is returned as the position in LINE
!                             of the last character of AIMAG(MA).

!  ZMPRNT(MA)           Print MA on unit KW using current format.

!  ZMPWR(MA,MB,MC)      MC = MA ** MB

!  ZMREAD(KREAD,MA)     MA   is returned after reading one (possibly multi-line)
!                            ZM number on unit KREAD.
!                            This routine reads numbers written by ZMWRIT.

!  ZMREAL(MA,MBFM)      MBFM = REAL(MA)    Real part.

!  ZMRPWR(MA,IVAL,JVAL,MB)     MB = MA ** (IVAL/JVAL)

!  ZMSET(NPREC)         Set precision to the equivalent of a few more than NPREC
!                       base 10 digits.  This is now the same as FMSET, but is
!                       retained for compatibility with earlier versions of the
!                       package.

!  ZMSIN(MA,MB)         MB = SIN(MA)

!  ZMSINH(MA,MB)        MB = SINH(MA)

!  ZMSQR(MA,MB)         MB = MA*MA    Faster than ZMMPY.

!  ZMSQRT(MA,MB)        MB = SQRT(MA)

!  ZMST2M(STRING,MA)    MA = STRING
!                            Convert from character string to ZM.
!                            Often more convenient than ZMINP, which
!                            converts an array of CHARACTER*1 values.
!                            Example: CALL ZMST2M('123.4+5.67i',MA).

!  ZMSUB(MA,MB,MC)      MC = MA - MB

!  ZMTAN(MA,MB)         MB = TAN(MA)

!  ZMTANH(MA,MB)        MB = TANH(MA)

!  ZMWRIT(KWRITE,MA)    Write MA on unit KWRITE.  Multi-line numbers are
!                       formatted for automatic reading with ZMREAD.

!  ZMZ2M(ZVAL,MA)       MA = ZVAL
```
