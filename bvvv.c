
/******************************************************************************
*
*	Description:
*	___________
*
*       This program takes a reddening estimate from E(B-V) and calculates Teff from (B-V)o,
*       returning a dereddened observed SED, that expected from a star with input spectral
*       type, as well as a face-on reprocessing disk model from HSKV (1992).
*
*       Watts cm^-2 mu^-1 were taken from LAH as follows -
*
*       UBVRcIc from Bessell 1979, AJ, V91, 589.
*       JHKLM from Campins, Rieke, and Lebofsky 1985, AJ V90 896.
*       NQ from Rieke, Lebofsky, and Low 1985, AJ, V90 904.
*       and IRAS12,25,60, and 100 from IRAS Exp. Supp. (color-corrected!).
*
*       Intrinsic stellar SEDs were taken from -
*
*       Teff B8-K5 from Schmitd-Kaler.  Teff K5-M6 from Bessell (1991, AJ, V101, 662).
*       J-H, H-K, V-K from Bessell and Brett (1988, PASP, V100, 1134).
*       >K5 R-I, V-I Bessell (1991); <K5 R-I, V-I from Johnson (1966, ARAA, V4, 193).
*       With color corrections Johnson to Cousins from Besell (1979, AJ, V91, 589).
*       B-V from Johnson B8-K5; Bessell (1991) K7-M6.
*       U-V from Johnson B8-M6.
*       (V-12) from Waters et al. (1987) A&A,  V172, 225.
*
*       Extinction law taken from Rieke and Lebofsky (ApJ, 1985, V288, 618).
*
*       Reprocessing disk models were taken from Hillenbrand et al. (1992, APJ, 397, 613)
*       and presented here for spectral types A0, F0, G0, K0, K5, M0, and M3.
*
*	Call Sequence:
*	______________
*
*	> cc bvv.c -o bvv -lm
*
*	Parameters:
*	___________
*
*       Input star name, distance modulus, estimate of E(B-V),
*       U,B,V,Rc,Ic,J,H,K,L,M,N,12,Q,20,60, and 100  IN MAGNITUDES!  In order
*       to convert Janskies to magnitudes, need flux for zero-mag star:
*
*       Band        U     B      V          Rc        Ic        J
*       Fnu(Jy)     1810  4260   3640       3080      2550      1603
*       Band        H     K      L          M         N         Q
*       Fnu(Jy)     1075  667    288        170       36.0      9.4
*       Band        IRAS12     IRAS25    IRAS60    IRAS100
*       Fnu(Jy)     28.3       6.73      1.19      0.43
*
*	Returns:
*	________
*
*       Dereddened observed SED as well as expected stellar SED based on
*       spectral type, normalized to match at I-band.
*
*	Notes:
*	______
*
*       Based on the program bigsed.c by MRM.
*
*	By:
*	___
*
*       ALM 6-20-2005     Steward Observatory
*
*******************************************************************************/

#include <stdio.h>
#include <math.h>

main()

{    double dm, teff, x, lteff;
     double ub, bv, bcv, vico, rico, jho, jhco, hko, hkco, klo, vko, v12;
     double obs[16], repro[16], stellar[16], extinc[16], flux[16], obsSED[16], derobsSED[16], reproSED[16], starSED[16], lambda[16];
     double av, mbolv, llumv, Mv, ebv, bvo, mvo;
     int i;
     char name[7], type[5], scrchar[20];

/* Set constants   */

/* Units of Watts cm^-2 mu^-1 for a zero-magnitude star */

	flux[0] = 4.19e-12; flux[1] = 6.59e-12; flux[2] = 3.60e-12; flux[3] = 2.25e-12;
        flux[4] = 1.23e-12; flux[5] = 3.03e-13; flux[6] = 1.26e-13; flux[7]=4.06e-14;
	flux[8] = 6.89e-15; flux[9] = 2.21e-15; flux[10] = 9.60e-17; flux[11] = 5.90e-17;
	flux[12] = 6.40e-18; flux[13] = 3.23e-18; flux[14] = 9.92e-20; flux[15] = 1.29e-20;

/* Wavelength scale adopted U, B, V, Rc, Ic, J, H, K, L, M, N, IRAS12, Q,
                            IRAS25, IRAS60, IRAS100                            */

	lambda[0] = 0.36; lambda[1] = 0.44; lambda[2] = 0.55; lambda[3] = 0.64;
	lambda[4] = 0.79; lambda[5] = 1.26; lambda[6] = 1.60; lambda[7] = 2.22;
	lambda[8] = 3.54; lambda[9] = 4.80; lambda[10] = 10.6; lambda[11] = 12.0;
	lambda[12] = 21.0; lambda[13] = 25.0; lambda[14] = 60.0; lambda[15] = 100.0;

/*  Extinction law adopted as Alam/Av from Reike and Lebofsky UBV, Rc, Ic, JHKLMN12  */

	extinc[0] = 1.53; extinc[1] = 1.32; extinc[2] = 1.0; extinc[3] = 0.82,
	extinc[4] = 0.60; extinc[5] = 0.265; extinc[6] = 0.175; extinc[7] = 0.112;
	extinc[8] = 0.058; extinc[9] = 0.02; extinc[10] = 0.05; extinc[11] = 0.04;
	extinc[12] = 0.02; extinc[13] = 0.01; extinc[14] = 0.00; extinc[15] = 0.00;

/*  Read in sources, dm, ebv, U, B, V, Rc, Ic, Jc, Hc, Kc, L, M, N, 12, Q, 25, 60, 100 */

while (scanf("%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
name, &dm, &ebv, &obs[0], &obs[1], &obs[2], &obs[3], &obs[4], &obs[5], &obs[6], &obs[7], &obs[8],
&obs[9], &obs[10], &obs[11], &obs[12], &obs[13], &obs[14], &obs[15]) !=EOF)  {

        bv = obs[1] - obs[2];
        bvo = bv - ebv;
	if ( 0.00 < bvo < 1.48 ) {
        teff =   9556.859 -11029.78*bvo + 12990.49*bvo*bvo -9241.79*bvo*bvo*bvo + 2502.884*bvo*bvo*bvo*bvo;
	}
	x = teff;
	lteff = log(teff)*0.4343;

/*  Calculate Extinction   */

	av = (ebv)/0.32;

/*  Use the derived fits (references in header). */

	ub = 2.674702E1 - x*3.675133E-2 + x*x*1.996307E-5 - x*x*x*5.315118E-9 + x*x*x*x*7.388618E-13 - x*x*x*x*x*5.155066E-17 + x*x*x*x*x*x*1.426480E-21;

	bv = 9.026432E1 - 1.116759E-1*x + x*x*5.917468E-5 -  x*x*x*1.703182E-8 + x*x*x*x*2.867242E-12 - x*x*x*x*x*2.826478E-16 + x*x*x*x*x*x*1.513376E-20 - x*x*x*x*x*x*x*3.400976E-25;

        bcv = -101.3651 + x*8.971698E-2 + x*x*-3.290524E-5 + x*x*x*6.334656E-9 + x*x*x*x*-6.714475E-13 + x*x*x*x*x*3.705377E-17 + x*x*x*x*x*x*-8.309852E-22;

	rico = 2.195985E1 - 1.499330E-2*x + x*x*4.177120E-6 -  x*x*x*5.805917E-10 + x*x*x*x*3.992331E-14 - x*x*x*x*x*1.084173E-18;

	vico = 1.103107e2 - 1.147069e-1*x + x*x*5.191987e-5 -  x*x*x*1.296276e-8 + x*x*x*x*1.911211E-12 - x*x*x*x*x*1.657877E-16 + x*x*x*x*x*x*7.820634E-21 - x*x*x*x*x*x*x*1.546372E-25;

	vko = 231.6392 - 0.2627116*x + x*x*1.32161e-4 - x*x*x*3.773421e-8 + x*x*x*x*6.658074e-12 - x*x*x*x*x*7.432943e-16 + x*x*x*x*x*x*5.129903e-20 - x*x*x*x*x*x*x*2.001820e-24 + x*x*x*x*x*x*x*x*3.381241e-29;

	jho = 2.857311E2 + x*-4.551107E-1 + x*x*3.134182E-4 + x*x*x*-1.224157E-7 + x*x*x*x*2.994256E-11 + x*x*x*x*x*-4.764616E-15 + x*x*x*x*x*x*4.938796E-19 + x*x*x*x*x*x*x*-3.218916E-23 + x*x*x*x*x*x*x*x*1.197969E-27 + x*x*x*x*x*x*x*x*x*-1.940957E-32;

	hko = 2.788605E0 + x*-1.731798E-3 + x*x*4.424067E-7 + x*x*x*-5.647841E-11 + x*x*x*x*3.559478E-15 + x*x*x*x*x*-8.828952E-20;

	klo = 4.964018E0 - x*3.818125E-3 + x*x*1.191141E-6 - x*x*x*1.846119E-10 + x*x*x*x*1.411467E-14 - x*x*x*x*x*4.248989E-19;

	v12 =  1.257185E1 - x*3.235974E-3 + x*x*2.866275E-7 - x*x*x*8.819869E-12;

/*  Perform color transformations	*/

	jhco =  jho*0.911; hkco =  hko*0.971 - 0.02;


/*  Calculate absolute and dereddened Magnitudes       */

        mvo = obs[2] - av;
	Mv = obs[2] - av - dm;

/*  Calculate luminosities */

	mbolv = Mv + bcv; llumv = 1.89 - 0.4*mbolv;

/* Calculate stellar SED normalized at V-band  - NOTE DISTANCE NOT USED! - */

	stellar[0] = ub + bv + mvo; stellar[1] = bv + mvo;
	stellar[2] = mvo; stellar[3] = rico + mvo - vico;
	stellar[4] = mvo - vico;
	stellar[5] = stellar[2] - vko + hkco + jhco; stellar[6] = stellar[5] - jhco;
	stellar[7] = stellar[6] - hkco; stellar[8] = stellar[7] - klo;
	stellar[9] = stellar[8]; stellar[10] = stellar[9]; stellar[11] =  stellar[10];
	stellar[12] = stellar[11];  stellar[13] = stellar[12]; stellar[14] = stellar[13];
	stellar[15] = stellar[14];

/* Calculate star + face-on reprocessing disk from Hillenbrand et al. (1992) APJ V397 613  */

	if (12000.0 > teff && teff > 8600.0) {

        repro[0] = stellar[0]; repro[1] = stellar[1]; repro[2] = stellar[2];
        repro[3] = stellar[3] - 0.49; repro[4] = stellar[4] - 0.72;
        repro[5] = stellar[5] - 1.09; repro[6] = stellar[6] - 1.41;
        repro[7] = stellar[7] - 1.86; repro[8] = stellar[8] - 2.60;
        repro[9] = stellar[9] - 3.09; repro[10] = stellar[10] - 4.45;
        repro[11] = stellar[11] - 4.65; repro[12] = stellar[12] - 5.51;
        repro[13] = stellar[13] - 5.95; repro[14] = stellar[14] - 7.51;
        repro[15] = stellar[15] - 8.43;
	}

	if (8600.0 > teff && teff > 6400.0) {

        repro[0] = stellar[0]; repro[1] = stellar[1]; repro[2] = stellar[2];
        repro[3] = stellar[3] - 0.26; repro[4] = stellar[4] - 0.42;
        repro[5] = stellar[5] - 0.72; repro[6] = stellar[6] - 0.99;
        repro[7] = stellar[7] - 1.39; repro[8] = stellar[8] - 2.07;
        repro[9] = stellar[9] - 2.54; repro[10] = stellar[10] - 3.87;
        repro[11] = stellar[11] - 4.07; repro[12] = stellar[12] - 4.91;
        repro[13] = stellar[13] - 5.35; repro[14] = stellar[14] - 6.91;
        repro[15] = stellar[15] - 7.83;
	}

	if (6400.0 > teff && teff > 5800.0) {

        repro[0] = stellar[0]; repro[1] = stellar[1]; repro[2] = stellar[2];
        repro[3] = stellar[3] - 0.17; repro[4] = stellar[4] - 0.30;
        repro[5] = stellar[5] - 0.55; repro[6] = stellar[6] - 0.79;
        repro[7] = stellar[7] - 1.16; repro[8] = stellar[8] - 1.81;
        repro[9] = stellar[9] - 2.27; repro[10] = stellar[10] - 3.57;
        repro[11] = stellar[11] - 3.77; repro[12] = stellar[12] - 4.61;
        repro[13] = stellar[13] - 5.04; repro[14] = stellar[14] - 6.60;
        repro[15] = stellar[15] - 7.51;
	}

	if (5800.0 > teff && teff > 5000.0) {

        repro[0] = stellar[0]; repro[1] = stellar[1]; repro[2] = stellar[2];
        repro[3] = stellar[3] - 0.13; repro[4] = stellar[4] - 0.23;
        repro[5] = stellar[5] - 0.45; repro[6] = stellar[6] - 0.67;
        repro[7] = stellar[7] - 1.01; repro[8] = stellar[8] - 1.64;
        repro[9] = stellar[9] - 2.08; repro[10] = stellar[10] - 3.37;
        repro[11] = stellar[11] - 3.57; repro[12] = stellar[12] - 4.40;
        repro[13] = stellar[13] - 4.83; repro[14] = stellar[14] - 6.38;
        repro[15] = stellar[15] - 7.30;
	}

	if (5000.0 > teff && teff > 4000.0) {

	repro[0] = stellar[0]; repro[1] = stellar[1]; repro[2] = stellar[2];
	repro[3] = stellar[3] - 0.08; repro[4] = stellar[4] - 0.15;
	repro[5] = stellar[5] - 0.32; repro[6] = stellar[6] - 0.51;
	repro[7] = stellar[7] - 0.82; repro[8] = stellar[8] - 1.40;
	repro[9] = stellar[9] - 1.82; repro[10] = stellar[10] - 3.08;
	repro[11] = stellar[11] - 3.27; repro[12] = stellar[12] - 4.10;
	repro[13] = stellar[13] - 4.53; repro[14] = stellar[14] - 6.07;
	repro[15] = stellar[15] - 6.98;
	}

	if (4000.0 > teff && teff > 3600.0) {

        repro[0] = stellar[0]; repro[1] = stellar[1]; repro[2] = stellar[2];
        repro[3] = stellar[3] - 0.05; repro[4] = stellar[4] - 0.11;
        repro[5] = stellar[5] - 0.26; repro[6] = stellar[6] - 0.42;
        repro[7] = stellar[7] - 0.70; repro[8] = stellar[8] - 1.25;
        repro[9] = stellar[9] - 1.66; repro[10] = stellar[10] - 2.89;
        repro[11] = stellar[11] - 3.08; repro[12] = stellar[12] - 3.90;
        repro[13] = stellar[13] - 4.33; repro[14] = stellar[14] - 5.87;
        repro[15] = stellar[15] - 6.78;
	}

	if (3600.0 > teff && teff > 3000.0) {

        repro[0] = stellar[0]; repro[1] = stellar[1]; repro[2] = stellar[2];
        repro[3] = stellar[3] - 0.03; repro[4] = stellar[4] - 0.07;
        repro[5] = stellar[5] - 0.18; repro[6] = stellar[6] - 0.31;
        repro[7] = stellar[7] - 0.55; repro[8] = stellar[8] - 1.06;
        repro[9] = stellar[9] - 1.44; repro[10] = stellar[10] - 2.64;
        repro[11] = stellar[11] - 2.83; repro[12] = stellar[12] - 3.64;
        repro[13] = stellar[13] - 4.06; repro[14] = stellar[14] - 5.59;
        repro[15] = stellar[15] - 6.50;
	}

/*  Calculate the observed and expected SED  and print out the results.  */

	/*  Print out the header */

	printf("### %10s ************************************ \n", name);
	printf("# Teff = %6.0f.  log(L*) = %6.2f.  Av = %4.2f bvo = %4.2f dm = %4.2f \n", teff, llumv, av, bvo, dm);
	printf("###************************************ \n");
	printf("#Lambda      Fobs     deredF  F*      Frepro  \n");
        printf("#=========   ======   ======  =====   ======  \n");


	for (i = 13; i < 15; i++) {

			obsSED[i] = flux[i]/pow(10.0,obs[i]/2.5);
			derobsSED[i] = flux[i]/pow(10.0,(obs[i] - av*extinc[i])/2.5);
			starSED[i] = flux[i]/pow(10.0,stellar[i]/2.5);
			reproSED[i] = flux[i]/pow(10.0,repro[i]/2.5);

			printf("%6.2f %10.3e %10.3e %10.3e %10.3e \n", lambda[i], lambda[i]*obsSED[i], lambda[i]*derobsSED[i], lambda[i]*starSED[i], lambda[i]*reproSED[i]);
	}
     }
}
