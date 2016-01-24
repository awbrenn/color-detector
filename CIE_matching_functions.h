//
// Created by awbrenn on 1/24/16.
//

#ifndef COLOR_DETECTOR_CIE_MATCHING_FUNCTIONS_H_H
#define COLOR_DETECTOR_CIE_MATCHING_FUNCTIONS_H_H


struct point {
  double x,y;
};

struct point fx[]={
380.0, 0.0013928,
390.0, 0.0042563,
400.0, 0.0144416,
410.0, 0.0431500,
420.0, 0.1344313,
430.0, 0.2839478,
440.0, 0.3468990,
450.0, 0.3362033,
460.0, 0.2909090,
470.0, 0.1953823,
480.0, 0.0956345,
490.0, 0.0320108,
500.0, 0.0049194,
510.0, 0.0093197,
520.0, 0.0632681,
530.0, 0.1654689,
540.0, 0.2903863,
550.0, 0.4334710,
560.0, 0.5944697,
570.0, 0.7621835,
580.0, 0.9162317,
590.0, 1.0265858,
600.0, 1.0619702,
610.0, 1.0028369,
620.0, 0.8544512,
630.0, 0.6424434,
640.0, 0.4478554,
650.0, 0.2835212,
660.0, 0.1649259,
670.0, 0.0874030,
680.0, 0.0467659,
690.0, 0.0226734,
700.0, 0.0113525,
710.0, 0.0058147,
720.0, 0.0029073,
730.0, 0.0014398,
740.0, 0.0006922,
750.0, 0.0003323,
760.0, 0.0001661,
770.0, 0.0000831
};

struct point fy[]={
380.0, 0.00004,
390.0, 0.00012,
400.0, 0.0004,
410.0, 0.0012,
420.0, 0.004,
430.0, 0.0116,
440.0, 0.023,
450.0, 0.038,
460.0, 0.06,
470.0, 0.091,
480.0, 0.139,
490.0, 0.208,
500.0, 0.323,
510.0, 0.503,
520.0, 0.71,
530.0, 0.862,
540.0, 0.954,
550.0, 0.995,
560.0, 0.995,
570.0, 0.952,
580.0, 0.87,
590.0, 0.757,
600.0, 0.631,
610.0, 0.503,
620.0, 0.381,
630.0, 0.265,
640.0, 0.175,
650.0, 0.107,
660.0, 0.061,
670.0, 0.032,
680.0, 0.017,
690.0, 0.0082,
700.0, 0.0041,
710.0, 0.0021,
720.0, 0.00105,
730.0, 0.00052,
740.0, 0.00025,
750.0, 0.00012,
760.0, 0.00006,
770.0, 0.00003
};

struct point fz[]={
380.0, 0.0065672,
390.0, 0.0201134,
400.0, 0.0684916,
410.0, 0.2056500,
420.0, 0.6458823,
430.0, 1.3856115,
440.0, 1.7401926,
450.0, 1.7726892,
460.0, 1.6692929,
470.0, 1.2880121,
480.0, 0.8128409,
490.0, 0.4650738,
500.0, 0.2720063,
510.0, 0.1581680,
520.0, 0.0782549,
530.0, 0.0421426,
540.0, 0.0203624,
550.0, 0.0087671,
560.0, 0.0038558,
570.0, 0.0020595,
580.0, 0.0015374,
590.0, 0.0011600,
600.0, 0.0007621,
610.0, 0.0003765,
620.0, 0.0001977,
630.0, 0.0000907,
640.0, 0.0000311,
650.0, 0.0000039,
660.0, 0.0,
670.0, 0.0,
680.0, 0.0,
690.0, 0.0,
700.0, 0.0,
710.0, 0.0,
720.0, 0.0,
730.0, 0.0,
740.0, 0.0,
750.0, 0.0,
760.0, 0.0,
770.0, 0.0
};

#endif //COLOR_DETECTOR_CIE_MATCHING_FUNCTIONS_H_H
