/*
 * Copyright (c) 2003, 2007-11 Matteo Frigo
 * Copyright (c) 2003, 2007-11 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* This file was automatically generated --- DO NOT EDIT */
/* Generated on Wed Jul 27 06:12:26 EDT 2011 */

#include "codelet-dft.h"

#ifdef HAVE_FMA

/* Generated by: ../../../genfft/gen_notw.native -fma -reorder-insns -schedule-for-pipeline -compact -variables 4 -pipeline-latency 4 -n 5 -name n1_5 -include n.h */

/*
 * This function contains 32 FP additions, 18 FP multiplications,
 * (or, 14 additions, 0 multiplications, 18 fused multiply/add),
 * 37 stack variables, 4 constants, and 20 memory accesses
 */
#include "n.h"

static void n1_5(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DK(KP951056516, +0.951056516295153572116439333379382143405698634);
     DK(KP559016994, +0.559016994374947424102293417182819058860154590);
     DK(KP250000000, +0.250000000000000000000000000000000000000000000);
     DK(KP618033988, +0.618033988749894848204586834365638117720309180);
     {
	  INT i;
	  for (i = v; i > 0; i = i - 1, ri = ri + ivs, ii = ii + ivs, ro = ro + ovs, io = io + ovs, MAKE_VOLATILE_STRIDE(is), MAKE_VOLATILE_STRIDE(os)) {
	       E Tq, Ti, Tk, Tu, Tw, Tp, Tb, Tj, Tr, Tv;
	       {
		    E T1, Tl, Ts, Tt, T8, Ta, Te, Tm, Tn, Th, To, T9;
		    T1 = ri[0];
		    Tl = ii[0];
		    {
			 E T2, T3, T5, T6;
			 T2 = ri[WS(is, 1)];
			 T3 = ri[WS(is, 4)];
			 T5 = ri[WS(is, 2)];
			 T6 = ri[WS(is, 3)];
			 {
			      E Tc, T4, T7, Td, Tf, Tg;
			      Tc = ii[WS(is, 1)];
			      Ts = T2 - T3;
			      T4 = T2 + T3;
			      Tt = T5 - T6;
			      T7 = T5 + T6;
			      Td = ii[WS(is, 4)];
			      Tf = ii[WS(is, 2)];
			      Tg = ii[WS(is, 3)];
			      T8 = T4 + T7;
			      Ta = T4 - T7;
			      Te = Tc - Td;
			      Tm = Tc + Td;
			      Tn = Tf + Tg;
			      Th = Tf - Tg;
			 }
		    }
		    ro[0] = T1 + T8;
		    To = Tm + Tn;
		    Tq = Tm - Tn;
		    Ti = FMA(KP618033988, Th, Te);
		    Tk = FNMS(KP618033988, Te, Th);
		    io[0] = Tl + To;
		    T9 = FNMS(KP250000000, T8, T1);
		    Tu = FMA(KP618033988, Tt, Ts);
		    Tw = FNMS(KP618033988, Ts, Tt);
		    Tp = FNMS(KP250000000, To, Tl);
		    Tb = FMA(KP559016994, Ta, T9);
		    Tj = FNMS(KP559016994, Ta, T9);
	       }
	       Tr = FMA(KP559016994, Tq, Tp);
	       Tv = FNMS(KP559016994, Tq, Tp);
	       ro[WS(os, 2)] = FNMS(KP951056516, Tk, Tj);
	       ro[WS(os, 3)] = FMA(KP951056516, Tk, Tj);
	       ro[WS(os, 1)] = FMA(KP951056516, Ti, Tb);
	       ro[WS(os, 4)] = FNMS(KP951056516, Ti, Tb);
	       io[WS(os, 2)] = FMA(KP951056516, Tw, Tv);
	       io[WS(os, 3)] = FNMS(KP951056516, Tw, Tv);
	       io[WS(os, 4)] = FMA(KP951056516, Tu, Tr);
	       io[WS(os, 1)] = FNMS(KP951056516, Tu, Tr);
	  }
     }
}

static const kdft_desc desc = { 5, "n1_5", {14, 0, 18, 0}, &GENUS, 0, 0, 0, 0 };

void X(codelet_n1_5) (planner *p) {
     X(kdft_register) (p, n1_5, &desc);
}

#else				/* HAVE_FMA */

/* Generated by: ../../../genfft/gen_notw.native -compact -variables 4 -pipeline-latency 4 -n 5 -name n1_5 -include n.h */

/*
 * This function contains 32 FP additions, 12 FP multiplications,
 * (or, 26 additions, 6 multiplications, 6 fused multiply/add),
 * 21 stack variables, 4 constants, and 20 memory accesses
 */
#include "n.h"

static void n1_5(const R *ri, const R *ii, R *ro, R *io, stride is, stride os, INT v, INT ivs, INT ovs)
{
     DK(KP250000000, +0.250000000000000000000000000000000000000000000);
     DK(KP587785252, +0.587785252292473129168705954639072768597652438);
     DK(KP951056516, +0.951056516295153572116439333379382143405698634);
     DK(KP559016994, +0.559016994374947424102293417182819058860154590);
     {
	  INT i;
	  for (i = v; i > 0; i = i - 1, ri = ri + ivs, ii = ii + ivs, ro = ro + ovs, io = io + ovs, MAKE_VOLATILE_STRIDE(is), MAKE_VOLATILE_STRIDE(os)) {
	       E T1, To, T8, Tt, T9, Ts, Te, Tp, Th, Tn;
	       T1 = ri[0];
	       To = ii[0];
	       {
		    E T2, T3, T4, T5, T6, T7;
		    T2 = ri[WS(is, 1)];
		    T3 = ri[WS(is, 4)];
		    T4 = T2 + T3;
		    T5 = ri[WS(is, 2)];
		    T6 = ri[WS(is, 3)];
		    T7 = T5 + T6;
		    T8 = T4 + T7;
		    Tt = T5 - T6;
		    T9 = KP559016994 * (T4 - T7);
		    Ts = T2 - T3;
	       }
	       {
		    E Tc, Td, Tl, Tf, Tg, Tm;
		    Tc = ii[WS(is, 1)];
		    Td = ii[WS(is, 4)];
		    Tl = Tc + Td;
		    Tf = ii[WS(is, 2)];
		    Tg = ii[WS(is, 3)];
		    Tm = Tf + Tg;
		    Te = Tc - Td;
		    Tp = Tl + Tm;
		    Th = Tf - Tg;
		    Tn = KP559016994 * (Tl - Tm);
	       }
	       ro[0] = T1 + T8;
	       io[0] = To + Tp;
	       {
		    E Ti, Tk, Tb, Tj, Ta;
		    Ti = FMA(KP951056516, Te, KP587785252 * Th);
		    Tk = FNMS(KP587785252, Te, KP951056516 * Th);
		    Ta = FNMS(KP250000000, T8, T1);
		    Tb = T9 + Ta;
		    Tj = Ta - T9;
		    ro[WS(os, 4)] = Tb - Ti;
		    ro[WS(os, 3)] = Tj + Tk;
		    ro[WS(os, 1)] = Tb + Ti;
		    ro[WS(os, 2)] = Tj - Tk;
	       }
	       {
		    E Tu, Tv, Tr, Tw, Tq;
		    Tu = FMA(KP951056516, Ts, KP587785252 * Tt);
		    Tv = FNMS(KP587785252, Ts, KP951056516 * Tt);
		    Tq = FNMS(KP250000000, Tp, To);
		    Tr = Tn + Tq;
		    Tw = Tq - Tn;
		    io[WS(os, 1)] = Tr - Tu;
		    io[WS(os, 3)] = Tw - Tv;
		    io[WS(os, 4)] = Tu + Tr;
		    io[WS(os, 2)] = Tv + Tw;
	       }
	  }
     }
}

static const kdft_desc desc = { 5, "n1_5", {26, 6, 6, 0}, &GENUS, 0, 0, 0, 0 };

void X(codelet_n1_5) (planner *p) {
     X(kdft_register) (p, n1_5, &desc);
}

#endif				/* HAVE_FMA */
