*  Copyright 2010 Pim Schellart. All rights reserved.
*  
*  Redistribution and use in source and binary forms, with or
*  without modification, are permitted provided that the following
*  conditions are met:
*  
*     1. Redistributions of source code must retain the above
*        copyright notice, this list of conditions and the following
*        disclaimer.
*  
*     2. Redistributions in binary form must reproduce the above
*        copyright notice, this list of conditions and the following
*        disclaimer in the documentation and/or other materials
*        provided with the distribution.
*  
*  THIS SOFTWARE IS PROVIDED BY PIM SCHELLART ``AS IS'' AND ANY
*  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
*  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
*  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PIM SCHELLART OR
*  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
*  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
*  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
*  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
*  OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
*  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
*  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*  
*  The views and conclusions contained in the software and documentation
*  are those of the authors and should not be interpreted as representing
*  official policies, either expressed or implied, of Pim Schellart.

      SUBROUTINE LOMBSCARGLE(T, X, W, P, NT, NW)
      IMPLICIT NONE

Cf2py intent(in) T(NT)
Cf2py intent(in) X(NT)
Cf2py intent(in) W(NW)
Cf2py intent(out) P(NW)
Cf2py intent(hide) NT
Cf2py intent(hide) NW

* 
*  Input arguments
* 
      INTEGER NT, NW
      DOUBLE PRECISION T(NT), X(NT), W(NW)
      DOUBLE PRECISION P(NW)

* 
*  Purpose
*  =======
* 
*  Computes the Lomb-Scargle periodogram as developed by Lomb (1976)
*  and further extended by Scargle (1982) to find, and test the
*  significance of weak periodic signals with uneven temporal sampling.
* 
*  This subroutine calculates the periodogram using a slightly
*  modified algorithm due to Townsend (2010) which allows the
*  periodogram to be calculated using only a single pass through
*  the input samples.
*  This requires Nw(2Nt+3) trigonometric function evaluations (where
*  Nw is the number of frequencies and Nt the number of input samples),
*  giving a factor of ~2 speed increase over the straightforward
*  implementation.
* 
*  Arguments
*  =========
* 
*  T   (input) double precision array, dimension (NT)
*      Sample times
* 
*  X   (input) double precision array, dimension (NT)
*      Measurement values
* 
*  W   (input) double precision array, dimension (NT)
*      Angular frequencies for output periodogram
* 
*  P   (output) double precision array, dimension (NW)
*      Lomb-Scargle periodogram
* 
*  NT (input) integer
*      Dimension of input arrays
* 
*  NW (output) integer
*      Dimension of output array
* 
*  Further details
*  ===============
* 
*  P(i) takes a value of A^2*N/4 for a harmonic signal with
*  frequency w(i).
* 

* 
*  Local variables
* 
      INTEGER I, J
      DOUBLE PRECISION C, S, XC, XS, CC, SS, CS
      DOUBLE PRECISION TAU, C_TAU, S_TAU, C_TAU2, S_TAU2, CS_TAU

      DO I = 1, NW

          XC = 0.
          XS = 0.
          CC = 0.
          SS = 0.
          CS = 0.

          DO J = 1, NT

              C = COS(W(I) * T(J))
              S = SIN(W(I) * T(J))

              XC = XC + X(J) * C
              XS = XS + X(J) * S
              CC = CC + C * C
              SS = SS + S * S
              CS = CS + C * S

          END DO

          TAU = ATAN(2 * CS / (CC - SS)) / (2 * W(I))
          C_TAU = COS(W(I) * TAU)
          S_TAU = SIN(W(I) * TAU)
          C_TAU2 = C_TAU * C_TAU
          S_TAU2 = S_TAU * S_TAU
          CS_TAU = 2 * C_TAU * S_TAU

          P(i) = 0.5 * (
     $          ((C_TAU * XC + S_TAU * XS)**2 /
     $          (C_TAU2 * CC + CS_TAU * CS + S_TAU2 * SS)) +
     $          ((C_TAU * XS - S_TAU * XC)**2 /
     $          (C_TAU2 * SS - CS_TAU * CS + S_TAU2 * CC)))

      END DO

      END SUBROUTINE LOMBSCARGLE